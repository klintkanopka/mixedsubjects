#' @title DiP++ Estimator
#' @description Power-tuned Difference-in-Predictions estimator for ATE.

#' DiP++ Estimator
#'
#' Computes the DiP++ (power-tuned difference-in-predictions) estimator for
#' the average treatment effect (ATE). This estimator uses paired predictions
#' S^(1) and S^(0) for each unlabeled unit, with a single tuning parameter
#' lambda estimated via cross-fitting.
#'
#' @param formula_or_data Either an msd_data object created by \code{\link{msd_data}},
#'   or a formula of the form \code{outcome ~ treatment | pred_treated + pred_control}.
#' @param data If \code{formula_or_data} is a formula, this should be either:
#'   an msd_data object, a combined dataframe, or NULL (if using observed/unobserved).
#' @param observed If using formula with separate dataframes, the observed data.
#' @param unobserved If using formula with separate dataframes, the unobserved data.
#' @param n_folds Number of folds for cross-fitting (default 2)
#' @param conf_level Confidence level for the confidence interval (default 0.95)
#' @param seed Random seed for fold splitting (optional)
#'
#' @return An msd_result object containing:
#'   \item{estimate}{Point estimate of the ATE}
#'   \item{variance}{Estimated variance (delta-method)}
#'   \item{se}{Standard error}
#'   \item{ci_lower, ci_upper}{Confidence interval bounds}
#'   \item{method}{Name of the estimation method}
#'   \item{lambda}{Estimated tuning parameter (single value)}
#'
#' @details
#' The DiP++ estimator is:
#' \deqn{\hat{\tau}^{DiP++}(\lambda) = \frac{\lambda}{|\mathcal{U}|}
#'   \sum_{i \in \mathcal{U}} (S_i^{(1)} - S_i^{(0)}) +
#'   \frac{1}{n_1}\sum_{i \in \mathcal{O}_1}(Y_i - \lambda S_i^{(1)}) -
#'   \frac{1}{n_0}\sum_{i \in \mathcal{O}_0}(Y_i - \lambda S_i^{(0)})}
#'
#' @note
#' DiP++ requires BOTH S0 and S1 predictions for ALL units.
#' For arm-specific tuning, use \code{\link{msd_dt_dip}}.
#'
#' @examples
#' # Using msd_data object
#' set.seed(123)
#' n <- 100
#' obs_df <- data.frame(
#'   Y = rnorm(n),
#'   D = rep(c(1, 0), each = n/2)
#' )
#' obs_df$Y <- obs_df$Y + 0.3 * obs_df$D
#' obs_df$S1 <- 0.5 * obs_df$Y + rnorm(n, 0, 0.5)
#' obs_df$S0 <- 0.5 * obs_df$Y + rnorm(n, 0, 0.5) - 0.1
#'
#' unobs_df <- data.frame(
#'   S0 = rnorm(200, 0, 0.5),
#'   S1 = rnorm(200, 0.2, 0.5),
#'   D = rep(c(1, 0), each = 100)
#' )
#'
#' msd <- msd_data(observed = obs_df, unobserved = unobs_df)
#' result <- msd_dip_pp(msd)
#'
#' # Using formula interface
#' result2 <- msd_dip_pp(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
#'
#' @export
msd_dip_pp <- function(formula_or_data, data = NULL, observed = NULL,
                       unobserved = NULL, n_folds = 2, conf_level = 0.95,
                       seed = NULL) {

  # Resolve flexible input to msd_data
  msd <- resolve_msd_data(formula_or_data, data, observed, unobserved)

  if (!msd$has_both_predictions) {
    stop("DiP++ requires both S0 and S1 predictions for all units.")
  }

  if (is.null(msd$unobserved) || msd$m == 0) {
    stop("DiP++ requires unlabeled data. Use msd_dim() for labeled-only estimation.")
  }

  # Extract data
  obs <- msd$observed
  unobs <- msd$unobserved

  # Get sample sizes
  n1 <- sum(obs$D == 1)
  n0 <- sum(obs$D == 0)
  m <- nrow(unobs)

  # Create fold assignments for each arm
  fold_info <- split_by_arm_fold(msd, n_folds, seed)

  # Storage for fold estimates and lambdas
  fold_estimates <- numeric(n_folds)
  fold_lambdas <- numeric(n_folds)

  # Cross-fitting loop
  for (k in 1:n_folds) {
    # Get indices for this fold and opposite folds
    treated_in_k <- fold_info$treated_idx[fold_info$treated_folds == k]
    treated_not_k <- fold_info$treated_idx[fold_info$treated_folds != k]
    control_in_k <- fold_info$control_idx[fold_info$control_folds == k]
    control_not_k <- fold_info$control_idx[fold_info$control_folds != k]

    # Estimate lambda on opposite folds
    lambda_k <- estimate_lambda_dip_pp(obs, treated_not_k, control_not_k,
                                        unobs, m)

    # Compute estimate on fold k using lambda from opposite folds
    est_k <- compute_dip_pp_estimate(obs, unobs, treated_in_k, control_in_k,
                                      lambda_k)

    fold_estimates[k] <- est_k
    fold_lambdas[k] <- lambda_k
  }

  # Average across folds (equal weights)
  estimate <- mean(fold_estimates)
  lambda <- mean(fold_lambdas)

  # Compute variance using delta method
  variance <- compute_dip_pp_variance(msd, lambda, n_folds)

  # Standard error
  se <- sqrt(variance)

  # Confidence interval
  ci <- compute_ci(estimate, se, conf_level)

  # Return result
  new_msd_result(
    estimate = estimate,
    variance = variance,
    se = se,
    ci_lower = ci["ci_lower"],
    ci_upper = ci["ci_upper"],
    method = paste0("DiP++ (cross-fit, K=", n_folds, ")"),
    lambda = lambda,
    n1 = n1,
    n0 = n0,
    m = m,
    conf_level = conf_level,
    additional = list(
      fold_lambdas = fold_lambdas,
      fold_estimates = fold_estimates
    )
  )
}

#' Estimate lambda for DiP++ (single lambda)
#' @keywords internal
estimate_lambda_dip_pp <- function(obs, treated_idx, control_idx, unobs, m) {
  # Get data
  Y1 <- obs$Y[treated_idx]
  Y0 <- obs$Y[control_idx]
  S1_obs_1 <- obs$S1[treated_idx]
  S0_obs_0 <- obs$S0[control_idx]

  n1 <- length(Y1)
  n0 <- length(Y0)

  var_S_diff <- var(unobs$S1 - unobs$S0)
  var_S1_obs <- var(S1_obs_1)
  var_S0_obs <- var(S0_obs_0)
  cov_YS_1 <- cov(Y1, S1_obs_1)
  cov_YS_0 <- cov(Y0, S0_obs_0)

  C <- var_S_diff / m + var_S1_obs / n1 + var_S0_obs / n0
  D <- cov_YS_1 / n1 + cov_YS_0 / n0

  if (C > 0) {
    lambda <- D / C
  } else {
    lambda <- 0
  }

  return(lambda)
}

#' Compute DiP++ estimate for a fold
#' @keywords internal
compute_dip_pp_estimate <- function(obs, unobs, treated_idx, control_idx, lambda) {
  Y1 <- obs$Y[treated_idx]
  Y0 <- obs$Y[control_idx]
  S1_obs_1 <- obs$S1[treated_idx]
  S0_obs_0 <- obs$S0[control_idx]

  mean_S_diff_unobs <- mean(unobs$S1 - unobs$S0)

  estimate <- lambda * mean_S_diff_unobs +
              mean(Y1 - lambda * S1_obs_1) -
              mean(Y0 - lambda * S0_obs_0)

  return(estimate)
}

#' Compute DiP++ variance using delta method
#'
#' The unobserved component var_U = lambda^2 * Var(S1-S0) / m is shared across
#' all folds and is NOT divided by K. The labeled components' K factors cancel
#' when averaging across folds.
#' @keywords internal
compute_dip_pp_variance <- function(data, lambda, n_folds) {
  obs <- data$observed
  unobs <- data$unobserved

  n1 <- sum(obs$D == 1)
  n0 <- sum(obs$D == 0)
  m <- nrow(unobs)

  Y1 <- obs$Y[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]
  S1_obs_1 <- obs$S1[obs$D == 1]
  S0_obs_0 <- obs$S0[obs$D == 0]

  # Unobserved component (shared across folds, not divided by K)
  var_U <- lambda^2 * var(unobs$S1 - unobs$S0) / m

  # Labeled components per arm
  var_1_labeled <- (var(Y1) + lambda^2 * var(S1_obs_1) - 2 * lambda * cov(Y1, S1_obs_1)) / n1
  var_0_labeled <- (var(Y0) + lambda^2 * var(S0_obs_0) - 2 * lambda * cov(Y0, S0_obs_0)) / n0

  # Total: unobserved (no /K) + labeled
  variance <- var_U + var_1_labeled + var_0_labeled

  return(max(variance, 0))
}

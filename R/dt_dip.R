#' @title D-T DiP (Doubly-Tuned Difference-in-Predictions) Estimator
#' @description D-T DiP estimator with arm-specific tuning for ATE.

#' D-T DiP (Doubly-Tuned Difference-in-Predictions) Estimator
#'
#' Computes the D-T DiP estimator for the average treatment effect (ATE).
#' This estimator uses paired predictions S^(1) and S^(0) for each unlabeled
#' unit, with arm-specific tuning parameters (lambda_1 and lambda_0) estimated
#' via cross-fitting.
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
#'   \item{lambda}{Vector of arm-specific tuning parameters (lambda_1, lambda_0)}
#'
#' @details
#' The D-T DiP estimator is:
#' \deqn{\hat{\tau}^{D-T DiP} = \frac{1}{|\mathcal{U}|}
#'   \sum_{i \in \mathcal{U}} (\lambda_1 S_i^{(1)} - \lambda_0 S_i^{(0)}) +
#'   \frac{1}{n_1}\sum_{i \in \mathcal{O}_1}(Y_i - \lambda_1 S_i^{(1)}) -
#'   \frac{1}{n_0}\sum_{i \in \mathcal{O}_0}(Y_i - \lambda_0 S_i^{(0)})}
#'
#' Each lambda_d is chosen to minimize the variance in arm d:
#' \deqn{\lambda_d^* = \frac{Cov(Y(d), S^{(d)})}{Var(S^{(d)})}}
#'
#' The tuning parameters are estimated via cross-fitting:
#' \enumerate{
#'   \item Split labeled data into K folds
#'   \item For each fold k, estimate (lambda_1, lambda_0) on opposite folds
#'   \item Compute the fold-k estimate using estimated lambdas
#'   \item Average across folds with equal weights
#' }
#'
#' @note
#' D-T DiP requires BOTH S0 and S1 predictions for ALL units.
#' This corresponds to 2 predictions per unlabeled unit.
#'
#' D-T DiP combines the benefits of:
#' - DiP: exploiting positive correlation between S^(1) and S^(0)
#' - D-T: arm-specific tuning for heterogeneous prediction quality
#'
#' @examples
#' # Create sample data with both predictions
#' set.seed(123)
#' n <- 100
#' obs_df <- data.frame(
#'   Y = rnorm(n),
#'   D = rep(c(1, 0), each = n/2)
#' )
#' obs_df$Y <- obs_df$Y + 0.3 * obs_df$D
#' obs_df$S1 <- 0.6 * obs_df$Y + rnorm(n, 0, 0.4)
#' obs_df$S0 <- 0.4 * obs_df$Y + rnorm(n, 0, 0.5)
#'
#' unobs_df <- data.frame(
#'   S0 = rnorm(300, 0, 0.5),
#'   S1 = rnorm(300, 0.2, 0.4),
#'   D = rep(c(1, 0), 150)
#' )
#' # Add correlation between S0 and S1
#' unobs_df$S1 <- unobs_df$S1 + 0.5 * unobs_df$S0
#'
#' msd <- msd_data(observed = obs_df, unobserved = unobs_df)
#' result <- msd_dt_dip(msd)
#'
#' # Using formula interface
#' result2 <- msd_dt_dip(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
#'
#' @export
msd_dt_dip <- function(formula_or_data, data = NULL, observed = NULL,
                       unobserved = NULL, n_folds = 2, conf_level = 0.95,
                       seed = NULL) {

  # Resolve flexible input to msd_data
  msd <- resolve_msd_data(formula_or_data, data, observed, unobserved)

  if (!msd$has_both_predictions) {
    stop("D-T DiP requires both S0 and S1 predictions for all units.")
  }

  if (is.null(msd$unobserved) || msd$m == 0) {
    stop("D-T DiP requires unlabeled data. Use msd_dim() for labeled-only estimation.")
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
  fold_lambda1 <- numeric(n_folds)
  fold_lambda0 <- numeric(n_folds)

  # Cross-fitting loop
  for (k in 1:n_folds) {
    # Get indices for this fold and opposite folds
    treated_in_k <- fold_info$treated_idx[fold_info$treated_folds == k]
    treated_not_k <- fold_info$treated_idx[fold_info$treated_folds != k]
    control_in_k <- fold_info$control_idx[fold_info$control_folds == k]
    control_not_k <- fold_info$control_idx[fold_info$control_folds != k]

    # Estimate arm-specific lambdas on opposite folds
    lambdas_k <- estimate_lambda_dt_dip(obs, treated_not_k, control_not_k,
                                         unobs, m)

    # Compute estimate on fold k using lambdas from opposite folds
    est_k <- compute_dt_dip_estimate(obs, unobs, treated_in_k, control_in_k,
                                      lambdas_k$lambda1, lambdas_k$lambda0)

    fold_estimates[k] <- est_k
    fold_lambda1[k] <- lambdas_k$lambda1
    fold_lambda0[k] <- lambdas_k$lambda0
  }

  # Average across folds (equal weights)
  estimate <- mean(fold_estimates)
  lambda1 <- mean(fold_lambda1)
  lambda0 <- mean(fold_lambda0)

  # Compute variance using delta method
  variance <- compute_dt_dip_variance(msd, lambda1, lambda0, n_folds)

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
    method = paste0("D-T DiP (cross-fit, K=", n_folds, ")"),
    lambda = c(lambda1, lambda0),
    n1 = n1,
    n0 = n0,
    m = m,
    conf_level = conf_level,
    additional = list(
      fold_lambda1 = fold_lambda1,
      fold_lambda0 = fold_lambda0,
      fold_estimates = fold_estimates,
      cov_S1_S0 = cov(unobs$S1, unobs$S0)
    )
  )
}

#' Estimate arm-specific lambdas for D-T DiP
#' @keywords internal
estimate_lambda_dt_dip <- function(obs, treated_idx, control_idx, unobs, m) {
  # Get data for treatment arm
  Y1 <- obs$Y[treated_idx]
  S1_obs_1 <- obs$S1[treated_idx]
  n1 <- length(Y1)

  # Get data for control arm
  Y0 <- obs$Y[control_idx]
  S0_obs_0 <- obs$S0[control_idx]
  n0 <- length(Y0)

  # For D-T DiP, lambda_d minimizes the arm-specific variance contribution
  # The formulas involve both the labeled and unlabeled variance terms
  #
  # For arm 1: lambda_1^* = cov(Y1, S1) / var(S1) approximately
  # (simplified - full formula includes unlabeled pool contribution)

  var_S1 <- var(S1_obs_1)
  cov_YS_1 <- cov(Y1, S1_obs_1)

  if (var_S1 > 0) {
    # Use regression coefficient as lambda
    lambda1 <- cov_YS_1 / var_S1
  } else {
    lambda1 <- 0
  }

  var_S0 <- var(S0_obs_0)
  cov_YS_0 <- cov(Y0, S0_obs_0)

  if (var_S0 > 0) {
    lambda0 <- cov_YS_0 / var_S0
  } else {
    lambda0 <- 0
  }

  list(lambda1 = lambda1, lambda0 = lambda0)
}

#' Compute D-T DiP estimate for a fold
#' @keywords internal
compute_dt_dip_estimate <- function(obs, unobs, treated_idx, control_idx,
                                     lambda1, lambda0) {
  # Get fold data
  Y1 <- obs$Y[treated_idx]
  Y0 <- obs$Y[control_idx]
  S1_obs_1 <- obs$S1[treated_idx]
  S0_obs_0 <- obs$S0[control_idx]

  # Get unlabeled weighted prediction difference
  # D-T DiP uses: lambda1 * S1 - lambda0 * S0 (not lambda * (S1 - S0))
  mean_weighted_diff <- mean(lambda1 * unobs$S1 - lambda0 * unobs$S0)

  # D-T DiP estimate:
  # tau = mean(lambda1*S1 - lambda0*S0)_U +
  #       mean(Y - lambda1*S1)_O1 - mean(Y - lambda0*S0)_O0
  estimate <- mean_weighted_diff +
              mean(Y1 - lambda1 * S1_obs_1) -
              mean(Y0 - lambda0 * S0_obs_0)

  return(estimate)
}

#' Compute D-T DiP variance using delta method
#'
#' Based on the variance formula from the paper's Supplementary Information.
#' @keywords internal
compute_dt_dip_variance <- function(data, lambda1, lambda0, n_folds) {
  obs <- data$observed
  unobs <- data$unobserved

  # Sample sizes
  n1 <- sum(obs$D == 1)
  n0 <- sum(obs$D == 0)
  m <- nrow(unobs)

  # Per-fold sizes
  n1_k <- n1 / n_folds
  n0_k <- n0 / n_folds

  # Get data
  Y1 <- obs$Y[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]
  S1_obs_1 <- obs$S1[obs$D == 1]
  S0_obs_0 <- obs$S0[obs$D == 0]
  S1_unobs <- unobs$S1
  S0_unobs <- unobs$S0

  # Variance components from observed data
  var_Y1 <- var(Y1)
  var_Y0 <- var(Y0)
  var_S1_obs <- var(S1_obs_1)
  var_S0_obs <- var(S0_obs_0)
  cov_YS_1 <- cov(Y1, S1_obs_1)
  cov_YS_0 <- cov(Y0, S0_obs_0)

  # Variance components from unlabeled data
  var_S1_unobs <- var(S1_unobs)
  var_S0_unobs <- var(S0_unobs)
  cov_S1_S0_unobs <- cov(S1_unobs, S0_unobs)

  # Unlabeled component: Var(lambda1*S1 - lambda0*S0) / m
  var_weighted_diff <- lambda1^2 * var_S1_unobs +
                       lambda0^2 * var_S0_unobs -
                       2 * lambda1 * lambda0 * cov_S1_S0_unobs
  var_U <- var_weighted_diff / m

  # Arm 1 labeled component
  var_1 <- var_Y1 / n1_k +
           lambda1^2 * var_S1_obs / n1_k -
           2 * lambda1 * cov_YS_1 / n1_k

  # Arm 0 labeled component
  var_0 <- var_Y0 / n0_k +
           lambda0^2 * var_S0_obs / n0_k -
           2 * lambda0 * cov_YS_0 / n0_k

  # Cross-fold covariance term (shared unlabeled pool)
  # This is the additional variance from reusing U across folds
  # Approximately: Var(lambda1*S1 - lambda0*S0) / m
  cross_fold_cov <- var_weighted_diff / m

  # Total variance
  # For cross-fit with K=2 folds, variance is approximately:
  # (1/2) * MSE_fold + (1/2) * cross_fold_cov
  variance <- (var_U + var_1 + var_0) / n_folds + cross_fold_cov / (2 * n_folds)

  return(max(variance, 0))  # Ensure non-negative
}

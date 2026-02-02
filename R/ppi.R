#' @title PPI++ Estimator
#' @description Power-tuned Prediction-Powered Inference estimator for ATE.

#' PPI++ Estimator
#'
#' Computes the PPI++ (power-tuned prediction-powered inference) estimator for
#' the average treatment effect (ATE). This estimator uses a single tuning
#' parameter lambda that is estimated via cross-fitting to minimize variance.
#'
#' @param formula_or_data Either an msd_data object created by \code{\link{msd_data}},
#'   or a formula of the form \code{outcome ~ treatment | prediction}.
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
#' The PPI++ estimator uses a single tuning parameter lambda across both arms:
#' \deqn{\hat{\mu}_d^{PPI}(\lambda) = \bar{Y}_{\mathcal{O}_d} +
#'   \lambda(\bar{S}^{(d)}_{\mathcal{U}_d} - \bar{S}^{(d)}_{\mathcal{O}_d})}
#'
#' The tuning parameter is estimated via cross-fitting:
#' \enumerate{
#'   \item Split labeled data into K folds
#'   \item For each fold k, estimate lambda on the opposite folds
#'   \item Compute the fold-k estimate using the estimated lambda
#'   \item Average across folds
#' }
#'
#' Lambda is chosen to minimize the combined variance across arms.
#'
#' @note
#' PPI++ uses a single lambda across arms, unlike D-T which uses arm-specific
#' tuning parameters. For arm-specific tuning, use \code{\link{msd_dt}}.
#'
#' @examples
#' # Create sample data
#' set.seed(123)
#' n <- 100
#' obs_df <- data.frame(
#'   Y = rnorm(n),
#'   S0 = rnorm(n, 0, 0.5),
#'   S1 = rnorm(n, 0.2, 0.5),
#'   D = rep(c(1, 0), each = n/2)
#' )
#' obs_df$Y <- obs_df$Y + 0.3 * obs_df$D
#' obs_df$S1[obs_df$D == 1] <- obs_df$S1[obs_df$D == 1] + 0.5 * obs_df$Y[obs_df$D == 1]
#' obs_df$S0[obs_df$D == 0] <- obs_df$S0[obs_df$D == 0] + 0.5 * obs_df$Y[obs_df$D == 0]
#'
#' unobs_df <- data.frame(
#'   S0 = rnorm(200, 0, 0.5),
#'   S1 = rnorm(200, 0.2, 0.5),
#'   D = rep(c(1, 0), each = 100)
#' )
#'
#' msd <- msd_data(observed = obs_df, unobserved = unobs_df)
#' result <- msd_ppi(msd)
#'
#' # Using formula interface
#' result2 <- msd_ppi(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
#'
#' @export
msd_ppi <- function(formula_or_data, data = NULL, observed = NULL,
                    unobserved = NULL, n_folds = 2, conf_level = 0.95,
                    seed = NULL) {

  # Resolve flexible input to msd_data
  msd <- resolve_msd_data(formula_or_data, data, observed, unobserved)

  if (!msd$has_S0 && !msd$has_S1) {
    stop("PPI++ requires predictions. Neither S0 nor S1 found in data.")
  }

  if (is.null(msd$unobserved) || msd$m == 0) {
    stop("PPI++ requires unlabeled data. Use msd_dim() for labeled-only estimation.")
  }

  # Extract data
  obs <- msd$observed
  unobs <- msd$unobserved

  # Get sample sizes
  n1 <- sum(obs$D == 1)
  n0 <- sum(obs$D == 0)
  m1 <- sum(unobs$D == 1)
  m0 <- sum(unobs$D == 0)

  if (m1 == 0 || m0 == 0) {
    stop("PPI++ requires unlabeled units in both treatment arms.")
  }

  # Create fold assignments for each arm
  fold_info <- split_by_arm_fold(msd, n_folds, seed)

  # Storage for fold estimates
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
    lambda_k <- estimate_lambda_ppi(obs, treated_not_k, control_not_k,
                                    unobs, m1, m0)

    # Compute estimate on fold k using lambda from opposite folds
    est_k <- compute_ppi_estimate(obs, unobs, treated_in_k, control_in_k,
                                   lambda_k)

    fold_estimates[k] <- est_k
    fold_lambdas[k] <- lambda_k
  }

  # Average across folds (equal weights)
  estimate <- mean(fold_estimates)
  lambda <- mean(fold_lambdas)

  # Compute variance using delta method
  variance <- compute_ppi_variance(msd, lambda, n_folds)

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
    method = paste0("PPI++ (cross-fit, K=", n_folds, ")"),
    lambda = lambda,
    n1 = n1,
    n0 = n0,
    m = m1 + m0,
    conf_level = conf_level,
    additional = list(
      m1 = m1,
      m0 = m0,
      fold_lambdas = fold_lambdas,
      fold_estimates = fold_estimates
    )
  )
}

#' Estimate lambda for PPI++ (single lambda across arms)
#' @keywords internal
estimate_lambda_ppi <- function(obs, treated_idx, control_idx, unobs, m1, m0) {
  # Get data for estimation
  Y1 <- obs$Y[treated_idx]
  Y0 <- obs$Y[control_idx]
  S1_obs <- obs$S1[treated_idx]
  S0_obs <- obs$S0[control_idx]

  n1 <- length(Y1)
  n0 <- length(Y0)

  # Lambda minimizes: sum_d [ lambda^2 * C_d - 2 * lambda * D_d ]
  # where C_d = var(S) * (1/m_d + 1/n_d)
  #       D_d = cov(Y, S) / n_d

  # Compute C and D for each arm
  var_S1 <- var(S1_obs)
  var_S0 <- var(S0_obs)
  cov_YS_1 <- cov(Y1, S1_obs)
  cov_YS_0 <- cov(Y0, S0_obs)

  C1 <- var_S1 * (1/m1 + 1/n1)
  C0 <- var_S0 * (1/m0 + 1/n0)
  D1 <- cov_YS_1 / n1
  D0 <- cov_YS_0 / n0

  # Combined optimal lambda: (D1 + D0) / (C1 + C0)
  C_total <- C1 + C0
  D_total <- D1 + D0

  if (C_total > 0) {
    lambda <- D_total / C_total
  } else {
    lambda <- 0
  }

  return(lambda)
}

#' Compute PPI++ estimate for a fold
#' @keywords internal
compute_ppi_estimate <- function(obs, unobs, treated_idx, control_idx, lambda) {
  # Get fold data
  Y1 <- obs$Y[treated_idx]
  Y0 <- obs$Y[control_idx]
  S1_obs <- obs$S1[treated_idx]
  S0_obs <- obs$S0[control_idx]

  # Get unlabeled means by arm
  S1_unobs_mean <- mean(unobs$S1[unobs$D == 1])
  S0_unobs_mean <- mean(unobs$S0[unobs$D == 0])

  # PPI estimate for each arm
  # mu_d = Y_bar + lambda * (S_U_bar - S_O_bar)
  mu_1 <- mean(Y1) + lambda * (S1_unobs_mean - mean(S1_obs))
  mu_0 <- mean(Y0) + lambda * (S0_unobs_mean - mean(S0_obs))

  return(mu_1 - mu_0)
}

#' Compute PPI++ variance using delta method
#' @keywords internal
compute_ppi_variance <- function(data, lambda, n_folds) {
  obs <- data$observed
  unobs <- data$unobserved

  # Sample sizes
  n1 <- sum(obs$D == 1)
  n0 <- sum(obs$D == 0)
  m1 <- sum(unobs$D == 1)
  m0 <- sum(unobs$D == 0)

  # Per-fold sizes
  n1_k <- n1 / n_folds
  n0_k <- n0 / n_folds

  # Get data
  Y1 <- obs$Y[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]
  S1_obs <- obs$S1[obs$D == 1]
  S0_obs <- obs$S0[obs$D == 0]
  S1_unobs <- unobs$S1[unobs$D == 1]
  S0_unobs <- unobs$S0[unobs$D == 0]

  # Variance components
  var_Y1 <- var(Y1)
  var_Y0 <- var(Y0)
  var_S1_obs <- var(S1_obs)
  var_S0_obs <- var(S0_obs)
  var_S1_unobs <- var(S1_unobs)
  var_S0_unobs <- var(S0_unobs)
  cov_YS_1 <- cov(Y1, S1_obs)
  cov_YS_0 <- cov(Y0, S0_obs)

  # Arm 1 variance (using per-fold sample size for labeled part)
  var_1 <- var_Y1 / n1_k +
           lambda^2 * var_S1_obs / n1_k +
           lambda^2 * var_S1_unobs / m1 -
           2 * lambda * cov_YS_1 / n1_k

  # Arm 0 variance
  var_0 <- var_Y0 / n0_k +
           lambda^2 * var_S0_obs / n0_k +
           lambda^2 * var_S0_unobs / m0 -
           2 * lambda * cov_YS_0 / n0_k

  # Total variance (accounting for cross-fit averaging)
  # For K=2 folds with equal weights, the variance is approximately halved
  variance <- (var_1 + var_0) / n_folds

  return(max(variance, 0))  # Ensure non-negative
}

#' @title D-T (Doubly-Tuned) Estimator
#' @description D-T estimator with arm-specific tuning parameters for ATE.

#' D-T (Doubly-Tuned) Estimator
#'
#' Computes the Doubly-Tuned (D-T) estimator for the average treatment
#' effect (ATE). This estimator uses arm-specific tuning parameters (lambda_1
#' and lambda_0) estimated via cross-fitting.
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
#'   \item{lambda}{Vector of arm-specific tuning parameters (lambda_1, lambda_0)}
#'
#' @details
#' The D-T estimator uses arm-specific tuning parameters:
#' \deqn{\hat{\mu}_d^{D-T}(\lambda_d) = \bar{Y}_{\mathcal{O}_d} +
#'   \lambda_d(\bar{S}^{(d)}_{\mathcal{U}_d} - \bar{S}^{(d)}_{\mathcal{O}_d})}
#'
#' Each lambda_d is chosen to minimize the variance in arm d:
#' \deqn{\lambda_d^* = \frac{Cov(Y(d), S^{(d)}) / n_d}{Var(S^{(d)})(1/m_d + 1/n_d)}}
#'
#' The tuning parameters are estimated via cross-fitting to avoid bias.
#'
#' @note
#' D-T differs from PPI++ by using separate tuning parameters for each arm,
#' which can improve efficiency when the prediction quality differs between
#' treatment and control.
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
#' result <- msd_dt(msd)
#'
#' # Using formula interface
#' result2 <- msd_dt(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
#'
#' @export
msd_dt <- function(formula_or_data, data = NULL, observed = NULL,
                   unobserved = NULL, n_folds = 2, conf_level = 0.95,
                   seed = NULL) {

  # Resolve flexible input to msd_data
  msd <- resolve_msd_data(formula_or_data, data, observed, unobserved)

  if (!msd$has_S0 && !msd$has_S1) {
    stop("D-T requires predictions. Neither S0 nor S1 found in data.")
  }

  if (is.null(msd$unobserved) || msd$m == 0) {
    stop("D-T requires unlabeled data. Use msd_dim() for labeled-only estimation.")
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
    stop("D-T requires unlabeled units in both treatment arms.")
  }

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
    lambdas_k <- estimate_lambda_dt(obs, treated_not_k, control_not_k,
                                     unobs, m1, m0)

    # Compute estimate on fold k using lambdas from opposite folds
    est_k <- compute_dt_estimate(obs, unobs, treated_in_k, control_in_k,
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
  variance <- compute_dt_variance(msd, lambda1, lambda0, n_folds)

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
    method = paste0("D-T (Doubly-Tuned, cross-fit, K=", n_folds, ")"),
    lambda = c(lambda1, lambda0),
    n1 = n1,
    n0 = n0,
    m = m1 + m0,
    conf_level = conf_level,
    additional = list(
      m1 = m1,
      m0 = m0,
      fold_lambda1 = fold_lambda1,
      fold_lambda0 = fold_lambda0,
      fold_estimates = fold_estimates
    )
  )
}

#' Estimate arm-specific lambdas for D-T
#' @noRd
estimate_lambda_dt <- function(obs, treated_idx, control_idx, unobs, m1, m0) {
  # Get data for treatment arm
  Y1 <- obs$Y[treated_idx]
  S1_obs <- obs$S1[treated_idx]
  n1 <- length(Y1)

  # Get data for control arm
  Y0 <- obs$Y[control_idx]
  S0_obs <- obs$S0[control_idx]
  n0 <- length(Y0)

  # Lambda_1: minimize variance in arm 1
  # lambda_1^* = D_1 / C_1 where
  # C_1 = var(S1) * (1/m1 + 1/n1)
  # D_1 = cov(Y1, S1) / n1
  var_S1 <- var(S1_obs)
  cov_YS_1 <- cov(Y1, S1_obs)
  C1 <- var_S1 * (1/m1 + 1/n1)
  D1 <- cov_YS_1 / n1

  if (C1 > 0) {
    lambda1 <- D1 / C1
  } else {
    lambda1 <- 0
  }

  # Lambda_0: minimize variance in arm 0
  var_S0 <- var(S0_obs)
  cov_YS_0 <- cov(Y0, S0_obs)
  C0 <- var_S0 * (1/m0 + 1/n0)
  D0 <- cov_YS_0 / n0

  if (C0 > 0) {
    lambda0 <- D0 / C0
  } else {
    lambda0 <- 0
  }

  list(lambda1 = lambda1, lambda0 = lambda0)
}

#' Compute D-T estimate for a fold
#' @noRd
compute_dt_estimate <- function(obs, unobs, treated_idx, control_idx,
                                 lambda1, lambda0) {
  # Get fold data
  Y1 <- obs$Y[treated_idx]
  Y0 <- obs$Y[control_idx]
  S1_obs <- obs$S1[treated_idx]
  S0_obs <- obs$S0[control_idx]

  # Get unlabeled means by arm
  S1_unobs_mean <- mean(unobs$S1[unobs$D == 1])
  S0_unobs_mean <- mean(unobs$S0[unobs$D == 0])

  # D-T estimate for each arm (with arm-specific lambda)
  mu_1 <- mean(Y1) + lambda1 * (S1_unobs_mean - mean(S1_obs))
  mu_0 <- mean(Y0) + lambda0 * (S0_unobs_mean - mean(S0_obs))

  return(mu_1 - mu_0)
}

#' Compute D-T variance using delta method
#'
#' See compute_ppi_variance for the derivation. The labeled terms' K factors
#' cancel when averaging across folds; the unobserved term is shared across
#' folds and is NOT divided by K.
#' @noRd
compute_dt_variance <- function(data, lambda1, lambda0, n_folds) {
  obs <- data$observed
  unobs <- data$unobserved

  # Sample sizes
  n1 <- sum(obs$D == 1)
  n0 <- sum(obs$D == 0)
  m1 <- sum(unobs$D == 1)
  m0 <- sum(unobs$D == 0)

  # Get data
  Y1 <- obs$Y[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]
  S1_obs <- obs$S1[obs$D == 1]
  S0_obs <- obs$S0[obs$D == 0]
  S1_unobs <- unobs$S1[unobs$D == 1]
  S0_unobs <- unobs$S0[unobs$D == 0]

  # Labeled component per arm: [Var(Y) + lambda^2*Var(S) - 2*lambda*Cov(Y,S)] / n_d
  var_1_labeled <- (var(Y1) + lambda1^2 * var(S1_obs) - 2 * lambda1 * cov(Y1, S1_obs)) / n1
  var_0_labeled <- (var(Y0) + lambda0^2 * var(S0_obs) - 2 * lambda0 * cov(Y0, S0_obs)) / n0

  # Unobserved component per arm: lambda_d^2 * Var(S_d) / m_d (shared, not divided by K)
  var_1_unobs <- lambda1^2 * var(S1_unobs) / m1
  var_0_unobs <- lambda0^2 * var(S0_unobs) / m0

  # Total variance
  variance <- (var_1_labeled + var_1_unobs) + (var_0_labeled + var_0_unobs)

  return(max(variance, 0))  # Ensure non-negative
}

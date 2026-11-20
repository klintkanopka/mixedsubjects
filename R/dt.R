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
#'   \item{variance}{Estimated variance of the cross-fit estimator (see the
#'     Variance section)}
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
#' With cross-fitting, the arm mean is computed on each fold using a tuning
#' estimated on the other folds' labeled data, and the fold-specific arm means
#' are combined with weights equal to each fold's share of the arm's labeled
#' units, \eqn{\pi_d^{(k)} = n_d^{(k)} / n_d}. The weights are returned in
#' \code{additional$fold_weights}. They equal \eqn{1/K} whenever an arm's
#' labeled sample size is divisible by the number of folds.
#'
#' The tuning parameters are estimated via cross-fitting to avoid bias.
#'
#' @section Variance:
#' The reported variance is an analytic estimate of the finite-sample variance
#' of the cross-fit estimator, which averages K fold estimates that share the
#' unlabeled pool. It has three parts: (i) the variance of each fold estimate,
#' evaluated at the tuning actually used in that fold with that fold's own
#' sample moments (these are independent of the tuning, so estimation of the
#' tuning is accounted for without approximation); (ii) the covariance between
#' fold estimates induced by the shared unlabeled pool; and (iii) a
#' second-order cross-fold term that arises because the tuning used in one
#' fold is estimated on the data that the other fold is evaluated on. Part
#' (iii) is the product of the covariances between each tuning estimate and the
#' labeled prediction means of the fold it was estimated on; it is estimated in
#' closed form from third central moments, is zero when outcomes and
#' predictions are jointly symmetric (for example Gaussian), and is of order
#' \eqn{n^{-2}}. The components are returned in
#' \code{additional$variance_components}. The fixed-tuning plug-in used in
#' version 1.0.0, which omits the variance of the estimated tuning and part
#' (iii), is returned in \code{additional$variance_fixed_lambda}.
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

  # Storage for fold-by-arm estimates, lambdas and fold weights
  fold_mu1 <- numeric(n_folds)
  fold_mu0 <- numeric(n_folds)
  fold_lambda1 <- numeric(n_folds)
  fold_lambda0 <- numeric(n_folds)
  w1 <- numeric(n_folds)
  w0 <- numeric(n_folds)

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

    # Arm means on fold k using lambdas from opposite folds
    mu_k <- compute_dt_arm_means(obs, unobs, treated_in_k, control_in_k,
                                 lambdas_k$lambda1, lambdas_k$lambda0)

    fold_mu1[k] <- mu_k$mu1
    fold_mu0[k] <- mu_k$mu0
    fold_lambda1[k] <- lambdas_k$lambda1
    fold_lambda0[k] <- lambdas_k$lambda0

    # Each fold is weighted by its share of the arm's labeled units
    w1[k] <- length(treated_in_k) / n1
    w0[k] <- length(control_in_k) / n0
  }

  # Combine folds within arm: tau = sum_k (w_1k mu_1k - w_0k mu_0k)
  fold_estimates <- fold_mu1 - fold_mu0
  estimate <- sum(w1 * fold_mu1) - sum(w0 * fold_mu0)
  lambda1 <- sum(w1 * fold_lambda1)
  lambda0 <- sum(w0 * fold_lambda0)

  # Cross-fit variance (SI section S2.3)
  vc <- compute_dt_crossfit_variance(msd, fold_info, fold_lambda1, fold_lambda0,
                                     w1, w0)
  variance <- vc$variance

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
      fold_estimates = fold_estimates,
      fold_weights = list(treated = w1, control = w0),
      variance_components = vc[c("fold_variances", "cross_unlabeled", "cross_fold")],
      variance_fixed_lambda = compute_dt_variance(msd, lambda1, lambda0, n_folds)
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

#' Compute the D-T arm means for a fold
#' @return list(mu1, mu0)
#' @noRd
compute_dt_arm_means <- function(obs, unobs, treated_idx, control_idx,
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
  list(
    mu1 = mean(Y1) + lambda1 * (S1_unobs_mean - mean(S1_obs)),
    mu0 = mean(Y0) + lambda0 * (S0_unobs_mean - mean(S0_obs))
  )
}

#' Cross-fit variance for D-T
#'
#' Arm-specific tunings lambda_d = (cov_d / N_d) / (var_d * a_d) depend on arm
#' d only, so Cov(lambda_d, Sbar_e) = 0 for e != d and
#' Cov(lambda_d, Sbar_d) ~ g_d / (var_d * a_d), with
#' g_d = kappa_YSS,d / N_d^2 - lambda_d * a_d * kappa_SSS,d / N_d.
#' @noRd
compute_dt_crossfit_variance <- function(data, fold_info, fold_lambda1,
                                         fold_lambda0, w1 = NULL, w0 = NULL) {
  obs <- data$observed
  unobs <- data$unobserved
  K <- fold_info$n_folds

  Y1 <- obs$Y[obs$D == 1]; S1 <- obs$S1[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]; S0 <- obs$S0[obs$D == 0]
  S1_unobs <- unobs$S1[unobs$D == 1]
  S0_unobs <- unobs$S0[unobs$D == 0]
  m1 <- length(S1_unobs)
  m0 <- length(S0_unobs)

  Sigma_U <- diag(c(safe_var(S1_unobs, safe_var(S1)) / m1,
                    safe_var(S0_unobs, safe_var(S0)) / m0))

  k1 <- arm_third_moments(Y1, S1)
  k0 <- arm_third_moments(Y0, S0)
  v <- c(safe_var(S1), safe_var(S0))
  cv <- c(stats::cov(Y1, S1), stats::cov(Y0, S0))
  kYSS <- c(k1$kappa_YSS, k0$kappa_YSS)
  kSSS <- c(k1$kappa_SSS, k0$kappa_SSS)

  G <- zero_G(K)
  for (k in seq_len(K)) {
    N <- tuning_sizes(fold_info, k)
    a <- 1 / c(m1, m0) + 1 / N
    b <- v * a
    ok <- N > 0 & is.finite(b) & b > 0
    lam <- ifelse(ok, (cv / N) / b, 0)
    g <- kYSS / N^2 - lam * a * kSSS / N
    G[[k]] <- diag(ifelse(ok, g / b, 0))
  }

  crossfit_variance(obs, fold_info, fold_lambda1, fold_lambda0, Sigma_U, G,
                    w1 = w1, w0 = w0)
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

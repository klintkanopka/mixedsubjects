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
#'   \item{variance}{Estimated variance of the cross-fit estimator (see the
#'     Variance section)}
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
#' @section Tuning:
#' The tuning parameter used in each fold is estimated from the other folds'
#' labeled data only. The unlabeled pool enters the tuning rule only through
#' its size. In version 1.0.0 the variance of \eqn{S^{(1)} - S^{(0)}} in the
#' tuning rule was estimated from the unlabeled pool; it is now estimated from
#' the labeled tuning units, so that the tuning is independent of the data it
#' is applied to and the cross-fit estimator is unbiased in finite samples.
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
    lambda_k <- estimate_lambda_dip_pp(obs, treated_not_k, control_not_k, m)

    # Compute estimate on fold k using lambda from opposite folds
    est_k <- compute_dip_pp_estimate(obs, unobs, treated_in_k, control_in_k,
                                      lambda_k)

    fold_estimates[k] <- est_k
    fold_lambdas[k] <- lambda_k
  }

  # Average across folds (equal weights)
  estimate <- mean(fold_estimates)
  lambda <- mean(fold_lambdas)

  # Cross-fit variance (SI section S2.5)
  vc <- compute_dip_pp_crossfit_variance(msd, fold_info, fold_lambdas)
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
    method = paste0("DiP++ (cross-fit, K=", n_folds, ")"),
    lambda = lambda,
    n1 = n1,
    n0 = n0,
    m = m,
    conf_level = conf_level,
    additional = list(
      fold_lambdas = fold_lambdas,
      fold_estimates = fold_estimates,
      variance_components = vc[c("fold_variances", "cross_unlabeled", "cross_fold")],
      variance_fixed_lambda = compute_dip_pp_variance(msd, lambda, n_folds)
    )
  )
}

#' Estimate lambda for DiP++ (single lambda)
#'
#' Uses labeled data only. The variance of the within-unit prediction
#' difference S1 - S0, which enters with weight 1/m, is estimated from the
#' labeled tuning units of both arms (every labeled unit carries both
#' predictions), not from the unlabeled pool. The tuning is therefore
#' independent of the unlabeled pool and of the fold it is applied to, which
#' is what finite-sample unbiasedness of the cross-fit estimator requires.
#'
#' @param obs Observed data frame
#' @param treated_idx,control_idx Row indices of the labeled tuning units
#' @param m Size of the unlabeled pool (enters only as the known factor 1/m)
#' @noRd
estimate_lambda_dip_pp <- function(obs, treated_idx, control_idx, m) {
  # Get data
  Y1 <- obs$Y[treated_idx]
  Y0 <- obs$Y[control_idx]
  S1_obs_1 <- obs$S1[treated_idx]
  S0_obs_0 <- obs$S0[control_idx]

  n1 <- length(Y1)
  n0 <- length(Y0)

  cov_S <- labeled_pred_cov(obs, c(treated_idx, control_idx))
  var_S_diff <- cov_S[1, 1] + cov_S[2, 2] - 2 * cov_S[1, 2]
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
#' @noRd
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

#' Cross-fit variance for DiP++
#'
#' The unlabeled term uses the pooled within-unit difference S1 - S0, so
#' Sigma_U is the full covariance matrix of (S1, S0) over the unlabeled pool
#' divided by m. The tuning is lambda = D / C with D = sum_d cov_d / N_d and
#' C = var(S1 - S0) / m + sum_d var_d / N_d, all moments from labeled data,
#' giving Cov(lambda, Sbar_e) ~ g_e / C with
#' g_e = (kappa_YSS,e - lambda * kappa_SSS,e) / N_e^2, up to terms of relative
#' order n / m from the 1/m-weighted second moment.
#' @noRd
compute_dip_pp_crossfit_variance <- function(data, fold_info, fold_lambdas) {
  obs <- data$observed
  unobs <- data$unobserved
  K <- fold_info$n_folds
  m <- nrow(unobs)

  Y1 <- obs$Y[obs$D == 1]; S1 <- obs$S1[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]; S0 <- obs$S0[obs$D == 0]

  Sigma_U <- if (m > 1) {
    stats::cov(cbind(unobs$S1, unobs$S0)) / m
  } else {
    matrix(0, 2, 2)
  }
  # Labeled estimate of Var(S1 - S0), matching the tuning rule
  cov_S <- labeled_pred_cov(obs, seq_len(nrow(obs)))
  var_S_diff <- cov_S[1, 1] + cov_S[2, 2] - 2 * cov_S[1, 2]

  k1 <- arm_third_moments(Y1, S1)
  k0 <- arm_third_moments(Y0, S0)
  v <- c(safe_var(S1), safe_var(S0))
  cv <- c(stats::cov(Y1, S1), stats::cov(Y0, S0))
  kYSS <- c(k1$kappa_YSS, k0$kappa_YSS)
  kSSS <- c(k1$kappa_SSS, k0$kappa_SSS)

  G <- zero_G(K)
  for (k in seq_len(K)) {
    N <- tuning_sizes(fold_info, k)
    C <- var_S_diff / m + sum(v / N)
    if (all(N > 0) && is.finite(C) && C > 0) {
      lam <- sum(cv / N) / C
      g <- (kYSS - lam * kSSS) / N^2
      G[[k]] <- rbind(g / C, g / C)
    }
  }

  crossfit_variance(obs, fold_info, fold_lambdas, fold_lambdas, Sigma_U, G)
}


#' Compute DiP++ variance using delta method
#'
#' The unobserved component var_U = lambda^2 * Var(S1-S0) / m is shared across
#' all folds and is NOT divided by K. The labeled components' K factors cancel
#' when averaging across folds.
#' @noRd
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

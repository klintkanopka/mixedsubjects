# Variance of the cross-fit (tuned) estimators
#
# Internal helpers shared by msd_ppi(), msd_dt(), msd_dip_pp() and
# msd_dt_dip(). Every tuned estimator in the package has the form
#
#   tau_hat = sum_k (w_1k * mu_1k - w_0k * mu_0k),
#   mu_dk   = Ybar_dk + lambda_dk * (Sbar_dU - Sbar_dk),
#
# where lambda_dk is the tuning for arm d used in fold k, estimated on the
# other folds (PPI++ and DiP++ set lambda_1k = lambda_0k), and w_dk are fold
# weights that sum to one within each arm. PPI++, DiP++ and D-T DiP weight
# folds equally (w_dk = 1/K); D-T weights each fold by its share of the arm's
# labeled units (w_dk = n_dk / n_d). The variance is
#
#   Var(tau_hat) = sum_k sum_d w_dk^2 * Var(Y_d - lambda_dk S_d) / n_dk
#                  + sum_{j,k} E[ l_j' Sigma_U l_k ]
#                  + sum_{j != k} R_jk,           l_k = (w_1k lambda_1k, -w_0k lambda_0k)
#
# with Sigma_U the covariance matrix of the unlabeled means (Sbar_1U, Sbar_0U).
# R_jk is the second-order term that arises because the tuning used in fold j
# is estimated on data that include fold k, and vice versa:
#
#   R_jk = sum_{d,e} s_d s_e w_dj w_ek Cov(lambda_dj, Sbar_ek) Cov(lambda_ek, Sbar_dj),
#
# s_1 = 1, s_0 = -1 (exact for K = 2, leading order for K > 2). Everything
# except R_jk is estimated without any approximation of the tuning
# distribution: each fold's variance is evaluated at the realized tuning with
# that fold's own moments, which are independent of that tuning. R_jk is
# estimated in closed form from third central moments (delta method). See
# Supplementary Information section S2.

#' Unbiased estimate of a third central cross-moment
#'
#' Returns the k-statistic estimating E[(x - Ex)(y - Ey)(z - Ez)], or 0 when
#' fewer than three observations are available.
#' @noRd
third_moment <- function(x, y, z) {
  n <- length(x)
  if (n < 3) return(0)
  k3 <- n / ((n - 1) * (n - 2)) *
    sum((x - mean(x)) * (y - mean(y)) * (z - mean(z)))
  if (is.finite(k3)) k3 else 0
}

#' Arm-level third moments entering the cross-fold term
#'
#' @return list with kappa_YSS = k3(Y, S, S) and kappa_SSS = k3(S, S, S)
#' @noRd
arm_third_moments <- function(Y, S) {
  list(kappa_YSS = third_moment(Y, S, S), kappa_SSS = third_moment(S, S, S))
}

#' Sample variance / covariance that fall back to a supplied value
#' @noRd
safe_var <- function(x, fallback = 0) {
  if (length(x) < 2) return(fallback)
  v <- stats::var(x)
  if (is.finite(v)) v else fallback
}

#' Assemble the cross-fit variance estimate
#'
#' @param obs Observed data frame (columns Y, D, S1, S0)
#' @param fold_info Fold assignments from split_by_arm_fold()
#' @param lambda1,lambda0 Length-K vectors: tuning applied in each fold to the
#'   treated and control predictions
#' @param Sigma_U 2 x 2 estimated covariance matrix of the unlabeled prediction
#'   means (Sbar_1U, Sbar_0U)
#' @param G List of K 2 x 2 matrices. `G[[k]][d, e]` estimates the covariance
#'   between the tuning for arm d used in fold k and the labeled prediction
#'   mean of arm e in any one of the other folds (rows/cols ordered treated,
#'   control).
#' @param w1,w0 Length-K fold weights for the treated and control arms, each
#'   summing to one. Default is equal weights 1/K.
#' @return list(variance, fold_variances, cross_unlabeled, cross_fold).
#'   `fold_variances` holds each fold's contribution (its labeled part plus
#'   its own unlabeled part, fold weights applied), so that
#'   `variance = sum(fold_variances) + cross_unlabeled + cross_fold`.
#' @noRd
crossfit_variance <- function(obs, fold_info, lambda1, lambda0, Sigma_U, G,
                              w1 = NULL, w0 = NULL) {
  K <- fold_info$n_folds
  sgn <- c(1, -1)
  if (is.null(w1)) w1 <- rep(1 / K, K)
  if (is.null(w0)) w0 <- rep(1 / K, K)

  Y1_all <- obs$Y[fold_info$treated_idx]
  S1_all <- obs$S1[fold_info$treated_idx]
  Y0_all <- obs$Y[fold_info$control_idx]
  S0_all <- obs$S0[fold_info$control_idx]

  V_fold <- numeric(K)
  for (k in seq_len(K)) {
    in1 <- fold_info$treated_folds == k
    in0 <- fold_info$control_folds == k
    r1 <- Y1_all[in1] - lambda1[k] * S1_all[in1]
    r0 <- Y0_all[in0] - lambda0[k] * S0_all[in0]

    # Fold-k residual variances are independent of the fold-k tuning. A fold
    # with a single unit in an arm borrows that arm's full-sample variance.
    v1 <- safe_var(r1, fallback = safe_var(Y1_all - lambda1[k] * S1_all))
    v0 <- safe_var(r0, fallback = safe_var(Y0_all - lambda0[k] * S0_all))

    l_k <- c(w1[k] * lambda1[k], -w0[k] * lambda0[k])
    V_fold[k] <- w1[k]^2 * v1 / max(sum(in1), 1) +
      w0[k]^2 * v0 / max(sum(in0), 1) +
      drop(t(l_k) %*% Sigma_U %*% l_k)
  }

  cross_U <- 0
  cross_G <- 0
  if (K > 1) {
    for (j in seq_len(K)) {
      for (k in seq_len(K)) {
        if (j == k) next
        l_j <- c(w1[j] * lambda1[j], -w0[j] * lambda0[j])
        l_k <- c(w1[k] * lambda1[k], -w0[k] * lambda0[k])
        cross_U <- cross_U + drop(t(l_j) %*% Sigma_U %*% l_k)
        # sum_{d,e} s_d s_e w_dj w_ek G_j[d, e] G_k[e, d]
        W <- outer(c(w1[j], w0[j]), c(w1[k], w0[k]))
        cross_G <- cross_G + sum(outer(sgn, sgn) * W * G[[j]] * t(G[[k]]))
      }
    }
  }

  variance <- sum(V_fold) + cross_U + cross_G

  list(
    variance = max(variance, 0),
    fold_variances = V_fold,
    cross_unlabeled = cross_U,
    cross_fold = cross_G
  )
}

#' Labeled-sample covariance matrix of the paired predictions (S1, S0)
#'
#' Every labeled unit carries both predictions, and under random assignment
#' and random labeling their joint distribution is the same in both arms and
#' in the unlabeled pool. The DiP tuning rules use this labeled estimate so
#' that the tuning does not depend on the unlabeled pool.
#'
#' @param obs Observed data frame (columns S1, S0)
#' @param idx Row indices of the labeled units to use (both arms pooled)
#' @return 2 x 2 covariance matrix ordered (S1, S0); zeros if fewer than two
#'   units or if it cannot be computed
#' @noRd
labeled_pred_cov <- function(obs, idx) {
  if (length(idx) < 2) return(matrix(0, 2, 2))
  cv <- stats::cov(cbind(obs$S1[idx], obs$S0[idx]))
  if (all(is.finite(cv))) cv else matrix(0, 2, 2)
}

#' Sizes of the tuning sample (all folds but k), by arm
#' @noRd
tuning_sizes <- function(fold_info, k) {
  c(N1 = sum(fold_info$treated_folds != k),
    N0 = sum(fold_info$control_folds != k))
}

#' Zero matrices for the cross-fold term (used when it cannot be estimated)
#' @noRd
zero_G <- function(K) {
  replicate(K, matrix(0, 2, 2), simplify = FALSE)
}

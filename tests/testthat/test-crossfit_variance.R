test_that("third_moment is the unbiased k-statistic and handles tiny samples", {
  x <- c(1, 2, 4, 8, 16)
  n <- length(x)
  expect_equal(third_moment(x, x, x),
               n / ((n - 1) * (n - 2)) * sum((x - mean(x))^3))
  expect_equal(third_moment(c(1, 2), c(1, 2), c(1, 2)), 0)
  # symmetric sample -> zero third moment
  expect_equal(third_moment(c(-2, -1, 0, 1, 2), c(-2, -1, 0, 1, 2), c(-2, -1, 0, 1, 2)), 0)
})

test_that("cross-fit variance equals the sum of its reported components", {
  td <- make_test_data()
  for (f in list(msd_ppi, msd_dt, msd_dip_pp, msd_dt_dip)) {
    for (K in c(2, 3, 5)) {
      r <- f(td$msd, n_folds = K, seed = 7)
      vc <- r$additional$variance_components
      expect_equal(r$variance,
                   sum(vc$fold_variances) + vc$cross_unlabeled + vc$cross_fold)
      expect_length(vc$fold_variances, K)
      expect_true(is.finite(r$additional$variance_fixed_lambda))
    }
  }
})

test_that("PPI++ two-fold variance matches the closed form computed by hand", {
  td <- make_test_data()
  r <- msd_ppi(td$msd, seed = 11)
  obs <- td$msd$observed
  un <- td$msd$unobserved
  fi <- split_by_arm_fold(td$msd, 2, 11)
  lam <- r$additional$fold_lambdas
  m1 <- sum(un$D == 1); m0 <- sum(un$D == 0)
  VU <- var(un$S1[un$D == 1]) / m1 + var(un$S0[un$D == 0]) / m0

  Y1 <- obs$Y[fi$treated_idx]; S1 <- obs$S1[fi$treated_idx]
  Y0 <- obs$Y[fi$control_idx]; S0 <- obs$S0[fi$control_idx]

  # (i) fold-wise fixed-lambda variances at the realized tunings
  Vk <- sapply(1:2, function(k) {
    a <- fi$treated_folds == k; b <- fi$control_folds == k
    var(Y1[a] - lam[k] * S1[a]) / sum(a) + var(Y0[b] - lam[k] * S0[b]) / sum(b) +
      lam[k]^2 * VU
  })
  # (iii) Gamma for the tuning used in fold k (estimated on the other fold)
  gam <- sapply(1:2, function(k) {
    N1 <- sum(fi$treated_folds != k); N0 <- sum(fi$control_folds != k)
    a1 <- 1 / m1 + 1 / N1; a0 <- 1 / m0 + 1 / N0
    b <- var(S1) * a1 + var(S0) * a0
    l <- (cov(Y1, S1) / N1 + cov(Y0, S0) / N0) / b
    g1 <- third_moment(Y1, S1, S1) / N1^2 - l * a1 * third_moment(S1, S1, S1) / N1
    g0 <- third_moment(Y0, S0, S0) / N0^2 - l * a0 * third_moment(S0, S0, S0) / N0
    (g1 - g0) / b
  })
  expected <- sum(Vk) / 4 + lam[1] * lam[2] * VU / 2 + gam[1] * gam[2] / 2
  expect_equal(r$variance, expected)
  expect_equal(r$additional$variance_components$cross_fold, gam[1] * gam[2] / 2)
})

test_that("cross-fold term vanishes when third moments are zero", {
  # Mirror every labeled unit around the arm means: all third central moments are 0
  set.seed(5)
  h <- 30
  base <- data.frame(S1 = rnorm(h), S0 = rnorm(h), e = rnorm(h))
  mk <- function(d) {
    x <- rbind(base, -base)
    data.frame(D = d, S1 = x$S1, S0 = x$S0,
               Y = d * 0.3 + 0.7 * ifelse(d == 1, x$S1, x$S0) + x$e)
  }
  obs <- rbind(mk(1), mk(0))
  un <- data.frame(D = rep(c(1, 0), each = 150), S1 = rnorm(300), S0 = rnorm(300))
  msd <- msd_data(observed = obs, unobserved = un)
  for (f in list(msd_ppi, msd_dt, msd_dip_pp, msd_dt_dip)) {
    r <- f(msd, seed = 2)
    expect_equal(r$additional$variance_components$cross_fold, 0, tolerance = 1e-12)
  }
})

test_that("cross-fold term is non-negative for two equal folds of the doubly tuned estimator", {
  set.seed(9)
  n <- 80; m <- 400
  obs <- data.frame(D = rep(c(1, 0), each = n / 2), S1 = rexp(n))
  obs$S0 <- 0.8 * obs$S1 + rnorm(n, 0, 0.2)
  obs$Y <- obs$S1^2 + rnorm(n)
  un <- data.frame(D = rep(c(1, 0), each = m / 2), S1 = rexp(m))
  un$S0 <- 0.8 * un$S1 + rnorm(m, 0, 0.2)
  r <- msd_dt(msd_data(observed = obs, unobserved = un), seed = 1)
  expect_gte(r$additional$variance_components$cross_fold, 0)
  # skewed data: the new variance exceeds the fixed-lambda plug-in
  expect_gt(r$variance, r$additional$variance_fixed_lambda)
})

test_that("variance is finite when a fold has a single unit in an arm", {
  set.seed(3)
  obs <- data.frame(D = rep(c(1, 0), each = 7), S1 = rnorm(14), S0 = rnorm(14))
  obs$Y <- obs$S1 + rnorm(14)
  un <- data.frame(D = rep(c(1, 0), each = 20), S1 = rnorm(40), S0 = rnorm(40))
  msd <- msd_data(observed = obs, unobserved = un)
  for (f in list(msd_ppi, msd_dt, msd_dip_pp, msd_dt_dip)) {
    r <- f(msd, n_folds = 5, seed = 4)
    expect_true(is.finite(r$variance))
    expect_gt(r$se, 0)
  }
})

test_that("D-T weights folds by their share of each arm's labeled units", {
  # 51 treated and 49 control units: two folds cannot be equal in either arm
  set.seed(21)
  obs <- data.frame(D = c(rep(1, 51), rep(0, 49)), S1 = rexp(100))
  obs$S0 <- 0.7 * obs$S1 + rnorm(100, 0, 0.3)
  obs$Y <- 0.4 * obs$D + obs$S1^2 + rnorm(100)
  un <- data.frame(D = rep(c(1, 0), each = 150), S1 = rexp(300))
  un$S0 <- 0.7 * un$S1 + rnorm(300, 0, 0.3)
  msd <- msd_data(observed = obs, unobserved = un)

  r <- msd_dt(msd, seed = 5)
  fi <- split_by_arm_fold(msd, 2, 5)
  w1 <- as.numeric(table(fi$treated_folds)) / 51
  w0 <- as.numeric(table(fi$control_folds)) / 49
  expect_equal(r$additional$fold_weights$treated, w1)
  expect_equal(r$additional$fold_weights$control, w0)
  expect_false(isTRUE(all.equal(w1, c(0.5, 0.5))))

  # Point estimate: sum_k (w_1k mu_1k - w_0k mu_0k), computed by hand
  l1 <- r$additional$fold_lambda1; l0 <- r$additional$fold_lambda0
  Y1 <- obs$Y[fi$treated_idx]; S1 <- obs$S1[fi$treated_idx]
  Y0 <- obs$Y[fi$control_idx]; S0 <- obs$S0[fi$control_idx]
  S1U <- un$S1[un$D == 1]; S0U <- un$S0[un$D == 0]
  mu1 <- sapply(1:2, function(k) {
    a <- fi$treated_folds == k
    mean(Y1[a]) + l1[k] * (mean(S1U) - mean(S1[a]))
  })
  mu0 <- sapply(1:2, function(k) {
    b <- fi$control_folds == k
    mean(Y0[b]) + l0[k] * (mean(S0U) - mean(S0[b]))
  })
  expect_equal(r$estimate, sum(w1 * mu1) - sum(w0 * mu0))

  # Variance: the plug-in estimator of the SI, arm by arm, computed by hand
  arm_var <- function(Y, S, SU, folds, lam, w) {
    n <- length(Y); mU <- length(SU)
    VU <- var(SU) / mU
    own <- sum(sapply(1:2, function(k) {
      a <- folds == k
      w[k]^2 * (var(Y[a] - lam[k] * S[a]) / sum(a) + lam[k]^2 * VU)
    }))
    gam <- sapply(1:2, function(k) {
      N <- sum(folds != k)
      a <- 1 / mU + 1 / N
      b <- var(S) * a
      l <- (cov(Y, S) / N) / b
      (third_moment(Y, S, S) / N^2 - l * a * third_moment(S, S, S) / N) / b
    })
    own + 2 * w[1] * w[2] * (lam[1] * lam[2] * VU + gam[1] * gam[2])
  }
  expected <- arm_var(Y1, S1, S1U, fi$treated_folds, l1, w1) +
    arm_var(Y0, S0, S0U, fi$control_folds, l0, w0)
  expect_equal(r$variance, expected)
})

test_that("D-T fold weights are 1/K when arm sizes divide evenly", {
  td <- make_test_data()
  r <- msd_dt(td$msd, n_folds = 2, seed = 3)
  expect_equal(r$additional$fold_weights$treated, c(0.5, 0.5))
  expect_equal(r$additional$fold_weights$control, c(0.5, 0.5))
  expect_equal(r$estimate, mean(r$additional$fold_estimates))
})

test_that("DiP tunings use labeled data only", {
  td <- make_test_data()
  # Replace the unlabeled predictions with very different values (same size):
  # the tunings must not move, because the pool enters them only through m.
  un2 <- td$unobs_df
  un2$S1 <- 5 * un2$S1 + rexp(nrow(un2))
  un2$S0 <- -3 * un2$S0 + rnorm(nrow(un2), 0, 4)
  msd2 <- msd_data(observed = td$obs_df, unobserved = un2)

  a <- msd_dip_pp(td$msd, seed = 8); b <- msd_dip_pp(msd2, seed = 8)
  expect_equal(a$additional$fold_lambdas, b$additional$fold_lambdas)

  a <- msd_dt_dip(td$msd, seed = 8); b <- msd_dt_dip(msd2, seed = 8)
  expect_equal(a$additional$fold_lambda1, b$additional$fold_lambda1)
  expect_equal(a$additional$fold_lambda0, b$additional$fold_lambda0)
})

test_that("DiP++ tuning matches the labeled-only rule computed by hand", {
  td <- make_test_data()
  obs <- td$msd$observed
  m <- nrow(td$msd$unobserved)
  r <- msd_dip_pp(td$msd, seed = 13)
  fi <- split_by_arm_fold(td$msd, 2, 13)
  lam <- sapply(1:2, function(k) {
    t_idx <- fi$treated_idx[fi$treated_folds != k]
    c_idx <- fi$control_idx[fi$control_folds != k]
    idx <- c(t_idx, c_idx)
    a <- var(obs$S1[idx] - obs$S0[idx]) / m +
      var(obs$S1[t_idx]) / length(t_idx) + var(obs$S0[c_idx]) / length(c_idx)
    cc <- cov(obs$Y[t_idx], obs$S1[t_idx]) / length(t_idx) +
      cov(obs$Y[c_idx], obs$S0[c_idx]) / length(c_idx)
    cc / a
  })
  expect_equal(r$additional$fold_lambdas, lam)
})

test_that("D-T DiP tuning matches the labeled-only rule computed by hand", {
  td <- make_test_data()
  obs <- td$msd$observed
  m <- nrow(td$msd$unobserved)
  r <- msd_dt_dip(td$msd, seed = 17)
  fi <- split_by_arm_fold(td$msd, 2, 17)
  lam <- sapply(1:2, function(k) {
    t_idx <- fi$treated_idx[fi$treated_folds != k]
    c_idx <- fi$control_idx[fi$control_folds != k]
    idx <- c(t_idx, c_idx)
    N1 <- length(t_idx); N0 <- length(c_idx)
    M <- matrix(c(var(obs$S1[t_idx]) * (1 / m + 1 / N1),
                  -cov(obs$S1[idx], obs$S0[idx]) / m,
                  -cov(obs$S1[idx], obs$S0[idx]) / m,
                  var(obs$S0[c_idx]) * (1 / m + 1 / N0)), 2, 2)
    D <- c(cov(obs$Y[t_idx], obs$S1[t_idx]) / N1,
           cov(obs$Y[c_idx], obs$S0[c_idx]) / N0)
    solve(M, D)
  })
  expect_equal(r$additional$fold_lambda1, lam[1, ])
  expect_equal(r$additional$fold_lambda0, lam[2, ])
})

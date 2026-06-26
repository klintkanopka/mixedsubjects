test_that("D-T DiP returns msd_result class", {
  td <- make_test_data()
  result <- msd_dt_dip(td$msd, seed = 123)

  expect_s3_class(result, "msd_result")
})

test_that("D-T DiP lambda is length-2 vector", {
  td <- make_test_data()
  result <- msd_dt_dip(td$msd, seed = 123)

  expect_true(is.numeric(result$lambda))
  expect_length(result$lambda, 2)
})

test_that("D-T DiP estimate is finite and SE > 0", {
  td <- make_test_data()
  result <- msd_dt_dip(td$msd, seed = 123)

  expect_true(is.finite(result$estimate))
  expect_gt(result$se, 0)
})

# --- coupled lambda solver (regression for cross-term bug) ---

test_that("solve_dt_dip_lambda satisfies the coupled 2x2 system", {
  args <- list(
    var_S1_unobs = 1.3, var_S0_unobs = 0.9, cov_S1_S0_unobs = 0.5,
    var_S1_obs = 1.1, var_S0_obs = 0.8,
    cov_YS_1 = 0.7, cov_YS_0 = 0.4,
    n1 = 40, n0 = 60, m = 200
  )
  lam <- do.call(mixedsubjects:::solve_dt_dip_lambda, args)

  A1 <- args$var_S1_unobs / args$m + args$var_S1_obs / args$n1
  A0 <- args$var_S0_unobs / args$m + args$var_S0_obs / args$n0
  B  <- args$cov_S1_S0_unobs / args$m
  D1 <- args$cov_YS_1 / args$n1
  D0 <- args$cov_YS_0 / args$n0

  expect_equal(lam$lambda1 * A1 - lam$lambda0 * B, D1)
  expect_equal(lam$lambda0 * A0 - lam$lambda1 * B, D0)
})

test_that("solve_dt_dip_lambda reduces to independent solution when uncoupled", {
  # With zero cross-covariance the arms decouple to lambda_d = D_d / A_d
  lam <- mixedsubjects:::solve_dt_dip_lambda(
    var_S1_unobs = 1.2, var_S0_unobs = 0.9, cov_S1_S0_unobs = 0,
    var_S1_obs = 1.2, var_S0_obs = 0.9,
    cov_YS_1 = 0.6, cov_YS_0 = 0.3,
    n1 = 50, n0 = 50, m = 100
  )
  expect_equal(lam$lambda1, (0.6 / 50) / (1.2 / 100 + 1.2 / 50))
  expect_equal(lam$lambda0, (0.3 / 50) / (0.9 / 100 + 0.9 / 50))
})

test_that("solve_dt_dip_lambda returns zeros for degenerate (zero-variance) input", {
  lam <- mixedsubjects:::solve_dt_dip_lambda(
    var_S1_unobs = 0, var_S0_unobs = 0, cov_S1_S0_unobs = 0,
    var_S1_obs = 0, var_S0_obs = 0,
    cov_YS_1 = 0, cov_YS_0 = 0,
    n1 = 10, n0 = 10, m = 10
  )
  expect_equal(lam$lambda1, 0)
  expect_equal(lam$lambda0, 0)
})

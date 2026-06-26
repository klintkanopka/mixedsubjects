# --- bootstrap_variance tests ---

test_that("bootstrap_variance returns list with expected fields", {
  td <- make_test_data()
  result <- bootstrap_variance(td$msd, "dim", n_bootstrap = 50, seed = 42)

  expect_true(is.list(result))
  expect_true("estimate" %in% names(result))
  expect_true("variance" %in% names(result))
  expect_true("se" %in% names(result))
  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true("bootstrap_estimates" %in% names(result))
})

test_that("bootstrap SE > 0", {
  td <- make_test_data()
  result <- bootstrap_variance(td$msd, "dim", n_bootstrap = 50, seed = 42)

  expect_gt(result$se, 0)
})

test_that("bootstrap ci_lower < ci_upper", {
  td <- make_test_data()
  result <- bootstrap_variance(td$msd, "dim", n_bootstrap = 50, seed = 42)

  expect_lt(result$ci_lower, result$ci_upper)
})

test_that("bootstrap works for each estimator string", {
  td <- make_test_data()

  for (est in c("dim", "greg", "ppi", "dt", "dip", "dip_pp", "dt_dip")) {
    result <- bootstrap_variance(td$msd, est, n_bootstrap = 20, seed = 42)
    expect_gt(result$se, 0, label = paste("SE for", est))
  }
})

test_that("bootstrap errors on invalid estimator", {
  td <- make_test_data()
  expect_error(
    bootstrap_variance(td$msd, "invalid_name", n_bootstrap = 10),
    "arg"
  )
})

test_that("bootstrap errors on non-msd_data input", {
  expect_error(
    bootstrap_variance(data.frame(x = 1), "dim"),
    "must be an msd_data object"
  )
})

test_that("stratified resampling keeps both unobserved arms populated", {
  # With a tiny treated unobserved arm, non-stratified resampling would
  # occasionally draw zero treated units and fail arm-split estimators like
  # GREG. Arm-stratified resampling always preserves m1 treated / m0 control,
  # so every bootstrap sample succeeds.
  set.seed(7)
  n_obs <- 60
  D_obs <- rep(c(1, 0), each = n_obs / 2)
  Y_obs <- rnorm(n_obs) + 0.5 * D_obs
  obs_df <- data.frame(
    Y = Y_obs, D = D_obs,
    S0 = 0.4 * Y_obs + rnorm(n_obs, 0, 0.4),
    S1 = 0.6 * Y_obs + rnorm(n_obs, 0, 0.3)
  )
  # Heavily imbalanced unobserved pool: only 3 treated units
  D_unobs <- c(rep(1, 3), rep(0, 120))
  unobs_df <- data.frame(
    D = D_unobs,
    S0 = rnorm(length(D_unobs), 0, 0.5),
    S1 = rnorm(length(D_unobs), 0, 0.4)
  )
  msd <- msd_data(observed = obs_df, unobserved = unobs_df)

  expect_no_warning(
    result <- bootstrap_variance(msd, "greg", n_bootstrap = 200, seed = 1)
  )
  expect_equal(result$n_bootstrap, 200)
})

# --- estimate_all tests ---

test_that("estimate_all returns msd_summary data.frame", {
  td <- make_test_data()
  result <- estimate_all(td$msd)

  expect_s3_class(result, "msd_summary")
  expect_s3_class(result, "data.frame")
})

test_that("estimate_all includes DiM row always", {
  td <- make_test_data()
  result <- estimate_all(td$msd)

  expect_true(any(grepl("DiM", result$Estimator)))
})

test_that("estimate_all with full data includes all 7 estimators", {
  td <- make_test_data()
  result <- estimate_all(td$msd)

  expect_equal(nrow(result), 7)
})

test_that("estimate_all with obs-only data returns only DiM", {
  td <- make_obs_only_data()
  result <- estimate_all(td$msd)

  expect_equal(nrow(result), 1)
  expect_true(grepl("DiM", result$Estimator[1]))
})

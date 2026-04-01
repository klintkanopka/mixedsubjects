test_that("GREG returns msd_result class", {
  td <- make_test_data()
  result <- msd_greg(td$msd)

  expect_s3_class(result, "msd_result")
})

test_that("GREG lambda is always 1", {
  td <- make_test_data()
  result <- msd_greg(td$msd)

  expect_equal(result$lambda, 1)
})

test_that("GREG estimate is finite numeric", {
  td <- make_test_data()
  result <- msd_greg(td$msd)

  expect_true(is.finite(result$estimate))
})

test_that("GREG SE > 0", {
  td <- make_test_data()
  result <- msd_greg(td$msd)

  expect_gt(result$se, 0)
})

test_that("GREG CI lower < estimate < CI upper", {
  td <- make_test_data()
  result <- msd_greg(td$msd)

  expect_lt(result$ci_lower, result$estimate)
  expect_gt(result$ci_upper, result$estimate)
})

test_that("GREG recovers reasonable ATE on simulated data", {
  td <- make_test_data()
  result <- msd_greg(td$msd)

  expect_equal(result$estimate, td$true_tau, tolerance = 0.1)
})

test_that("GREG errors with no predictions", {
  td <- make_no_predictions_data()
  # Need to add unobserved data since no-predictions also has no unobserved
  unobs_df <- data.frame(D = c(1, 0))
  msd <- msd_data(observed = td$obs_df, unobserved = unobs_df)
  # This will still have no S0/S1 columns
  expect_error(
    msd_greg(msd),
    "requires predictions"
  )
})

test_that("GREG errors with no unlabeled data", {
  td <- make_obs_only_data()
  expect_error(
    msd_greg(td$msd),
    "requires unlabeled data"
  )
})

test_that("GREG errors when unlabeled in only one arm", {
  obs_df <- data.frame(
    Y = c(1.2, 1.4, 0.8, 0.6),
    S0 = c(1.0, 1.2, 0.7, 0.5),
    S1 = c(1.1, 1.3, 0.9, 0.7),
    D = c(1, 1, 0, 0)
  )
  unobs_df <- data.frame(
    S0 = c(1.1, 0.9),
    S1 = c(1.2, 1.0),
    D = c(1, 1)  # Only treatment arm
  )
  msd <- msd_data(observed = obs_df, unobserved = unobs_df)

  expect_error(
    msd_greg(msd),
    "unlabeled units in both treatment arms"
  )
})

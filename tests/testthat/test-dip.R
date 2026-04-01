test_that("DiP returns msd_result class", {
  td <- make_test_data()
  result <- msd_dip(td$msd)

  expect_s3_class(result, "msd_result")
})

test_that("DiP lambda is 1", {
  td <- make_test_data()
  result <- msd_dip(td$msd)

  expect_equal(result$lambda, 1)
})

test_that("DiP estimate is finite and SE > 0", {
  td <- make_test_data()
  result <- msd_dip(td$msd)

  expect_true(is.finite(result$estimate))
  expect_gt(result$se, 0)
})

test_that("DiP errors without both predictions", {
  # Create data with only S1 (no S0)
  obs_df <- data.frame(
    Y = c(1.2, 1.4, 0.8, 0.6),
    S1 = c(1.1, 1.3, 0.9, 0.7),
    D = c(1, 1, 0, 0)
  )
  unobs_df <- data.frame(
    S1 = c(1.2, 1.0),
    D = c(1, 0)
  )
  msd <- msd_data(observed = obs_df, unobserved = unobs_df)

  expect_error(
    msd_dip(msd),
    "requires both S0 and S1"
  )
})

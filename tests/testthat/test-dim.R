test_that("DiM estimate equals hand-computed mean difference", {
  obs_df <- data.frame(
    Y = c(3.0, 5.0, 7.0, 1.0, 2.0, 3.0),
    D = c(1, 1, 1, 0, 0, 0)
  )
  msd <- msd_data(observed = obs_df)
  result <- msd_dim(msd)

  expected_estimate <- mean(c(3.0, 5.0, 7.0)) - mean(c(1.0, 2.0, 3.0))
  expect_equal(result$estimate, expected_estimate, tolerance = 1e-10)
})

test_that("DiM variance equals var(Y1)/n1 + var(Y0)/n0", {
  obs_df <- data.frame(
    Y = c(3.0, 5.0, 7.0, 1.0, 2.0, 3.0),
    D = c(1, 1, 1, 0, 0, 0)
  )
  msd <- msd_data(observed = obs_df)
  result <- msd_dim(msd)

  Y1 <- c(3.0, 5.0, 7.0)
  Y0 <- c(1.0, 2.0, 3.0)
  expected_variance <- var(Y1) / 3 + var(Y0) / 3
  expect_equal(result$variance, expected_variance, tolerance = 1e-10)
})

test_that("DiM returns msd_result with correct fields", {
  td <- make_test_data()
  result <- msd_dim(td$msd)

  expect_s3_class(result, "msd_result")
  expect_true(is.numeric(result$estimate))
  expect_true(is.numeric(result$variance))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(is.character(result$method))
  expect_equal(result$m, 0)
  expect_null(result$lambda)
})

test_that("DiM works with formula interface Y ~ D", {
  obs_df <- data.frame(
    Y = c(3.0, 5.0, 1.0, 2.0),
    D = c(1, 1, 0, 0)
  )

  result <- msd_dim(Y ~ D, observed = obs_df)

  expect_s3_class(result, "msd_result")
  expected <- mean(c(3.0, 5.0)) - mean(c(1.0, 2.0))
  expect_equal(result$estimate, expected, tolerance = 1e-10)
})

test_that("DiM works with formula + observed dataframe", {
  df <- data.frame(
    response = c(4.0, 6.0, 1.0, 3.0),
    treated = c(1, 1, 0, 0)
  )

  result <- msd_dim(response ~ treated, observed = df)

  expect_s3_class(result, "msd_result")
  expected <- mean(c(4.0, 6.0)) - mean(c(1.0, 3.0))
  expect_equal(result$estimate, expected, tolerance = 1e-10)
})

test_that("DiM errors on non-msd_data, non-formula input", {
  expect_error(
    msd_dim("not_valid"),
    "must be an msd_data object or a formula"
  )
})

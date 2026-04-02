test_that("PPI returns msd_result class", {
  td <- make_test_data()
  result <- msd_ppi(td$msd, seed = 123)

  expect_s3_class(result, "msd_result")
})

test_that("PPI lambda is a single numeric value", {
  td <- make_test_data()
  result <- msd_ppi(td$msd, seed = 123)

  expect_true(is.numeric(result$lambda))
  expect_length(result$lambda, 1)
})

test_that("PPI estimate is finite and SE > 0", {
  td <- make_test_data()
  result <- msd_ppi(td$msd, seed = 123)

  expect_true(is.finite(result$estimate))
  expect_gt(result$se, 0)
})

test_that("PPI CI contains estimate", {
  td <- make_test_data()
  result <- msd_ppi(td$msd, seed = 123)

  expect_lt(result$ci_lower, result$estimate)
  expect_gt(result$ci_upper, result$estimate)
})

test_that("PPI works with custom n_folds", {
  td <- make_test_data()
  result <- msd_ppi(td$msd, n_folds = 3, seed = 123)

  expect_s3_class(result, "msd_result")
  expect_true(is.finite(result$estimate))
})

test_that("PPI works with custom conf_level", {
  td <- make_test_data()
  result_90 <- msd_ppi(td$msd, conf_level = 0.90, seed = 123)
  result_95 <- msd_ppi(td$msd, conf_level = 0.95, seed = 123)

  expect_equal(result_90$conf_level, 0.90)
  # 90% CI should be narrower than 95% CI
  width_90 <- result_90$ci_upper - result_90$ci_lower
  width_95 <- result_95$ci_upper - result_95$ci_lower
  expect_lt(width_90, width_95)
})

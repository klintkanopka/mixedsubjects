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

test_that("DiP++ returns msd_result class", {
  td <- make_test_data()
  result <- msd_dip_pp(td$msd, seed = 123)

  expect_s3_class(result, "msd_result")
})

test_that("DiP++ lambda is a single numeric", {
  td <- make_test_data()
  result <- msd_dip_pp(td$msd, seed = 123)

  expect_true(is.numeric(result$lambda))
  expect_length(result$lambda, 1)
})

test_that("DiP++ estimate is finite and SE > 0", {
  td <- make_test_data()
  result <- msd_dip_pp(td$msd, seed = 123)

  expect_true(is.finite(result$estimate))
  expect_gt(result$se, 0)
})

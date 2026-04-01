test_that("compute_ci returns correct bounds for known values", {
  ci <- compute_ci(estimate = 1, se = 0.5, conf_level = 0.95)

  z <- qnorm(0.975)
  expected_lower <- 1 - z * 0.5
  expected_upper <- 1 + z * 0.5

  expect_equal(ci[["ci_lower"]], expected_lower, tolerance = 1e-6)
  expect_equal(ci[["ci_upper"]], expected_upper, tolerance = 1e-6)
})

test_that("var_pop computes variance with n denominator", {
  x <- c(2, 4, 6)
  # mean = 4, deviations: -2, 0, 2, sum of squares = 8, /3 = 8/3
  expected <- 8 / 3

  expect_equal(var_pop(x), expected, tolerance = 1e-10)
})

test_that("cov_pop computes covariance with n denominator", {
  x <- c(1, 2, 3)
  y <- c(2, 4, 6)
  # mean_x = 2, mean_y = 4
  # deviations: x: -1, 0, 1; y: -2, 0, 2
  # products: 2, 0, 2; sum = 4; /3 = 4/3
  expected <- 4 / 3

  expect_equal(cov_pop(x, y), expected, tolerance = 1e-10)
})

test_that("create_folds returns vector of correct length with correct values", {
  folds <- create_folds(10, n_folds = 3, seed = 42)

  expect_length(folds, 10)
  expect_true(all(folds %in% 1:3))
})

test_that("create_folds is reproducible with seed", {
  folds1 <- create_folds(20, n_folds = 2, seed = 123)
  folds2 <- create_folds(20, n_folds = 2, seed = 123)

  expect_identical(folds1, folds2)
})

test_that("parse_msd_formula parses Y ~ D | S1 + S0 correctly", {
  parsed <- parse_msd_formula(Y ~ D | S1 + S0)

  expect_equal(parsed$outcome, "Y")
  expect_equal(parsed$treatment, "D")
  expect_equal(parsed$predictions$pred_treated, "S1")
  expect_equal(parsed$predictions$pred_control, "S0")
})

test_that("parse_msd_formula errors without | separator", {
  expect_error(
    parse_msd_formula(Y ~ D),
    "must contain '\\|' separator"
  )
})

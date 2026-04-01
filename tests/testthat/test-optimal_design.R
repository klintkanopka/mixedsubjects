test_that("optimal_design returns msd_design class", {
  td <- make_test_data()
  design <- optimal_design(td$msd, budget = 10000, cost_human = 10,
                           cost_prediction = 0.01)

  expect_s3_class(design, "msd_design")
})

test_that("optimal_n_obs >= min_observed", {
  td <- make_test_data()
  design <- optimal_design(td$msd, budget = 10000, cost_human = 10,
                           cost_prediction = 0.01, min_observed = 20)

  expect_gte(design$optimal_n_obs, 20)
})

test_that("budget_used <= budget", {
  td <- make_test_data()
  design <- optimal_design(td$msd, budget = 10000, cost_human = 10,
                           cost_prediction = 0.01)

  expect_lte(design$budget_used, 10000)
})

test_that("optimal_variance > 0", {
  td <- make_test_data()
  design <- optimal_design(td$msd, budget = 10000, cost_human = 10,
                           cost_prediction = 0.01)

  expect_gt(design$optimal_variance, 0)
})

test_that("optimal_estimator is a valid name", {
  td <- make_test_data()
  design <- optimal_design(td$msd, budget = 10000, cost_human = 10,
                           cost_prediction = 0.01)

  valid_names <- c("dim", "greg", "ppi", "dt", "dip", "dip_pp", "dt_dip")
  expect_true(design$optimal_estimator %in% valid_names)
})

# --- Error cases ---

test_that("optimal_design errors on non-msd_data", {
  expect_error(
    optimal_design(data.frame(x = 1), budget = 100, cost_human = 1,
                   cost_prediction = 0.01),
    "must be an msd_data object"
  )
})

test_that("optimal_design errors on budget <= 0", {
  td <- make_test_data()
  expect_error(
    optimal_design(td$msd, budget = -100, cost_human = 10,
                   cost_prediction = 0.01),
    "must be positive"
  )
})

test_that("optimal_design errors on cost_human <= 0", {
  td <- make_test_data()
  expect_error(
    optimal_design(td$msd, budget = 10000, cost_human = -1,
                   cost_prediction = 0.01),
    "must be positive"
  )
})

test_that("optimal_design errors on treatment_prob outside (0,1)", {
  td <- make_test_data()
  expect_error(
    optimal_design(td$msd, budget = 10000, cost_human = 10,
                   cost_prediction = 0.01, treatment_prob = 0),
    "treatment_prob must be between 0 and 1"
  )
  expect_error(
    optimal_design(td$msd, budget = 10000, cost_human = 10,
                   cost_prediction = 0.01, treatment_prob = 1),
    "treatment_prob must be between 0 and 1"
  )
})

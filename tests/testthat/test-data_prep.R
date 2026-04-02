test_that("msd_data Mode 1 (combined dataframe) creates valid msd_data", {
  set.seed(42)
  n <- 50
  combined <- data.frame(
    Y = c(rnorm(30), rep(NA, 20)),
    D = c(rep(1, 15), rep(0, 15), rep(1, 10), rep(0, 10)),
    S0 = rnorm(n),
    S1 = rnorm(n)
  )

  msd <- msd_data(data = combined)

  expect_s3_class(msd, "msd_data")
  expect_equal(msd$n1, 15)
  expect_equal(msd$n0, 15)
  expect_equal(msd$m, 20)
  expect_true(msd$has_S0)
  expect_true(msd$has_S1)
  expect_true(msd$has_both_predictions)
})

test_that("msd_data Mode 2 (separate dataframes) creates valid msd_data", {
  obs_df <- data.frame(
    Y = c(1.2, 0.8, 1.5, 0.9),
    S0 = c(1.0, 0.7, 1.3, 0.8),
    S1 = c(1.1, 0.9, 1.4, 1.0),
    D = c(1, 0, 1, 0)
  )
  unobs_df <- data.frame(
    S0 = c(1.1, 0.9),
    S1 = c(1.2, 1.0),
    D = c(1, 0)
  )

  msd <- msd_data(observed = obs_df, unobserved = unobs_df)

  expect_s3_class(msd, "msd_data")
  expect_equal(msd$n1, 2)
  expect_equal(msd$n0, 2)
  expect_equal(msd$m, 2)
})

test_that("default column names (Y, D, S0, S1) work", {
  obs_df <- data.frame(
    Y = c(1.0, 2.0, 3.0, 4.0),
    D = c(1, 1, 0, 0),
    S0 = c(0.9, 1.9, 2.8, 3.8),
    S1 = c(1.1, 2.1, 3.2, 4.2)
  )

  msd <- msd_data(observed = obs_df)

  expect_equal(msd$col_mapping$obs$outcome, "Y")
  expect_equal(msd$col_mapping$obs$treatment, "D")
  expect_equal(msd$col_mapping$obs$pred_control, "S0")
  expect_equal(msd$col_mapping$obs$pred_treated, "S1")
})

test_that("custom column names work when explicitly specified", {
  obs_df <- data.frame(
    outcome = c(1.0, 2.0, 3.0, 4.0),
    treatment = c(1, 1, 0, 0),
    pred_control = c(0.9, 1.9, 2.8, 3.8),
    pred_treated = c(1.1, 2.1, 3.2, 4.2)
  )

  msd <- msd_data(
    observed = obs_df,
    outcome = "outcome",
    treatment = "treatment",
    pred_control = "pred_control",
    pred_treated = "pred_treated"
  )

  expect_equal(msd$col_mapping$obs$outcome, "outcome")
  expect_equal(msd$col_mapping$obs$treatment, "treatment")
  expect_true(msd$has_S0)
  expect_true(msd$has_S1)
})

test_that("per-dataframe column name overrides work", {
  obs_df <- data.frame(
    my_outcome = c(1.2, 0.8),
    claude_pred_0 = c(1.0, 0.7),
    claude_pred_1 = c(1.1, 0.9),
    treat = c(1, 0)
  )
  unobs_df <- data.frame(
    s0_claude = c(1.1, 0.9),
    s1_claude = c(1.2, 1.0),
    D = c(1, 0)
  )

  msd <- msd_data(
    observed = obs_df,
    unobserved = unobs_df,
    obs_outcome = "my_outcome",
    obs_treatment = "treat",
    obs_pred_control = "claude_pred_0",
    obs_pred_treated = "claude_pred_1",
    unobs_treatment = "D",
    unobs_pred_control = "s0_claude",
    unobs_pred_treated = "s1_claude"
  )

  expect_s3_class(msd, "msd_data")
  expect_equal(msd$n1, 1)
  expect_equal(msd$n0, 1)
  expect_equal(msd$m, 2)
})

test_that("observed-only data works with m=0 and unobserved=NULL", {
  obs_df <- data.frame(
    Y = c(1.0, 2.0, 3.0, 4.0),
    D = c(1, 1, 0, 0),
    S0 = c(0.9, 1.9, 2.8, 3.8),
    S1 = c(1.1, 2.1, 3.2, 4.2)
  )

  msd <- msd_data(observed = obs_df)

  expect_equal(msd$m, 0)
  expect_null(msd$unobserved)
  expect_true(msd$has_S0)
  expect_true(msd$has_S1)
})

# --- Error cases ---

test_that("providing both data and observed errors", {
  df <- data.frame(Y = 1:4, D = c(1, 1, 0, 0))
  expect_error(
    msd_data(data = df, observed = df),
    "Provide either 'data' OR"
  )
})

test_that("providing neither data nor observed errors", {
  expect_error(
    msd_data(),
    "Must provide either 'data' or 'observed'"
  )
})

test_that("missing Y column errors with informative message", {
  df <- data.frame(x = 1:4, D = c(1, 1, 0, 0))
  expect_error(
    msd_data(observed = df),
    "Column 'Y' not found in data"
  )
})

test_that("missing treatment column errors", {
  df <- data.frame(Y = 1:4, x = c(1, 1, 0, 0))
  expect_error(
    msd_data(observed = df),
    "Column 'D' not found in data"
  )
})

test_that("non-binary D errors", {
  df <- data.frame(Y = 1:4, D = c(0, 1, 2, 3))
  expect_error(
    msd_data(observed = df),
    "Treatment indicator D must be binary"
  )
})

test_that("all treated (no control) errors", {
  df <- data.frame(Y = 1:4, D = c(1, 1, 1, 1))
  expect_error(
    msd_data(observed = df),
    "No control units found"
  )
})

test_that("all control (no treated) errors", {
  df <- data.frame(Y = 1:4, D = c(0, 0, 0, 0))
  expect_error(
    msd_data(observed = df),
    "No treated units found"
  )
})

test_that("missing Y values in observed data errors", {
  df <- data.frame(Y = c(1, NA, 3, 4), D = c(1, 1, 0, 0))
  # In Mode 2 (separate dfs), NAs in Y trigger error
  # Need to bypass the combined mode split behavior
  expect_error(
    msd_data(observed = df),
    "Observed data contains missing Y values"
  )
})

test_that("NA in treatment column errors", {
  df <- data.frame(Y = 1:4, D = c(1, NA, 0, 0))
  expect_error(
    msd_data(observed = df),
    "missing treatment assignments"
  )
})

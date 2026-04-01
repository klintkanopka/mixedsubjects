# Shared test fixtures for mixedsubjects tests
# Loaded automatically by testthat before any tests run

#' Generate full test data with observed + unobserved, both S0/S1
make_test_data <- function(seed = 42) {
  set.seed(seed)

  n_obs <- 100
  n_unobs <- 200
  true_tau <- 0.5

  # Observed data
  D_obs <- rep(c(1, 0), each = n_obs / 2)
  Y_obs <- rnorm(n_obs) + true_tau * D_obs

  # Predictions correlated with outcomes + common-mode prediction error
  # Common error captures shared LLM bias for the same unit across arms
  common_error_obs <- rnorm(n_obs, 0, 0.3)
  S1_obs <- 0.6 * Y_obs + common_error_obs + rnorm(n_obs, 0, 0.3)
  S0_obs <- 0.4 * Y_obs + common_error_obs + rnorm(n_obs, 0, 0.4)

  obs_df <- data.frame(Y = Y_obs, D = D_obs, S0 = S0_obs, S1 = S1_obs)

  # Unobserved data (same prediction pipeline as observed — Assumption C)
  D_unobs <- rep(c(1, 0), each = n_unobs / 2)
  latent_Y <- rnorm(n_unobs) + true_tau * D_unobs
  common_error_unobs <- rnorm(n_unobs, 0, 0.3)
  S1_unobs <- 0.6 * latent_Y + common_error_unobs + rnorm(n_unobs, 0, 0.3)
  S0_unobs <- 0.4 * latent_Y + common_error_unobs + rnorm(n_unobs, 0, 0.4)

  unobs_df <- data.frame(D = D_unobs, S0 = S0_unobs, S1 = S1_unobs)

  msd <- msd_data(observed = obs_df, unobserved = unobs_df)

  list(obs_df = obs_df, unobs_df = unobs_df, msd = msd, true_tau = true_tau)
}

#' Generate observed-only data (no unobserved) for DiM-only tests
make_obs_only_data <- function(seed = 42) {
  set.seed(seed)

  n_obs <- 100
  D_obs <- rep(c(1, 0), each = n_obs / 2)
  Y_obs <- rnorm(n_obs) + 0.5 * D_obs
  S0_obs <- 0.4 * Y_obs + rnorm(n_obs, 0, 0.5)
  S1_obs <- 0.6 * Y_obs + rnorm(n_obs, 0, 0.4)

  obs_df <- data.frame(Y = Y_obs, D = D_obs, S0 = S0_obs, S1 = S1_obs)
  msd <- msd_data(observed = obs_df)

  list(obs_df = obs_df, msd = msd)
}

#' Generate data with no prediction columns (Y and D only)
make_no_predictions_data <- function(seed = 42) {
  set.seed(seed)

  n_obs <- 100
  D_obs <- rep(c(1, 0), each = n_obs / 2)
  Y_obs <- rnorm(n_obs) + 0.5 * D_obs

  obs_df <- data.frame(Y = Y_obs, D = D_obs)
  msd <- msd_data(observed = obs_df)

  list(obs_df = obs_df, msd = msd)
}

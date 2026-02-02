# Test script for mixedsubjects package
# Run this script to verify all functions work correctly

# Load all package functions
pkg_dir <- "C:/Users/vanloon/OneDrive - Massachusetts Institute of Technology/MSD PO/mixedsubjects"
for (f in list.files(file.path(pkg_dir, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

cat("\n")
cat("========================================\n")
cat("Testing mixedsubjects package\n
")
cat("========================================\n\n")

# Set seed for reproducibility
set.seed(42)

# -----------------------------------------------------------------------------
# Simulate test data
# -----------------------------------------------------------------------------
cat("Generating simulated data...\n")

# True ATE
true_tau <- 0.5

# Sample sizes
n_obs <- 200   # Observed units
n_unobs <- 500 # Unobserved units

# Generate potential outcomes
# Y(0) ~ N(0, 1)
# Y(1) ~ N(tau, 1)
Y0_potential <- rnorm(n_obs + n_unobs, mean = 0, sd = 1)
Y1_potential <- rnorm(n_obs + n_unobs, mean = true_tau, sd = 1)

# Treatment assignment (balanced)
D <- c(rep(1, n_obs/2), rep(0, n_obs/2),  # Observed
       rep(1, n_unobs/2), rep(0, n_unobs/2))  # Unobserved

# Observed outcome (switching equation)
Y_full <- ifelse(D == 1, Y1_potential, Y0_potential)

# Generate predictions with correlation to true outcomes
# S^(d) = beta * Y(d) + noise
beta_1 <- 0.7  # Prediction quality for treatment
beta_0 <- 0.6  # Prediction quality for control
noise_sd <- 0.5

S1 <- beta_1 * Y1_potential + rnorm(n_obs + n_unobs, 0, noise_sd)
S0 <- beta_0 * Y0_potential + rnorm(n_obs + n_unobs, 0, noise_sd)

# Add positive correlation between S1 and S0 (for DiP benefit)
common_factor <- rnorm(n_obs + n_unobs, 0, 0.3)
S1 <- S1 + common_factor
S0 <- S0 + common_factor

# Split into observed and unobserved
obs_idx <- 1:n_obs
unobs_idx <- (n_obs + 1):(n_obs + n_unobs)

obs_df <- data.frame(
  Y = Y_full[obs_idx],
  S0 = S0[obs_idx],
  S1 = S1[obs_idx],
  D = D[obs_idx]
)

unobs_df <- data.frame(
  S0 = S0[unobs_idx],
  S1 = S1[unobs_idx],
  D = D[unobs_idx]
)

cat("  Observed units:", n_obs, "(n1 =", sum(obs_df$D), ", n0 =", sum(1-obs_df$D), ")\n")
cat("  Unobserved units:", n_unobs, "\n")
cat("  True ATE:", true_tau, "\n\n")

# -----------------------------------------------------------------------------
# Test msd_data()
# -----------------------------------------------------------------------------
cat("Testing msd_data()...\n")

# Test Mode 1: Combined dataframe
combined_df <- rbind(
  obs_df,
  data.frame(Y = NA, S0 = unobs_df$S0, S1 = unobs_df$S1, D = unobs_df$D)
)
msd_combined <- msd_data(data = combined_df)
cat("  Mode 1 (combined df): OK\n")

# Test Mode 2: Separate dataframes
msd_separate <- msd_data(observed = obs_df, unobserved = unobs_df)
cat("  Mode 2 (separate dfs): OK\n")

# Use the separate version for subsequent tests
msd <- msd_separate
print(msd)

# -----------------------------------------------------------------------------
# Test all 7 estimators
# -----------------------------------------------------------------------------
cat("\n")
cat("========================================\n")
cat("Testing all 7 estimators\n")
cat("========================================\n\n")

results <- list()

# 1. DiM
cat("1. Testing msd_dim()...\n")
tryCatch({
  results$dim <- msd_dim(msd)
  print(results$dim)
  cat("   PASSED\n\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n\n")
})

# 2. GREG
cat("2. Testing msd_greg()...\n")
tryCatch({
  results$greg <- msd_greg(msd)
  print(results$greg)
  cat("   PASSED\n\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n\n")
})

# 3. PPI++
cat("3. Testing msd_ppi()...\n")
tryCatch({
  results$ppi <- msd_ppi(msd, n_folds = 2)
  print(results$ppi)
  cat("   PASSED\n\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n\n")
})

# 4. D-T
cat("4. Testing msd_dt()...\n")
tryCatch({
  results$dt <- msd_dt(msd, n_folds = 2)
  print(results$dt)
  cat("   PASSED\n\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n\n")
})

# 5. DiP
cat("5. Testing msd_dip()...\n")
tryCatch({
  results$dip <- msd_dip(msd)
  print(results$dip)
  cat("   PASSED\n\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n\n")
})

# 6. DiP++
cat("6. Testing msd_dip_pp()...\n")
tryCatch({
  results$dip_pp <- msd_dip_pp(msd, n_folds = 2)
  print(results$dip_pp)
  cat("   PASSED\n\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n\n")
})

# 7. D-T DiP
cat("7. Testing msd_dt_dip()...\n")
tryCatch({
  results$dt_dip <- msd_dt_dip(msd, n_folds = 2)
  print(results$dt_dip)
  cat("   PASSED\n\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n\n")
})

# -----------------------------------------------------------------------------
# Summary comparison
# -----------------------------------------------------------------------------
cat("\n")
cat("========================================\n")
cat("Summary: All Estimators\n")
cat("========================================\n\n")

cat("True ATE:", true_tau, "\n\n")

summary_table <- data.frame(
  Estimator = character(),
  Estimate = numeric(),
  SE = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Covers_True = logical(),
  stringsAsFactors = FALSE
)

for (name in names(results)) {
  r <- results[[name]]
  covers <- (r$ci_lower <= true_tau) && (true_tau <= r$ci_upper)
  summary_table <- rbind(summary_table, data.frame(
    Estimator = name,
    Estimate = round(r$estimate, 4),
    SE = round(r$se, 4),
    CI_Lower = round(r$ci_lower, 4),
    CI_Upper = round(r$ci_upper, 4),
    Covers_True = covers,
    stringsAsFactors = FALSE
  ))
}

print(summary_table)

# -----------------------------------------------------------------------------
# Test estimate_all()
# -----------------------------------------------------------------------------
cat("\n")
cat("========================================\n")
cat("Testing estimate_all()\n")
cat("========================================\n")

tryCatch({
  all_est <- estimate_all(msd)
  print(all_est)
  cat("   PASSED\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n")
})

# -----------------------------------------------------------------------------
# Test optimal_design()
# -----------------------------------------------------------------------------
cat("\n")
cat("========================================\n")
cat("Testing optimal_design()\n")
cat("========================================\n\n")

tryCatch({
  design <- optimal_design(
    pilot_data = msd,
    budget = 10000,
    cost_human = 10,
    cost_prediction = 0.01,
    treatment_prob = 0.5
  )
  print(design)
  cat("   PASSED\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n")
  print(e)
})

# -----------------------------------------------------------------------------
# Test bootstrap_variance() (with fewer replications for speed)
# -----------------------------------------------------------------------------
cat("\n")
cat("========================================\n")
cat("Testing bootstrap_variance()\n")
cat("========================================\n\n")

tryCatch({
  boot_result <- bootstrap_variance(msd, "dt_dip", n_bootstrap = 100)
  cat("Bootstrap variance for D-T DiP:\n")
  cat("  Estimate:", round(boot_result$estimate, 4), "\n")
  cat("  Bootstrap SE:", round(boot_result$se, 4), "\n")
  cat("  Bootstrap 95% CI: [", round(boot_result$ci_lower, 4), ", ",
      round(boot_result$ci_upper, 4), "]\n", sep = "")
  cat("   PASSED\n")
}, error = function(e) {
  cat("   FAILED:", e$message, "\n")
})

# -----------------------------------------------------------------------------
# Final summary
# -----------------------------------------------------------------------------
cat("\n")
cat("========================================\n")
cat("TEST SUMMARY\n")
cat("========================================\n\n")

n_passed <- length(results) + 3  # estimators + estimate_all + optimal_design + bootstrap
cat("All tests completed!\n")
cat("Estimators tested:", length(results), "/ 7\n")
cat("\n")

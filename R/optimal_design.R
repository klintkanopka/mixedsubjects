#' @title Optimal Design for Mixed-Subjects Experiments
#' @description Find optimal budget allocation between human observations and predictions.

#' Optimal Design Selection
#'
#' Determines the optimal allocation of budget between human observations and
#' LLM predictions, and recommends the best estimator for the given pilot data.
#'
#' @param pilot_data An msd_data object from a pilot study
#' @param budget Total budget available (in dollars)
#' @param cost_human Cost per human observation (in dollars)
#' @param cost_prediction Cost per LLM prediction (in dollars)
#' @param treatment_prob Probability of treatment assignment (default 0.5)
#' @param estimators Which estimators to consider. Either "all" or a character
#'   vector with subset of: "dim", "greg", "ppi", "dt", "dip", "dip_pp", "dt_dip"
#' @param min_observed Minimum number of observed units required (default 20)
#' @param n_grid Number of grid points for optimization (default 100)
#'
#' @return An S3 object of class "msd_design" containing:
#'   \item{optimal_n_obs}{Recommended number of observed (human) units}
#'   \item{optimal_n_unobs}{Recommended number of unobserved (prediction) units}
#'   \item{optimal_estimator}{Recommended estimator}
#'   \item{optimal_variance}{Expected variance at the optimum}
#'   \item{optimal_se}{Expected standard error at the optimum}
#'   \item{budget_used}{Total budget used}
#'   \item{all_results}{Full grid search results for all estimators}
#'
#' @details
#' The function uses grid search to find the optimal (n_O, n_U) allocation
#' that minimizes expected variance given the budget constraint:
#'
#' \deqn{n_O \times cost_{human} + n_U \times (k \times cost_{prediction}) \leq budget}
#'
#' where k is the number of predictions per unit:
#' \itemize{
#'   \item k = 0 for DiM (no predictions)
#'   \item k = 1 for GREG, PPI++, D-T (one prediction per arm)
#'   \item k = 2 for DiP, DiP++, D-T DiP (both S^(0) and S^(1))
#' }
#'
#' The expected variance is computed using population moments estimated
#' from the pilot data.
#'
#' @examples
#' \dontrun{
#' # Pilot study data
#' pilot_obs <- data.frame(
#'   Y = rnorm(50),
#'   S0 = rnorm(50),
#'   S1 = rnorm(50),
#'   D = rep(c(1, 0), each = 25)
#' )
#' pilot_unobs <- data.frame(
#'   S0 = rnorm(100),
#'   S1 = rnorm(100),
#'   D = rep(c(1, 0), each = 50)
#' )
#' pilot <- msd_data(observed = pilot_obs, unobserved = pilot_unobs)
#'
#' # Find optimal design with $10,000 budget
#' design <- optimal_design(
#'   pilot_data = pilot,
#'   budget = 10000,
#'   cost_human = 10,      # $10 per human observation
#'   cost_prediction = 0.01 # $0.01 per prediction
#' )
#' print(design)
#' }
#'
#' @export
optimal_design <- function(pilot_data,
                            budget,
                            cost_human,
                            cost_prediction,
                            treatment_prob = 0.5,
                            estimators = "all",
                            min_observed = 20,
                            n_grid = 100) {
  # Validate inputs
  if (!inherits(pilot_data, "msd_data")) {
    stop("pilot_data must be an msd_data object")
  }

  if (budget <= 0 || cost_human <= 0 || cost_prediction < 0) {
    stop("budget and costs must be positive")
  }

  if (treatment_prob <= 0 || treatment_prob >= 1) {
    stop("treatment_prob must be between 0 and 1")
  }

  # Determine which estimators to consider
  all_estimators <- c("dim", "greg", "ppi", "dt", "dip", "dip_pp", "dt_dip")

  if (identical(estimators, "all")) {
    estimators <- all_estimators
  } else {
    estimators <- match.arg(estimators, all_estimators, several.ok = TRUE)
  }

  # Filter estimators based on data availability
  if (!pilot_data$has_both_predictions) {
    estimators <- setdiff(estimators, c("dip", "dip_pp", "dt_dip"))
  }
  if (!pilot_data$has_S0 && !pilot_data$has_S1) {
    estimators <- setdiff(estimators, c("greg", "ppi", "dt", "dip", "dip_pp", "dt_dip"))
  }

  if (length(estimators) == 0) {
    stop("No valid estimators for the given data.")
  }

  # Extract moments from pilot data
  moments <- extract_pilot_moments(pilot_data)

  # Prediction cost per unit by estimator type
  pred_costs <- list(
    dim = 0,
    greg = 1,
    ppi = 1,
    dt = 1,
    dip = 2,
    dip_pp = 2,
    dt_dip = 2
  )

  # Maximum possible observations
  max_obs <- floor(budget / cost_human)

  # Grid of n_obs values
  n_obs_grid <- seq(min_observed, max_obs, length.out = min(n_grid, max_obs - min_observed + 1))
  n_obs_grid <- unique(round(n_obs_grid))
  n_obs_grid <- n_obs_grid[n_obs_grid >= min_observed]

  # Storage for results
  results <- list()

  for (est in estimators) {
    n_preds <- pred_costs[[est]]
    cost_per_unobs <- n_preds * cost_prediction

    est_results <- data.frame(
      n_obs = integer(),
      n_unobs = integer(),
      variance = numeric(),
      budget_used = numeric()
    )

    for (n_obs in n_obs_grid) {
      # Calculate affordable n_unobs
      remaining_budget <- budget - n_obs * cost_human

      if (cost_per_unobs > 0) {
        n_unobs <- floor(remaining_budget / cost_per_unobs)
      } else {
        n_unobs <- 0  # DiM doesn't use predictions
      }

      if (n_unobs < 0) n_unobs <- 0

      # Calculate expected variance
      var_est <- expected_variance(
        estimator = est,
        n_obs = n_obs,
        n_unobs = n_unobs,
        treatment_prob = treatment_prob,
        moments = moments
      )

      budget_used <- n_obs * cost_human + n_unobs * cost_per_unobs

      est_results <- rbind(est_results, data.frame(
        n_obs = n_obs,
        n_unobs = n_unobs,
        variance = var_est,
        budget_used = budget_used
      ))
    }

    # Find optimal for this estimator
    if (nrow(est_results) > 0) {
      valid_rows <- !is.na(est_results$variance) & est_results$variance > 0
      if (any(valid_rows)) {
        best_idx <- which.min(est_results$variance[valid_rows])
        best_idx <- which(valid_rows)[best_idx]

        results[[est]] <- list(
          estimator = est,
          optimal_n_obs = est_results$n_obs[best_idx],
          optimal_n_unobs = est_results$n_unobs[best_idx],
          optimal_variance = est_results$variance[best_idx],
          budget_used = est_results$budget_used[best_idx],
          grid_results = est_results
        )
      }
    }
  }

  if (length(results) == 0) {
    stop("Could not find valid design for any estimator.")
  }

  # Find overall optimal
  min_vars <- sapply(results, function(r) r$optimal_variance)
  best_est <- names(which.min(min_vars))
  best_result <- results[[best_est]]

  # Create output object
  design <- structure(
    list(
      optimal_n_obs = best_result$optimal_n_obs,
      optimal_n_unobs = best_result$optimal_n_unobs,
      optimal_n1 = round(best_result$optimal_n_obs * treatment_prob),
      optimal_n0 = round(best_result$optimal_n_obs * (1 - treatment_prob)),
      optimal_estimator = best_est,
      optimal_variance = best_result$optimal_variance,
      optimal_se = sqrt(best_result$optimal_variance),
      budget = budget,
      budget_used = best_result$budget_used,
      cost_human = cost_human,
      cost_prediction = cost_prediction,
      treatment_prob = treatment_prob,
      all_results = results,
      pilot_moments = moments
    ),
    class = "msd_design"
  )

  return(design)
}

#' Extract population moments from pilot data
#' @keywords internal
extract_pilot_moments <- function(pilot_data) {
  obs <- pilot_data$observed
  unobs <- pilot_data$unobserved

  Y1 <- obs$Y[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]

  moments <- list(
    var_Y1 = var(Y1),
    var_Y0 = var(Y0),
    mean_Y1 = mean(Y1),
    mean_Y0 = mean(Y0)
  )

  # Prediction moments if available
  if (pilot_data$has_S1) {
    S1_obs_1 <- obs$S1[obs$D == 1]
    moments$var_S1 <- var(S1_obs_1)
    moments$cov_Y1_S1 <- cov(Y1, S1_obs_1)

    if (!is.null(unobs)) {
      moments$var_S1_unobs <- var(unobs$S1[unobs$D == 1])
    }
  }

  if (pilot_data$has_S0) {
    S0_obs_0 <- obs$S0[obs$D == 0]
    moments$var_S0 <- var(S0_obs_0)
    moments$cov_Y0_S0 <- cov(Y0, S0_obs_0)

    if (!is.null(unobs)) {
      moments$var_S0_unobs <- var(unobs$S0[unobs$D == 0])
    }
  }

  if (pilot_data$has_both_predictions && !is.null(unobs)) {
    moments$var_S_diff <- var(unobs$S1 - unobs$S0)
    moments$cov_S1_S0 <- cov(unobs$S1, unobs$S0)
  }

  return(moments)
}

#' Calculate expected variance for a given design
#' @keywords internal
expected_variance <- function(estimator, n_obs, n_unobs, treatment_prob, moments) {
  n1 <- round(n_obs * treatment_prob)
  n0 <- n_obs - n1
  m1 <- round(n_unobs * treatment_prob)
  m0 <- n_unobs - m1
  m <- n_unobs

  # Minimum sample sizes
  if (n1 < 2 || n0 < 2) return(NA)

  var_Y1 <- moments$var_Y1
  var_Y0 <- moments$var_Y0

  tryCatch({
    switch(estimator,
      "dim" = {
        # DiM variance
        var_Y1 / n1 + var_Y0 / n0
      },
      "greg" = {
        if (m1 < 1 || m0 < 1) return(NA)
        var_S1 <- moments$var_S1
        var_S0 <- moments$var_S0
        cov_YS_1 <- moments$cov_Y1_S1
        cov_YS_0 <- moments$cov_Y0_S0

        # Variance of residual Y - S
        var_resid_1 <- var_Y1 + var_S1 - 2 * cov_YS_1
        var_resid_0 <- var_Y0 + var_S0 - 2 * cov_YS_0

        var_S1 / m1 + var_S0 / m0 + var_resid_1 / n1 + var_resid_0 / n0
      },
      "ppi" = {
        if (m1 < 1 || m0 < 1) return(NA)

        # Estimate optimal lambda
        var_S1 <- moments$var_S1
        var_S0 <- moments$var_S0
        cov_YS_1 <- moments$cov_Y1_S1
        cov_YS_0 <- moments$cov_Y0_S0

        C1 <- var_S1 * (1/m1 + 1/n1)
        C0 <- var_S0 * (1/m0 + 1/n0)
        D1 <- cov_YS_1 / n1
        D0 <- cov_YS_0 / n0

        lambda <- (D1 + D0) / (C1 + C0)

        # Variance at optimal lambda
        var_1 <- var_Y1/n1 + lambda^2*var_S1*(1/m1 + 1/n1) - 2*lambda*cov_YS_1/n1
        var_0 <- var_Y0/n0 + lambda^2*var_S0*(1/m0 + 1/n0) - 2*lambda*cov_YS_0/n0

        var_1 + var_0
      },
      "dt" = {
        if (m1 < 1 || m0 < 1) return(NA)

        var_S1 <- moments$var_S1
        var_S0 <- moments$var_S0
        cov_YS_1 <- moments$cov_Y1_S1
        cov_YS_0 <- moments$cov_Y0_S0

        # Arm-specific optimal lambdas
        C1 <- var_S1 * (1/m1 + 1/n1)
        D1 <- cov_YS_1 / n1
        lambda1 <- if (C1 > 0) D1 / C1 else 0

        C0 <- var_S0 * (1/m0 + 1/n0)
        D0 <- cov_YS_0 / n0
        lambda0 <- if (C0 > 0) D0 / C0 else 0

        var_1 <- var_Y1/n1 + lambda1^2*var_S1*(1/m1 + 1/n1) - 2*lambda1*cov_YS_1/n1
        var_0 <- var_Y0/n0 + lambda0^2*var_S0*(1/m0 + 1/n0) - 2*lambda0*cov_YS_0/n0

        var_1 + var_0
      },
      "dip" = {
        if (m < 1) return(NA)

        var_S_diff <- moments$var_S_diff
        var_S1 <- moments$var_S1
        var_S0 <- moments$var_S0
        cov_YS_1 <- moments$cov_Y1_S1
        cov_YS_0 <- moments$cov_Y0_S0

        var_resid_1 <- var_Y1 + var_S1 - 2 * cov_YS_1
        var_resid_0 <- var_Y0 + var_S0 - 2 * cov_YS_0

        var_S_diff / m + var_resid_1 / n1 + var_resid_0 / n0
      },
      "dip_pp" = {
        if (m < 1) return(NA)

        var_S_diff <- moments$var_S_diff
        var_S1 <- moments$var_S1
        var_S0 <- moments$var_S0
        cov_YS_1 <- moments$cov_Y1_S1
        cov_YS_0 <- moments$cov_Y0_S0

        # Optimal lambda for DiP++
        C <- var_S_diff / m + var_S1 / n1 + var_S0 / n0
        D <- cov_YS_1 / n1 + cov_YS_0 / n0
        lambda <- if (C > 0) D / C else 0

        var_U <- lambda^2 * var_S_diff / m
        var_1 <- var_Y1/n1 + lambda^2*var_S1/n1 - 2*lambda*cov_YS_1/n1
        var_0 <- var_Y0/n0 + lambda^2*var_S0/n0 - 2*lambda*cov_YS_0/n0

        var_U + var_1 + var_0
      },
      "dt_dip" = {
        if (m < 1) return(NA)

        var_S1 <- moments$var_S1
        var_S0 <- moments$var_S0
        cov_YS_1 <- moments$cov_Y1_S1
        cov_YS_0 <- moments$cov_Y0_S0
        cov_S1_S0 <- moments$cov_S1_S0

        # TODO: Lambda bug — same issue as estimate_lambda_dt_dip() in dt_dip.R (discuss with team)
        #
        # Uses independent Cov(Y,S)/Var(S) but the correct solution requires
        # solving a coupled 2x2 system that accounts for Cov(S1,S0)/m and
        # Var(S)*(1/m + 1/n). See discussion point 6 in DISCUSSION_POINTS.txt.

        # Arm-specific lambdas
        lambda1 <- if (var_S1 > 0) cov_YS_1 / var_S1 else 0
        lambda0 <- if (var_S0 > 0) cov_YS_0 / var_S0 else 0

        # Unlabeled variance
        var_weighted_diff <- lambda1^2 * var_S1 + lambda0^2 * var_S0 -
                             2 * lambda1 * lambda0 * cov_S1_S0
        var_U <- var_weighted_diff / m

        var_1 <- var_Y1/n1 + lambda1^2*var_S1/n1 - 2*lambda1*cov_YS_1/n1
        var_0 <- var_Y0/n0 + lambda0^2*var_S0/n0 - 2*lambda0*cov_YS_0/n0

        var_U + var_1 + var_0
      },
      NA
    )
  }, error = function(e) NA)
}

#' Print method for msd_design
#' @export
print.msd_design <- function(x, digits = 4, ...) {
  cat("\n")
  cat("Optimal Mixed-Subjects Design\n")
  cat("==============================\n\n")

  cat("Recommended Design:\n")
  cat("  Estimator:        ", toupper(x$optimal_estimator), "\n")
  cat("  Observed units:   ", x$optimal_n_obs, "\n")
  cat("    - Treated:      ", x$optimal_n1, "\n")
  cat("    - Control:      ", x$optimal_n0, "\n")
  cat("  Unobserved units: ", x$optimal_n_unobs, "\n\n")

  cat("Expected Performance:\n")
  cat("  Variance: ", format(round(x$optimal_variance, digits), nsmall = digits), "\n")
  cat("  SE:       ", format(round(x$optimal_se, digits), nsmall = digits), "\n\n")

  cat("Budget:\n")
  cat("  Total:    $", format(x$budget, big.mark = ","), "\n", sep = "")
  cat("  Used:     $", format(round(x$budget_used, 2), big.mark = ","), "\n", sep = "")
  cat("  Human:    $", x$cost_human, "/observation\n", sep = "")
  cat("  Predict:  $", x$cost_prediction, "/prediction\n\n", sep = "")

  # Show comparison across estimators
  if (length(x$all_results) > 1) {
    cat("Comparison Across Estimators:\n")
    comparison <- data.frame(
      Estimator = sapply(x$all_results, function(r) toupper(r$estimator)),
      n_obs = sapply(x$all_results, function(r) r$optimal_n_obs),
      n_unobs = sapply(x$all_results, function(r) r$optimal_n_unobs),
      SE = sapply(x$all_results, function(r) round(sqrt(r$optimal_variance), digits))
    )
    print(comparison, row.names = FALSE)
    cat("\n")
  }

  invisible(x)
}

#' Summary method for msd_design
#' @export
summary.msd_design <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

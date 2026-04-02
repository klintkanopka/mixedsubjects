#' @title Variance Estimation for Mixed-Subjects Design
#' @description Variance estimation methods including fold-respecting bootstrap.

#' Bootstrap Variance Estimation
#'
#' Computes bootstrap variance estimates for MSD estimators using a
#' fold-respecting resampling procedure.
#'
#' @param data An msd_data object created by \code{\link{msd_data}}
#' @param estimator Character string specifying the estimator. One of:
#'   "dim", "greg", "ppi", "dt", "dip", "dip_pp", "dt_dip"
#' @param n_bootstrap Number of bootstrap replications (default 1000)
#' @param n_folds Number of folds for cross-fitting (default 2, used for
#'   tuned estimators)
#' @param conf_level Confidence level for confidence intervals (default 0.95)
#' @param seed Random seed for reproducibility (optional)
#'
#' @return A list containing:
#'   \item{estimate}{Point estimate from the original data}
#'   \item{variance}{Bootstrap variance estimate}
#'   \item{se}{Bootstrap standard error}
#'   \item{ci_lower, ci_upper}{Bootstrap percentile confidence interval}
#'   \item{bootstrap_estimates}{Vector of bootstrap estimates}
#'
#' @details
#' The fold-respecting bootstrap resamples within each stratum:
#' \enumerate{
#'   \item Resample observed treatment units with replacement
#'   \item Resample observed control units with replacement
#'   \item Resample unlabeled units with replacement
#'   \item Recompute the estimator on the bootstrap sample
#' }
#'
#' For cross-fit estimators, the fold assignments are regenerated for each
#' bootstrap replicate to properly account for the cross-fitting variance.
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' obs_df <- data.frame(
#'   Y = rnorm(100),
#'   S0 = rnorm(100),
#'   S1 = rnorm(100),
#'   D = rep(c(1, 0), each = 50)
#' )
#' unobs_df <- data.frame(
#'   S0 = rnorm(200),
#'   S1 = rnorm(200),
#'   D = rep(c(1, 0), each = 100)
#' )
#' msd <- msd_data(observed = obs_df, unobserved = unobs_df)
#'
#' # Bootstrap variance for D-T DiP
#' boot_result <- bootstrap_variance(msd, "dt_dip", n_bootstrap = 500)
#' print(boot_result)
#' }
#'
#' @export
bootstrap_variance <- function(data,
                                estimator = c("dim", "greg", "ppi", "dt",
                                              "dip", "dip_pp", "dt_dip"),
                                n_bootstrap = 1000,
                                n_folds = 2,
                                conf_level = 0.95,
                                seed = NULL) {
  # Validate inputs
  if (!inherits(data, "msd_data")) {
    stop("data must be an msd_data object created by msd_data()")
  }

  estimator <- match.arg(estimator)

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Get original estimate
  original_result <- run_estimator(data, estimator, n_folds)
  original_estimate <- original_result$estimate

  # Bootstrap storage
  boot_estimates <- numeric(n_bootstrap)

  # Extract data for resampling
  obs <- data$observed
  unobs <- data$unobserved

  n1 <- sum(obs$D == 1)
  n0 <- sum(obs$D == 0)
  m <- if (!is.null(unobs)) nrow(unobs) else 0

  treated_idx <- which(obs$D == 1)
  control_idx <- which(obs$D == 0)

  # Bootstrap loop
  for (b in 1:n_bootstrap) {
    # Resample within strata
    boot_treated <- sample(treated_idx, n1, replace = TRUE)
    boot_control <- sample(control_idx, n0, replace = TRUE)

    # Create bootstrap observed data
    boot_obs <- rbind(obs[boot_treated, ], obs[boot_control, ])

    # TODO: Bug — unobserved resampling not stratified by arm (discuss with team)
    #
    # For PPI++ and D-T estimators, unobserved data is split by arm (D==1 vs
    # D==0). Resampling the whole pool together can shift the m1/m0 ratio in
    # each bootstrap sample. Should stratify like the observed data:
    #   boot_unobs <- rbind(
    #     unobs[sample(which(unobs$D==1), m1, replace=TRUE), ],
    #     unobs[sample(which(unobs$D==0), m0, replace=TRUE), ]
    #   )
    # Not an issue for DiP-type estimators which pool all unobserved units.

    # Resample unlabeled
    if (m > 0) {
      boot_unobs_idx <- sample(1:m, m, replace = TRUE)
      boot_unobs <- unobs[boot_unobs_idx, ]
    } else {
      boot_unobs <- NULL
    }

    # Create bootstrap msd_data
    boot_data <- structure(
      list(
        observed = boot_obs,
        unobserved = boot_unobs,
        has_S0 = data$has_S0,
        has_S1 = data$has_S1,
        has_both_predictions = data$has_both_predictions,
        n1 = n1,
        n0 = n0,
        m = m
      ),
      class = "msd_data"
    )

    # Run estimator on bootstrap sample
    tryCatch({
      boot_result <- run_estimator(boot_data, estimator, n_folds)
      boot_estimates[b] <- boot_result$estimate
    }, error = function(e) {
      boot_estimates[b] <<- NA
    })
  }

  # Remove failed bootstrap samples
  valid_boots <- !is.na(boot_estimates)
  boot_estimates <- boot_estimates[valid_boots]
  n_valid <- length(boot_estimates)

  if (n_valid < n_bootstrap * 0.9) {
    warning(sprintf("Only %d of %d bootstrap samples succeeded.",
                    n_valid, n_bootstrap))
  }

  # Compute bootstrap statistics
  boot_variance <- var(boot_estimates)
  boot_se <- sd(boot_estimates)

  # Percentile confidence interval
  alpha <- 1 - conf_level
  ci_lower <- quantile(boot_estimates, alpha / 2, na.rm = TRUE)
  ci_upper <- quantile(boot_estimates, 1 - alpha / 2, na.rm = TRUE)

  list(
    estimate = original_estimate,
    variance = boot_variance,
    se = boot_se,
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    conf_level = conf_level,
    n_bootstrap = n_valid,
    bootstrap_estimates = boot_estimates,
    method = paste0("Bootstrap (", estimator, ", B=", n_valid, ")")
  )
}

#' Run a specific estimator
#' @keywords internal
run_estimator <- function(data, estimator, n_folds = 2) {
  switch(estimator,
    "dim" = msd_dim(data),
    "greg" = msd_greg(data),
    "ppi" = msd_ppi(data, n_folds = n_folds),
    "dt" = msd_dt(data, n_folds = n_folds),
    "dip" = msd_dip(data),
    "dip_pp" = msd_dip_pp(data, n_folds = n_folds),
    "dt_dip" = msd_dt_dip(data, n_folds = n_folds),
    stop("Unknown estimator: ", estimator)
  )
}

#' Compare variance estimates across methods
#'
#' Computes and compares variance estimates using both delta-method and
#' bootstrap for a given estimator.
#'
#' @param data An msd_data object
#' @param estimator Character string specifying the estimator
#' @param n_bootstrap Number of bootstrap replications
#' @param n_folds Number of folds for cross-fitting
#' @param seed Random seed
#'
#' @return A data frame comparing variance estimates
#'
#' @examples
#' \dontrun{
#' msd <- msd_data(observed = obs_df, unobserved = unobs_df)
#' comparison <- compare_variance_methods(msd, "dt_dip", n_bootstrap = 500)
#' print(comparison)
#' }
#'
#' @export
compare_variance_methods <- function(data,
                                      estimator,
                                      n_bootstrap = 1000,
                                      n_folds = 2,
                                      seed = NULL) {
  # Delta-method estimate
  delta_result <- run_estimator(data, estimator, n_folds)

  # Bootstrap estimate
  boot_result <- bootstrap_variance(data, estimator, n_bootstrap, n_folds,
                                     seed = seed)

  # Create comparison table
  comparison <- data.frame(
    Method = c("Delta-method", "Bootstrap"),
    Estimate = c(delta_result$estimate, boot_result$estimate),
    Variance = c(delta_result$variance, boot_result$variance),
    SE = c(delta_result$se, boot_result$se),
    CI_Lower = c(delta_result$ci_lower, boot_result$ci_lower),
    CI_Upper = c(delta_result$ci_upper, boot_result$ci_upper)
  )

  attr(comparison, "estimator") <- estimator
  attr(comparison, "n_bootstrap") <- n_bootstrap

  comparison
}

#' Estimate all available estimators
#'
#' Runs all applicable estimators on the data and returns a summary table.
#'
#' @param data An msd_data object
#' @param n_folds Number of folds for cross-fitting (default 2)
#' @param conf_level Confidence level (default 0.95)
#'
#' @return A data frame with estimates from all applicable estimators
#'
#' @examples
#' \dontrun{
#' msd <- msd_data(observed = obs_df, unobserved = unobs_df)
#' all_estimates <- estimate_all(msd)
#' print(all_estimates)
#' }
#'
#' @export
estimate_all <- function(data, n_folds = 2, conf_level = 0.95) {
  if (!inherits(data, "msd_data")) {
    stop("data must be an msd_data object")
  }

  results <- list()

  # DiM always available
  results$dim <- msd_dim(data, conf_level)

  # GREG-type if predictions available
  if ((data$has_S0 || data$has_S1) && data$m > 0) {
    tryCatch({
      results$greg <- msd_greg(data, conf_level)
    }, error = function(e) NULL)

    tryCatch({
      results$ppi <- msd_ppi(data, n_folds, conf_level)
    }, error = function(e) NULL)

    tryCatch({
      results$dt <- msd_dt(data, n_folds, conf_level)
    }, error = function(e) NULL)
  }

  # DiP-type if both predictions available
  if (data$has_both_predictions && data$m > 0) {
    tryCatch({
      results$dip <- msd_dip(data, conf_level)
    }, error = function(e) NULL)

    tryCatch({
      results$dip_pp <- msd_dip_pp(data, n_folds, conf_level)
    }, error = function(e) NULL)

    tryCatch({
      results$dt_dip <- msd_dt_dip(data, n_folds, conf_level)
    }, error = function(e) NULL)
  }

  # Create summary table
  summary_df <- data.frame(
    Estimator = character(),
    Estimate = numeric(),
    SE = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    stringsAsFactors = FALSE
  )

  for (name in names(results)) {
    r <- results[[name]]
    summary_df <- rbind(summary_df, data.frame(
      Estimator = r$method,
      Estimate = r$estimate,
      SE = r$se,
      CI_Lower = r$ci_lower,
      CI_Upper = r$ci_upper,
      stringsAsFactors = FALSE
    ))
  }

  attr(summary_df, "results") <- results
  attr(summary_df, "conf_level") <- conf_level

  class(summary_df) <- c("msd_summary", "data.frame")
  summary_df
}

#' Print method for msd_summary
#' @export
print.msd_summary <- function(x, digits = 4, ...) {
  cat("\n")
  cat("Mixed-Subjects Design: All Estimators\n")
  cat("======================================\n\n")

  conf_level <- attr(x, "conf_level")
  ci_pct <- paste0(round(conf_level * 100), "%")

  # Format the table
  x_print <- x
  x_print$Estimate <- round(x_print$Estimate, digits)
  x_print$SE <- round(x_print$SE, digits)
  x_print$CI_Lower <- round(x_print$CI_Lower, digits)
  x_print$CI_Upper <- round(x_print$CI_Upper, digits)

  names(x_print)[4:5] <- c(paste0(ci_pct, " CI Lower"),
                           paste0(ci_pct, " CI Upper"))

  print.data.frame(x_print, row.names = FALSE)
  cat("\n")

  invisible(x)
}

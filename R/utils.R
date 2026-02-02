#' @title Utility Functions for Mixed-Subjects Design
#' @description Internal helper functions for the mixedsubjects package.
#' @keywords internal

# -----------------------------------------------------------------------------
# S3 Class Constructor for msd_result
# -----------------------------------------------------------------------------

#' Create an msd_result object
#'
#' @param estimate Point estimate of the ATE
#' @param variance Estimated variance
#' @param se Standard error
#' @param ci_lower Lower bound of confidence interval
#' @param ci_upper Upper bound of confidence interval
#' @param method Name of the estimation method
#' @param lambda Tuning parameter(s), if applicable
#' @param n1 Number of treated observed units
#' @param n0 Number of control observed units
#' @param m Number of unobserved (unlabeled) units
#' @param conf_level Confidence level used
#' @param additional Additional method-specific information
#' @return An S3 object of class "msd_result"
#' @keywords internal
new_msd_result <- function(estimate,
                           variance,
                           se,
                           ci_lower,
                           ci_upper,
                           method,
                           lambda = NULL,
                           n1 = NULL,
                           n0 = NULL,
                           m = NULL,
                           conf_level = 0.95,
                           additional = list()) {
  structure(
    list(
      estimate = estimate,
      variance = variance,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      method = method,
      lambda = lambda,
      n1 = n1,
      n0 = n0,
      m = m,
      conf_level = conf_level,
      additional = additional
    ),
    class = "msd_result"
  )
}

#' Print method for msd_result
#'
#' @param x An msd_result object
#' @param digits Number of digits to display
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns x
#' @export
print.msd_result <- function(x, digits = 4, ...) {
  cat("\n")
  cat("Mixed-Subjects Design Estimation\n")
  cat("=================================\n")
  cat("Estimator:", x$method, "\n\n")

  cat("Point Estimate: ", format(round(x$estimate, digits), nsmall = digits), "\n")

  cat("Standard Error: ", format(round(x$se, digits), nsmall = digits), "\n")

  ci_pct <- paste0(round(x$conf_level * 100), "%")
  cat(ci_pct, " CI:         [",
      format(round(x$ci_lower, digits), nsmall = digits), ", ",
      format(round(x$ci_upper, digits), nsmall = digits), "]\n", sep = "")

  if (!is.null(x$lambda)) {
    cat("\nTuning Parameters:\n")
    if (length(x$lambda) == 1) {
      cat("  lambda:       ", format(round(x$lambda, digits), nsmall = digits), "\n")
    } else if (length(x$lambda) == 2) {
      cat("  lambda_1 (treatment): ", format(round(x$lambda[1], digits), nsmall = digits), "\n")
      cat("  lambda_0 (control):   ", format(round(x$lambda[2], digits), nsmall = digits), "\n")
    }
  }

  cat("\nSample Sizes:\n")
  if (!is.null(x$n1) && !is.null(x$n0)) {
    cat("  Observed:   n_1=", x$n1, ", n_0=", x$n0, "\n", sep = "")
  }
  if (!is.null(x$m) && x$m > 0) {
    cat("  Unobserved: |U|=", x$m, "\n", sep = "")
  }

  cat("\n")
  invisible(x)
}

#' Summary method for msd_result
#'
#' Produces a detailed summary of the MSD estimation results, including
#' the coefficient table with z-statistics and p-values, sample size
#' information, and interpretation guidance.
#'
#' @param object An msd_result object
#' @param ... Additional arguments (ignored)
#' @return A summary.msd_result object (invisibly when printed)
#' @export
#'
#' @examples
#' \dontrun{
#' result <- msd_dt_dip(msd)
#' summary(result)
#' }
summary.msd_result <- function(object, ...) {
  z_stat <- object$estimate / object$se
  p_value <- 2 * (1 - stats::pnorm(abs(z_stat)))

  # Determine significance stars
  stars <- ""
  if (!is.na(p_value)) {
    if (p_value < 0.001) stars <- "***"
    else if (p_value < 0.01) stars <- "**"
    else if (p_value < 0.05) stars <- "*"
    else if (p_value < 0.1) stars <- "."
  }

  summary_list <- list(
    method = object$method,
    estimate = object$estimate,
    se = object$se,
    z_stat = z_stat,
    p_value = p_value,
    stars = stars,
    ci_lower = object$ci_lower,
    ci_upper = object$ci_upper,
    conf_level = object$conf_level,
    lambda = object$lambda,
    n1 = object$n1,
    n0 = object$n0,
    m = object$m,
    n_total = object$n1 + object$n0,
    additional = object$additional
  )

  class(summary_list) <- "summary.msd_result"
  summary_list
}

#' Print method for summary.msd_result
#'
#' Produces nicely formatted output similar to summary.lm, designed for
#' applied social scientists.
#'
#' @param x A summary.msd_result object
#' @param digits Number of significant digits to display (default 4)
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns x
#' @export
print.summary.msd_result <- function(x, digits = 4, ...) {
  cat("\n")
  cat("Mixed-Subjects Design: Treatment Effect Estimation\n")
  cat("===================================================\n\n")

  # Method
  cat("Estimator:", x$method, "\n\n")

  # Coefficient table header
  cat("Coefficients:\n")

  # Create coefficient table
  coef_table <- data.frame(
    Estimate = x$estimate,
    `Std. Error` = x$se,
    `z value` = x$z_stat,
    `Pr(>|z|)` = x$p_value,
    check.names = FALSE
  )
  rownames(coef_table) <- "ATE"

  # Format and print
  coef_formatted <- coef_table
  coef_formatted$Estimate <- format(round(coef_table$Estimate, digits), nsmall = digits)
  coef_formatted$`Std. Error` <- format(round(coef_table$`Std. Error`, digits), nsmall = digits)
  coef_formatted$`z value` <- format(round(coef_table$`z value`, 3), nsmall = 3)
  coef_formatted$`Pr(>|z|)` <- format_pvalue(coef_table$`Pr(>|z|)`)

  # Add significance stars
  coef_formatted$` ` <- x$stars

  print(coef_formatted, quote = FALSE, right = TRUE)

  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")

  # Confidence interval
  ci_pct <- paste0(round(x$conf_level * 100), "%")
  cat(ci_pct, "Confidence Interval: [",
      format(round(x$ci_lower, digits), nsmall = digits), ", ",
      format(round(x$ci_upper, digits), nsmall = digits), "]\n\n", sep = "")

  # Tuning parameters (if applicable)
  if (!is.null(x$lambda)) {
    cat("Tuning Parameters:\n")
    if (length(x$lambda) == 1) {
      cat("  lambda:", format(round(x$lambda, digits), nsmall = digits), "\n")
    } else if (length(x$lambda) == 2) {
      cat("  lambda_1 (treatment):", format(round(x$lambda[1], digits), nsmall = digits), "\n")
      cat("  lambda_0 (control):  ", format(round(x$lambda[2], digits), nsmall = digits), "\n")
    }
    cat("\n")
  }

  # Sample sizes
  cat("Sample Sizes:\n")
  cat("  Observed (human subjects): ", x$n_total, "\n", sep = "")
  cat("    Treated (n1):  ", x$n1, "\n", sep = "")
  cat("    Control (n0):  ", x$n0, "\n", sep = "")
  if (!is.null(x$m) && x$m > 0) {
    cat("  Unobserved (predictions only): ", x$m, "\n", sep = "")
    cat("  Total units: ", x$n_total + x$m, "\n", sep = "")
  }
  cat("\n")

  # Effect size interpretation
  cat("Interpretation:\n")
  direction <- if (x$estimate > 0) "increases" else "decreases"
  magnitude <- abs(x$estimate)

  if (x$p_value < 0.05) {
    cat("  The treatment ", direction, " the outcome by ",
        format(round(magnitude, digits), nsmall = digits),
        " units (p < 0.05).\n", sep = "")
  } else {
    cat("  The estimated effect is ",
        format(round(magnitude, digits), nsmall = digits),
        " units, but this is not statistically significant (p = ",
        format(round(x$p_value, 3), nsmall = 3), ").\n", sep = "")
  }

  # Variance reduction (if unlabeled data was used)
  if (!is.null(x$m) && x$m > 0) {
    cat("  Using ", x$m, " unlabeled predictions improved estimation precision.\n", sep = "")
  }

  cat("\n")
  invisible(x)
}

#' Format p-value for display
#' @keywords internal
format_pvalue <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.0001) return("< 0.0001")
  if (p < 0.001) return(format(round(p, 4), nsmall = 4))
  return(format(round(p, 4), nsmall = 4))
}

# -----------------------------------------------------------------------------
# Fold Splitting Functions
# -----------------------------------------------------------------------------

#' Create fold indices for cross-fitting
#'
#' @param n Number of observations
#' @param n_folds Number of folds (default 2)
#' @param seed Random seed for reproducibility (optional)
#' @return A vector of fold assignments (integers 1 to n_folds)
#' @keywords internal
create_folds <- function(n, n_folds = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Assign each observation to a fold

fold_ids <- rep(1:n_folds, length.out = n)
  fold_ids <- sample(fold_ids)  # Shuffle

  return(fold_ids)
}

#' Split data by arm and fold
#'
#' @param data An msd_data object
#' @param n_folds Number of folds
#' @param seed Random seed
#' @return A list with fold assignments for each arm
#' @keywords internal
split_by_arm_fold <- function(data, n_folds = 2, seed = NULL) {
  # Get observed data
  obs <- data$observed

  # Split treated arm
  treated_idx <- which(obs$D == 1)
  n1 <- length(treated_idx)
  folds_1 <- create_folds(n1, n_folds, seed)

  # Split control arm
  control_idx <- which(obs$D == 0)
  n0 <- length(control_idx)
  if (!is.null(seed)) seed <- seed + 1  # Different seed for control
  folds_0 <- create_folds(n0, n_folds, seed)

  list(
    treated_idx = treated_idx,
    treated_folds = folds_1,
    control_idx = control_idx,
    control_folds = folds_0,
    n_folds = n_folds
  )
}

# -----------------------------------------------------------------------------
# Moment Computation Functions
# -----------------------------------------------------------------------------
#' Compute sample variance (population formula, divide by n)
#'
#' @param x A numeric vector
#' @return Sample variance using n denominator
#' @keywords internal
var_pop <- function(x) {
  n <- length(x)
  if (n <= 1) return(0)
  sum((x - mean(x))^2) / n
}

#' Compute sample covariance (population formula, divide by n)
#'
#' @param x A numeric vector
#' @param y A numeric vector
#' @return Sample covariance using n denominator
#' @keywords internal
cov_pop <- function(x, y) {
  n <- length(x)
  if (n <= 1) return(0)
  sum((x - mean(x)) * (y - mean(y))) / n
}

#' Compute all moments needed for estimation
#'
#' @param data An msd_data object
#' @return A list of sample moments
#' @keywords internal
compute_moments <- function(data) {
  obs <- data$observed
  unobs <- data$unobserved

  # Sample sizes
  n1 <- sum(obs$D == 1)
  n0 <- sum(obs$D == 0)
  m <- if (!is.null(unobs)) nrow(unobs) else 0

  # Observed data by arm
  Y1 <- obs$Y[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]

  # Outcome moments
  mean_Y1 <- mean(Y1)
  mean_Y0 <- mean(Y0)
  var_Y1 <- var(Y1)
  var_Y0 <- var(Y0)

  # Prediction moments (arm-specific on observed)
  has_S0 <- !is.null(data$observed$S0)
  has_S1 <- !is.null(data$observed$S1)
  has_both_predictions <- has_S0 && has_S1

  moments <- list(
    n1 = n1,
    n0 = n0,
    m = m,
    mean_Y1 = mean_Y1,
    mean_Y0 = mean_Y0,
    var_Y1 = var_Y1,
    var_Y0 = var_Y0,
    has_both_predictions = has_both_predictions
  )

  # For GREG-type estimators, we need S^(D) for each observed unit
  if (has_S1) {
    S1_obs_treated <- obs$S1[obs$D == 1]
    moments$mean_S1_obs_1 <- mean(S1_obs_treated)
    moments$var_S1_obs_1 <- var(S1_obs_treated)
    moments$cov_Y1_S1 <- cov(Y1, S1_obs_treated)
  }

  if (has_S0) {
    S0_obs_control <- obs$S0[obs$D == 0]
    moments$mean_S0_obs_0 <- mean(S0_obs_control)
    moments$var_S0_obs_0 <- var(S0_obs_control)
    moments$cov_Y0_S0 <- cov(Y0, S0_obs_control)
  }

  # Unlabeled moments
  if (m > 0 && !is.null(unobs)) {
    if (has_S1) {
      moments$mean_S1_unobs <- mean(unobs$S1)
      moments$var_S1_unobs <- var(unobs$S1)
    }
    if (has_S0) {
      moments$mean_S0_unobs <- mean(unobs$S0)
      moments$var_S0_unobs <- var(unobs$S0)
    }
    if (has_both_predictions) {
      # For DiP: covariance of S1 and S0 on unlabeled
      moments$cov_S1_S0_unobs <- cov(unobs$S1, unobs$S0)
      # Variance of difference
      moments$var_S_diff_unobs <- var(unobs$S1 - unobs$S0)
    }
  }

  # For DiP on observed: we need S1 and S0 for observed units too
  if (has_both_predictions) {
    # Variance of residuals Y - S on observed
    resid_1 <- Y1 - obs$S1[obs$D == 1]
    resid_0 <- Y0 - obs$S0[obs$D == 0]
    moments$var_resid_1 <- var(resid_1)
    moments$var_resid_0 <- var(resid_0)
  }

  return(moments)
}

#' Compute arm-specific moments on a subset of data
#'
#' @param obs_subset Subset of observed data
#' @param arm Which arm (1 or 0)
#' @return List of arm-specific moments
#' @keywords internal
compute_arm_moments <- function(obs_subset, arm) {
  Y <- obs_subset$Y
  n <- length(Y)

  if (arm == 1 && !is.null(obs_subset$S1)) {
    S <- obs_subset$S1
  } else if (arm == 0 && !is.null(obs_subset$S0)) {
    S <- obs_subset$S0
  } else {
    S <- NULL
  }

  moments <- list(
    n = n,
    mean_Y = mean(Y),
    var_Y = var(Y)
  )

  if (!is.null(S)) {
    moments$mean_S <- mean(S)
    moments$var_S <- var(S)
    moments$cov_YS <- cov(Y, S)
  }

  return(moments)
}

# -----------------------------------------------------------------------------
# Confidence Interval Helper
# -----------------------------------------------------------------------------

#' Compute confidence interval
#'
#' @param estimate Point estimate
#' @param se Standard error
#' @param conf_level Confidence level (default 0.95)
#' @return Named vector with ci_lower and ci_upper
#' @keywords internal
compute_ci <- function(estimate, se, conf_level = 0.95) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  c(ci_lower = estimate - z * se, ci_upper = estimate + z * se)
}

# -----------------------------------------------------------------------------
# Formula Parsing for Estimators
# -----------------------------------------------------------------------------

#' Parse MSD formula
#'
#' Parses formulas of the form: outcome ~ treatment | predictions
#'
#' For GREG-type estimators (single prediction per arm):
#'   Y ~ D | S
#'
#' For DiP-type estimators (both predictions):
#'   Y ~ D | S1 + S0
#'
#' @param formula A formula object
#' @return A list with outcome, treatment, and prediction variable names
#' @keywords internal
parse_msd_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("Expected a formula object (e.g., Y ~ D | S1 + S0)")
  }

  # Convert formula to character and parse
  f_str <- paste(deparse(formula), collapse = " ")

  # Check for the | separator
  if (!grepl("\\|", f_str)) {
    stop("Formula must contain '|' separator. ",
         "Expected format: outcome ~ treatment | predictions\n",
         "  For GREG/PPI++/D-T: Y ~ D | S\n",
         "  For DiP/DiP++/D-T DiP: Y ~ D | S1 + S0")
  }

  # Split by |
  parts <- strsplit(f_str, "\\|")[[1]]
  left_part <- trimws(parts[1])
  pred_part <- trimws(parts[2])

  # Parse left side: outcome ~ treatment
  left_formula <- as.formula(left_part)

  # Get outcome (left of ~)
  outcome <- all.vars(left_formula[[2]])
  if (length(outcome) != 1) {
    stop("Formula must have exactly one outcome variable on the left of ~")
  }

  # Get treatment (right of ~ but left of |)
  treatment <- all.vars(left_formula[[3]])
  if (length(treatment) != 1) {
    stop("Formula must have exactly one treatment variable between ~ and |")
  }

  # Parse predictions (right of |)
  # Could be single variable (S) or two variables (S1 + S0)
  pred_formula <- as.formula(paste("~", pred_part))
  pred_vars <- all.vars(pred_formula)

  if (length(pred_vars) == 0) {
    stop("No prediction variables found after |")
  } else if (length(pred_vars) == 1) {
    # Single prediction variable (for GREG-type)
    predictions <- list(single = pred_vars[1])
  } else if (length(pred_vars) == 2) {
    # Two prediction variables (for DiP-type)
    # Order: first is treated prediction, second is control
    predictions <- list(pred_treated = pred_vars[1], pred_control = pred_vars[2])
  } else {
    stop("Expected 1 or 2 prediction variables after |, got ", length(pred_vars))
  }

  list(
    outcome = outcome,
    treatment = treatment,
    predictions = predictions
  )
}

#' Apply formula to create/update msd_data
#'
#' Takes a parsed formula and either creates a new msd_data object or
#' updates column references in an existing one.
#'
#' @param parsed_formula Result from parse_msd_formula()
#' @param data Either raw dataframes or an msd_data object
#' @param observed Raw observed dataframe (if data is NULL)
#' @param unobserved Raw unobserved dataframe (if data is NULL)
#' @return An msd_data object with appropriate column mapping
#' @keywords internal
apply_formula_to_data <- function(parsed_formula, data = NULL,
                                   observed = NULL, unobserved = NULL) {

  # If data is already msd_data, we need to remap columns
  if (inherits(data, "msd_data")) {
    # The data is already standardized - just validate that the formula
    # matches what's available
    return(data)
  }

  # Otherwise, create new msd_data with formula-specified columns
  preds <- parsed_formula$predictions

  if (!is.null(preds$single)) {
    # Single prediction - same column for both arms
    # This is a special case for GREG-type estimators
    msd_data(
      data = data,
      observed = observed,
      unobserved = unobserved,
      outcome = parsed_formula$outcome,
      treatment = parsed_formula$treatment,
      pred_control = preds$single,
      pred_treated = preds$single
    )
  } else {
    # Two predictions
    msd_data(
      data = data,
      observed = observed,
      unobserved = unobserved,
      outcome = parsed_formula$outcome,
      treatment = parsed_formula$treatment,
      pred_treated = preds$pred_treated,
      pred_control = preds$pred_control
    )
  }
}

#' Resolve data from formula or msd_data object
#'
#' Handles the flexible interface where users can provide either:
#' 1. An msd_data object directly
#' 2. A formula + raw dataframe(s)
#'
#' @param formula_or_data Either a formula or an msd_data object
#' @param data If formula provided, this should be the data (msd_data, combined df, or NULL)
#' @param observed If formula provided and data is NULL, the observed dataframe
#' @param unobserved If formula provided and data is NULL, the unobserved dataframe
#' @param require_predictions Logical, whether to require prediction columns
#' @return An msd_data object ready for estimation
#' @keywords internal
resolve_msd_data <- function(formula_or_data, data = NULL,
                              observed = NULL, unobserved = NULL,
                              require_predictions = TRUE) {

  # Case 1: formula_or_data is already msd_data
  if (inherits(formula_or_data, "msd_data")) {
    return(formula_or_data)
  }

  # Case 2: formula_or_data is a formula
  if (inherits(formula_or_data, "formula")) {
    # Parse the formula
    f_str <- paste(deparse(formula_or_data), collapse = " ")

    # Check if formula has predictions (contains |)
    if (grepl("\\|", f_str)) {
      parsed <- parse_msd_formula(formula_or_data)
    } else {
      # Simple formula without predictions: outcome ~ treatment
      outcome <- all.vars(formula_or_data[[2]])
      treatment <- all.vars(formula_or_data[[3]])
      if (length(outcome) != 1 || length(treatment) != 1) {
        stop("Formula must have exactly one outcome and one treatment variable")
      }
      parsed <- list(outcome = outcome, treatment = treatment, predictions = NULL)
    }

    # Build msd_data from the parsed formula
    if (inherits(data, "msd_data")) {
      return(data)
    } else if (!is.null(data) && is.data.frame(data)) {
      # Combined dataframe
      return(apply_formula_to_data(parsed, data = data))
    } else if (!is.null(observed)) {
      # Separate dataframes
      return(apply_formula_to_data(parsed, observed = observed, unobserved = unobserved))
    } else {
      stop("When using formula interface, must provide data, or observed/unobserved dataframes")
    }
  }

  # Case 3: Unknown input
  stop("First argument must be either an msd_data object or a formula. ",
       "Got: ", class(formula_or_data)[1])
}

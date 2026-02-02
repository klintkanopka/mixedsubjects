#' @title Data Preparation for Mixed-Subjects Design
#' @description Create and validate data objects for MSD estimation.

#' Create an MSD data object
#'
#' Creates a validated data object for mixed-subjects design estimation.
#' Accepts either a single combined dataframe or two separate dataframes
#' for observed and unobserved units. Column names can be auto-detected,
#' explicitly specified, or provided separately for each dataframe.
#'
#' @param data A combined dataframe containing both observed and unobserved units.
#'   Observed units have non-missing Y values; unobserved units have Y = NA.
#'   If provided, \code{observed} and \code{unobserved} should be NULL.
#' @param observed A dataframe of observed (labeled) units with columns for
#'   outcome Y, predictions S0/S1, and treatment D. If provided, \code{data}
#'   should be NULL.
#' @param unobserved A dataframe of unobserved (unlabeled) units with columns
#'   for predictions S0/S1 and treatment D (no Y column needed).
#' @param outcome Name of the outcome column, or NULL for auto-detection.
#'   Common patterns detected: "Y", "y", "outcome", "response", "dependent".
#' @param treatment Name of the treatment column, or NULL for auto-detection.
#'   Common patterns detected: "D", "d", "treatment", "treat", "treated", "T", "W", "w", "Z", "z".
#' @param pred_control Name of the control prediction column, or NULL for auto-detection.
#'   Common patterns detected: "S0", "s0", "S_0", "S^0", "pred_control", "pred_0",
#'   "prediction_control", "prediction_0", "yhat_control", "yhat_0".
#' @param pred_treated Name of the treatment prediction column, or NULL for auto-detection.
#'   Common patterns detected: "S1", "s1", "S_1", "S^1", "pred_treated", "pred_1",
#'   "prediction_treated", "prediction_1", "yhat_treated", "yhat_1".
#' @param obs_outcome Outcome column name for observed data only (overrides \code{outcome}).
#' @param obs_treatment Treatment column name for observed data only (overrides \code{treatment}).
#' @param obs_pred_control Control prediction column for observed data only.
#' @param obs_pred_treated Treatment prediction column for observed data only.
#' @param unobs_treatment Treatment column name for unobserved data only.
#' @param unobs_pred_control Control prediction column for unobserved data only.
#' @param unobs_pred_treated Treatment prediction column for unobserved data only.
#'
#' @return An S3 object of class "msd_data" containing:
#'   \item{observed}{Dataframe of observed units with standardized columns Y, S0, S1, D}
#'   \item{unobserved}{Dataframe of unobserved units with standardized columns S0, S1, D (or NULL)}
#'   \item{has_S0}{Logical indicating if S0 predictions are available}
#'   \item{has_S1}{Logical indicating if S1 predictions are available}
#'   \item{has_both_predictions}{Logical indicating if both S0 and S1 are available}
#'   \item{n1}{Number of treated observed units}
#'   \item{n0}{Number of control observed units}
#'   \item{m}{Number of unobserved units}
#'   \item{col_mapping}{List of original column names used}
#'
#' @details
#' The function supports flexible column name specification:
#'
#' \strong{Auto-detection:}
#' If column names are not specified (NULL), the function will search for
#' common patterns. For example, columns named "outcome", "Y", "response",
#' or "dependent" will be recognized as the outcome variable.
#'
#' \strong{Global specification:}
#' Use \code{outcome}, \code{treatment}, \code{pred_control}, \code{pred_treated}
#' to set column names that apply to both observed and unobserved dataframes.
#'
#' \strong{Per-dataframe specification:}
#' Use \code{obs_*} and \code{unobs_*} arguments when column names differ
#' between observed and unobserved dataframes. These override global settings.
#'
#' \strong{Two input modes:}
#'
#' \emph{Mode 1: Single combined dataframe}
#' Provide a single dataframe via the \code{data} argument. The function will
#' automatically split it into observed and unobserved based on whether Y is NA.
#'
#' \emph{Mode 2: Separate dataframes}
#' Provide two separate dataframes via \code{observed} and \code{unobserved}.
#'
#' @examples
#' # Auto-detection with standard column names
#' obs_df <- data.frame(
#'   Y = c(1.2, 0.8, 1.5, 0.9),
#'   S0 = c(1.0, 0.7, 1.3, 0.8),
#'   S1 = c(1.1, 0.9, 1.4, 1.0),
#'   D = c(1, 0, 1, 0)
#' )
#' msd <- msd_data(observed = obs_df)
#'
#' # Custom column names (same in both dataframes)
#' obs_df2 <- data.frame(
#'   response = c(1.2, 0.8, 1.5, 0.9),
#'   pred_ctrl = c(1.0, 0.7, 1.3, 0.8),
#'   pred_trt = c(1.1, 0.9, 1.4, 1.0),
#'   treated = c(1, 0, 1, 0)
#' )
#' unobs_df2 <- data.frame(
#'   pred_ctrl = c(1.1, 0.9),
#'   pred_trt = c(1.2, 1.0),
#'   treated = c(1, 0)
#' )
#' msd2 <- msd_data(
#'   observed = obs_df2,
#'   unobserved = unobs_df2,
#'   outcome = "response",
#'   treatment = "treated",
#'   pred_control = "pred_ctrl",
#'   pred_treated = "pred_trt"
#' )
#'
#' # Different column names in observed vs unobserved
#' obs_df3 <- data.frame(
#'   outcome = c(1.2, 0.8),
#'   claude_pred_0 = c(1.0, 0.7),
#'   claude_pred_1 = c(1.1, 0.9),
#'   treatment = c(1, 0)
#' )
#' unobs_df3 <- data.frame(
#'   s0_claude = c(1.1, 0.9),
#'   s1_claude = c(1.2, 1.0),
#'   D = c(1, 0)
#' )
#' msd3 <- msd_data(
#'   observed = obs_df3,
#'   unobserved = unobs_df3,
#'   obs_outcome = "outcome",
#'   obs_treatment = "treatment",
#'   obs_pred_control = "claude_pred_0",
#'   obs_pred_treated = "claude_pred_1",
#'   unobs_treatment = "D",
#'   unobs_pred_control = "s0_claude",
#'   unobs_pred_treated = "s1_claude"
#' )
#'
#' @export
msd_data <- function(data = NULL,
                     observed = NULL,
                     unobserved = NULL,
                     # Global column names (apply to both)
                     outcome = NULL,
                     treatment = NULL,
                     pred_control = NULL,
                     pred_treated = NULL,
                     # Per-dataframe overrides for observed
                     obs_outcome = NULL,
                     obs_treatment = NULL,
                     obs_pred_control = NULL,
                     obs_pred_treated = NULL,
                     # Per-dataframe overrides for unobserved
                     unobs_treatment = NULL,
                     unobs_pred_control = NULL,
                     unobs_pred_treated = NULL) {

  # Validate input mode

has_combined <- !is.null(data)
  has_separate <- !is.null(observed)

  if (has_combined && has_separate) {
    stop("Provide either 'data' OR ('observed' and optionally 'unobserved'), not both.")
  }

  if (!has_combined && !has_separate) {
    stop("Must provide either 'data' or 'observed' argument.")
  }

  # Resolve column names for observed data
  if (has_combined) {
    obs_source <- data
  } else {
    obs_source <- observed
  }

  # Get effective column names for observed (per-df overrides global)
  y_col <- obs_outcome %||% outcome
  d_col_obs <- obs_treatment %||% treatment
  s0_col_obs <- obs_pred_control %||% pred_control
  s1_col_obs <- obs_pred_treated %||% pred_treated

  # Auto-detect if not specified
  y_col <- detect_column(obs_source, y_col, "outcome")
  d_col_obs <- detect_column(obs_source, d_col_obs, "treatment")
  s0_col_obs <- detect_column(obs_source, s0_col_obs, "pred_control", required = FALSE)
  s1_col_obs <- detect_column(obs_source, s1_col_obs, "pred_treated", required = FALSE)

  # Get effective column names for unobserved
  unobs_source <- if (has_combined) data else unobserved

  if (!is.null(unobs_source)) {
    d_col_unobs <- unobs_treatment %||% treatment
    s0_col_unobs <- unobs_pred_control %||% pred_control
    s1_col_unobs <- unobs_pred_treated %||% pred_treated

    # For unobserved, try to match observed column names first, then auto-detect
    d_col_unobs <- detect_column(unobs_source, d_col_unobs, "treatment",
                                  fallback = d_col_obs)
    s0_col_unobs <- detect_column(unobs_source, s0_col_unobs, "pred_control",
                                   required = FALSE, fallback = s0_col_obs)
    s1_col_unobs <- detect_column(unobs_source, s1_col_unobs, "pred_treated",
                                   required = FALSE, fallback = s1_col_obs)
  } else {
    d_col_unobs <- s0_col_unobs <- s1_col_unobs <- NULL
  }

  # Store column mapping for reference
  col_mapping <- list(
    obs = list(
      outcome = y_col,
      treatment = d_col_obs,
      pred_control = s0_col_obs,
      pred_treated = s1_col_obs
    ),
    unobs = list(
      treatment = d_col_unobs,
      pred_control = s0_col_unobs,
      pred_treated = s1_col_unobs
    )
  )

  # Process based on input mode
  if (has_combined) {
    result <- process_combined_data(data, y_col, s0_col_obs, s1_col_obs, d_col_obs)
  } else {
    result <- process_separate_data(observed, unobserved,
                                     y_col, d_col_obs, s0_col_obs, s1_col_obs,
                                     d_col_unobs, s0_col_unobs, s1_col_unobs)
  }

  # Add column mapping to result
  result$col_mapping <- col_mapping

  return(result)
}

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Detect column name with auto-detection
#'
#' @param df Dataframe to search
#' @param specified User-specified column name (or NULL)
#' @param type Type of column: "outcome", "treatment", "pred_control", "pred_treated"
#' @param required If TRUE, error when not found; if FALSE, return NULL
#' @param fallback Column name to try if auto-detection fails
#' @return Column name or NULL
#' @keywords internal
detect_column <- function(df, specified, type, required = TRUE, fallback = NULL) {
  # If user specified a name, validate and return it

  if (!is.null(specified)) {
    if (specified %in% names(df)) {
      return(specified)
    } else if (required) {
      stop("Column '", specified, "' not found in data. ",
           "Available columns: ", paste(names(df), collapse = ", "))
    } else {
      return(NULL)
    }
  }

  # Try fallback first (e.g., same name as observed)
  if (!is.null(fallback) && fallback %in% names(df)) {
    return(fallback)
  }

  # Auto-detection patterns
  patterns <- switch(type,
    outcome = c("^Y$", "^y$", "^outcome$", "^response$", "^dependent$",
                "^Outcome$", "^Response$", "^DV$", "^dv$"),
    treatment = c("^D$", "^d$", "^treatment$", "^treat$", "^treated$",
                  "^Treatment$", "^Treat$", "^Treated$",
                  "^T$", "^W$", "^w$", "^Z$", "^z$", "^assignment$"),
    pred_control = c("^S0$", "^s0$", "^S_0$", "^S\\^0$", "^S0_", "^s0_",
                     "^pred_control$", "^pred_0$", "^pred0$",
                     "^prediction_control$", "^prediction_0$",
                     "^yhat_control$", "^yhat_0$", "^yhat0$",
                     "^S_control$", "^s_control$", "^Scontrol$",
                     "_S0$", "_s0$", "_control$", "control_pred"),
    pred_treated = c("^S1$", "^s1$", "^S_1$", "^S\\^1$", "^S1_", "^s1_",
                     "^pred_treated$", "^pred_1$", "^pred1$",
                     "^prediction_treated$", "^prediction_1$",
                     "^yhat_treated$", "^yhat_1$", "^yhat1$",
                     "^S_treated$", "^s_treated$", "^Streated$",
                     "_S1$", "_s1$", "_treated$", "treated_pred"),
    character(0)
  )

  # Search for matching column
  col_names <- names(df)
  for (pattern in patterns) {
    matches <- grep(pattern, col_names, value = TRUE, ignore.case = FALSE)
    if (length(matches) == 1) {
      return(matches[1])
    } else if (length(matches) > 1) {
      # Multiple matches - use first but warn
      message("Multiple columns match pattern for ", type, ": ",
              paste(matches, collapse = ", "), ". Using '", matches[1], "'.")
      return(matches[1])
    }
  }

  # Not found
  if (required) {
    stop("Could not auto-detect ", type, " column. ",
         "Please specify explicitly. Available columns: ",
         paste(col_names, collapse = ", "))
  }

  return(NULL)
}

#' Process combined dataframe
#' @keywords internal
process_combined_data <- function(data, y_col, s0_col, s1_col, d_col) {
  # Check prediction columns
  has_S0 <- !is.null(s0_col) && s0_col %in% names(data)
  has_S1 <- !is.null(s1_col) && s1_col %in% names(data)

  # Split into observed and unobserved based on Y being NA
  is_observed <- !is.na(data[[y_col]])

  obs_data <- data[is_observed, , drop = FALSE]
  unobs_data <- data[!is_observed, , drop = FALSE]

  # Create standardized observed dataframe
  observed <- data.frame(
    Y = obs_data[[y_col]],
    D = obs_data[[d_col]]
  )

  if (has_S0) observed$S0 <- obs_data[[s0_col]]
  if (has_S1) observed$S1 <- obs_data[[s1_col]]

  # Create standardized unobserved dataframe (if any)
  if (nrow(unobs_data) > 0) {
    unobserved <- data.frame(D = unobs_data[[d_col]])
    if (has_S0) unobserved$S0 <- unobs_data[[s0_col]]
    if (has_S1) unobserved$S1 <- unobs_data[[s1_col]]
  } else {
    unobserved <- NULL
  }

  # Validate and create msd_data object
  validate_and_create(observed, unobserved, has_S0, has_S1)
}

#' Process separate dataframes
#' @keywords internal
process_separate_data <- function(observed_df, unobserved_df,
                                   y_col, d_col_obs, s0_col_obs, s1_col_obs,
                                   d_col_unobs, s0_col_unobs, s1_col_unobs) {
  # Check prediction columns in observed
  has_S0 <- !is.null(s0_col_obs) && s0_col_obs %in% names(observed_df)
  has_S1 <- !is.null(s1_col_obs) && s1_col_obs %in% names(observed_df)

  # Create standardized observed dataframe
  observed <- data.frame(
    Y = observed_df[[y_col]],
    D = observed_df[[d_col_obs]]
  )

  if (has_S0) observed$S0 <- observed_df[[s0_col_obs]]
  if (has_S1) observed$S1 <- observed_df[[s1_col_obs]]

  # Process unobserved if provided
  if (!is.null(unobserved_df) && nrow(unobserved_df) > 0) {
    if (is.null(d_col_unobs) || !d_col_unobs %in% names(unobserved_df)) {
      stop("Treatment column not found in unobserved data. ",
           "Available columns: ", paste(names(unobserved_df), collapse = ", "))
    }

    unobserved <- data.frame(D = unobserved_df[[d_col_unobs]])

    # Check for predictions in unobserved
    if (has_S0) {
      if (!is.null(s0_col_unobs) && s0_col_unobs %in% names(unobserved_df)) {
        unobserved$S0 <- unobserved_df[[s0_col_unobs]]
      } else {
        warning("S0 predictions not found in unobserved data. ",
                "Looked for column: ", s0_col_unobs %||% "(auto-detect failed)")
        has_S0 <- FALSE
      }
    }

    if (has_S1) {
      if (!is.null(s1_col_unobs) && s1_col_unobs %in% names(unobserved_df)) {
        unobserved$S1 <- unobserved_df[[s1_col_unobs]]
      } else {
        warning("S1 predictions not found in unobserved data. ",
                "Looked for column: ", s1_col_unobs %||% "(auto-detect failed)")
        has_S1 <- FALSE
      }
    }
  } else {
    unobserved <- NULL
  }

  # Validate and create msd_data object
  validate_and_create(observed, unobserved, has_S0, has_S1)
}

#' Validate data and create msd_data object
#' @keywords internal
validate_and_create <- function(observed, unobserved, has_S0, has_S1) {
  # Validate observed data
  if (nrow(observed) == 0) {
    stop("No observed units found in data.")
  }

  if (any(is.na(observed$Y))) {
    stop("Observed data contains missing Y values.")
  }

  if (any(is.na(observed$D))) {
    stop("Observed data contains missing treatment assignments.")
  }

  # Check treatment indicator is binary
  if (!all(observed$D %in% c(0, 1))) {
    stop("Treatment indicator D must be binary (0 or 1).")
  }

  # Count by arm
  n1 <- sum(observed$D == 1)
  n0 <- sum(observed$D == 0)

  if (n1 == 0) {
    stop("No treated units found in observed data.")
  }

  if (n0 == 0) {
    stop("No control units found in observed data.")
  }

  # Count unobserved
  m <- if (!is.null(unobserved)) nrow(unobserved) else 0

  # Create S3 object
  result <- structure(
    list(
      observed = observed,
      unobserved = unobserved,
      has_S0 = has_S0,
      has_S1 = has_S1,
      has_both_predictions = has_S0 && has_S1,
      n1 = n1,
      n0 = n0,
      m = m
    ),
    class = "msd_data"
  )

  return(result)
}

#' Print method for msd_data
#'
#' @param x An msd_data object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns x
#' @export
print.msd_data <- function(x, ...) {
  cat("\n")
  cat("Mixed-Subjects Design Data\n")
  cat("==========================\n\n")

  cat("Sample Sizes:\n")
  cat("  Observed (labeled):    ", x$n1 + x$n0, "\n")
  cat("    - Treated (D=1):     ", x$n1, "\n")
  cat("    - Control (D=0):     ", x$n0, "\n")
  cat("  Unobserved (unlabeled):", x$m, "\n\n")

  cat("Predictions Available:\n")
  cat("  S0 (control arm):  ", ifelse(x$has_S0, "Yes", "No"), "\n")
  cat("  S1 (treatment arm):", ifelse(x$has_S1, "Yes", "No"), "\n")

  # Show column mapping if available
  if (!is.null(x$col_mapping)) {
    cat("\nColumn Mapping (original names):\n")
    cat("  Observed:  Y=", x$col_mapping$obs$outcome %||% "NA",
        ", D=", x$col_mapping$obs$treatment %||% "NA",
        ", S0=", x$col_mapping$obs$pred_control %||% "NA",
        ", S1=", x$col_mapping$obs$pred_treated %||% "NA", "\n", sep = "")
    if (x$m > 0) {
      cat("  Unobserved: D=", x$col_mapping$unobs$treatment %||% "NA",
          ", S0=", x$col_mapping$unobs$pred_control %||% "NA",
          ", S1=", x$col_mapping$unobs$pred_treated %||% "NA", "\n", sep = "")
    }
  }

  cat("\nAvailable Estimators:\n")
  cat("  - DiM:     Yes (no predictions needed)\n")

  if (x$has_S0 && x$has_S1) {
    cat("  - GREG:    Yes\n")
    cat("  - PPI++:   Yes\n")
    cat("  - D-T:     Yes\n")
    cat("  - DiP:     Yes (both S0 and S1 available)\n")
    cat("  - DiP++:   Yes\n")
    cat("  - D-T DiP: Yes\n")
  } else if (x$has_S0 || x$has_S1) {
    cat("  - GREG:    Yes (arm-specific predictions)\n")
    cat("  - PPI++:   Yes (arm-specific predictions)\n")
    cat("  - D-T:     Yes (arm-specific predictions)\n")
    cat("  - DiP:     No (requires both S0 and S1)\n")
    cat("  - DiP++:   No (requires both S0 and S1)\n")
    cat("  - D-T DiP: No (requires both S0 and S1)\n")
  } else {
    cat("  - GREG:    No (no predictions)\n")
    cat("  - PPI++:   No (no predictions)\n")
    cat("  - D-T:     No (no predictions)\n")
    cat("  - DiP:     No (no predictions)\n")
    cat("  - DiP++:   No (no predictions)\n")
    cat("  - D-T DiP: No (no predictions)\n")
  }

  cat("\n")
  invisible(x)
}

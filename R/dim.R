#' @title Difference-in-Means Estimator
#' @description Classical difference-in-means estimator for ATE.

#' Difference-in-Means Estimator
#'
#' Computes the classical difference-in-means estimator for the average
#' treatment effect (ATE). This estimator uses only observed (labeled) data
#' and does not incorporate any predictions.
#'
#' @param formula_or_data Either an msd_data object created by \code{\link{msd_data}},
#'   or a formula of the form \code{outcome ~ treatment} (predictions not needed for DiM).
#' @param data If \code{formula_or_data} is a formula, this should be either:
#'   an msd_data object, a combined dataframe, or NULL (if using observed/unobserved).
#' @param observed If using formula with separate dataframes, the observed data.
#' @param unobserved If using formula with separate dataframes, the unobserved data.
#' @param conf_level Confidence level for the confidence interval (default 0.95)
#'
#' @return An msd_result object containing:
#'   \item{estimate}{Point estimate of the ATE: mean(Y|D=1) - mean(Y|D=0)}
#'   \item{variance}{Estimated variance: var(Y|D=1)/n1 + var(Y|D=0)/n0}
#'   \item{se}{Standard error}
#'   \item{ci_lower, ci_upper}{Confidence interval bounds}
#'   \item{method}{Name of the estimation method}
#'
#' @details
#' The difference-in-means estimator is:
#' \deqn{\hat{\tau}^{DiM} = \bar{Y}_1 - \bar{Y}_0}
#'
#' where \eqn{\bar{Y}_d} is the sample mean of outcomes in arm \eqn{d}.
#'
#' The variance is estimated as:
#' \deqn{\widehat{Var}(\hat{\tau}^{DiM}) = \frac{s^2_{Y(1)}}{n_1} + \frac{s^2_{Y(0)}}{n_0}}
#'
#' where \eqn{s^2_{Y(d)}} is the sample variance of outcomes in arm \eqn{d}.
#'
#' @examples
#' # Using msd_data object (standard interface)
#' obs_df <- data.frame(
#'   Y = c(1.2, 1.4, 0.8, 0.6, 1.1, 0.9, 1.3, 0.7),
#'   S0 = c(1.0, 1.2, 0.7, 0.5, 1.0, 0.8, 1.1, 0.6),
#'   S1 = c(1.1, 1.3, 0.9, 0.7, 1.1, 0.9, 1.2, 0.8),
#'   D = c(1, 1, 0, 0, 1, 0, 1, 0)
#' )
#' msd <- msd_data(observed = obs_df)
#' 
#' result <- msd_dim(msd)
#' print(result)
#'
#' # Using formula interface with custom column names
#' df <- data.frame(
#'   response = c(1.2, 1.4, 0.8, 0.6),
#'   treated = c(1, 1, 0, 0)
#' )
#' result2 <- msd_dim(response ~ treated, observed = df)
#'
#' @export
msd_dim <- function(formula_or_data, data = NULL, observed = NULL,
                    unobserved = NULL, conf_level = 0.95) {

  # Handle flexible interface
  if (inherits(formula_or_data, "formula")) {
    # Parse formula - for DiM we only need outcome and treatment
    f_str <- paste(deparse(formula_or_data), collapse = " ")

    # DiM doesn't require |, so handle both cases
    if (grepl("\\|", f_str)) {
      parsed <- parse_msd_formula(formula_or_data)
    } else {
      # Simple formula: outcome ~ treatment
      outcome <- all.vars(formula_or_data[[2]])
      treatment <- all.vars(formula_or_data[[3]])
      if (length(outcome) != 1 || length(treatment) != 1) {
        stop("Formula must have exactly one outcome and one treatment variable")
      }
      parsed <- list(outcome = outcome, treatment = treatment, predictions = NULL)
    }

    # Create msd_data from formula
    if (inherits(data, "msd_data")) {
      msd <- data
    } else if (!is.null(data) && is.data.frame(data)) {
      msd <- msd_data(data = data, outcome = parsed$outcome, treatment = parsed$treatment)
    } else if (!is.null(observed)) {
      msd <- msd_data(observed = observed, unobserved = unobserved,
                      outcome = parsed$outcome, treatment = parsed$treatment)
    } else {
      stop("When using formula, must provide data or observed dataframe")
    }
  } else if (inherits(formula_or_data, "msd_data")) {
    msd <- formula_or_data
  } else {
    stop("First argument must be an msd_data object or a formula")
  }

  # Extract observed data
  obs <- msd$observed

  # Split by treatment arm
  Y1 <- obs$Y[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]

  n1 <- length(Y1)
  n0 <- length(Y0)

  # Compute point estimate
  mean_Y1 <- mean(Y1)
  mean_Y0 <- mean(Y0)
  estimate <- mean_Y1 - mean_Y0

  # Compute variance (using sample variance with n-1 denominator)
  var_Y1 <- var(Y1)
  var_Y0 <- var(Y0)
  variance <- var_Y1 / n1 + var_Y0 / n0

  # Standard error
  se <- sqrt(variance)

  # Confidence interval
  ci <- compute_ci(estimate, se, conf_level)

  # Return result
  new_msd_result(
    estimate = estimate,
    variance = variance,
    se = se,
    ci_lower = ci["ci_lower"],
    ci_upper = ci["ci_upper"],
    method = "Difference-in-Means (DiM)",
    lambda = NULL,
    n1 = n1,
    n0 = n0,
    m = 0,  # DiM doesn't use unlabeled data
    conf_level = conf_level
  )
}

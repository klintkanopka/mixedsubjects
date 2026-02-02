#' @title DiP (Difference-in-Predictions) Estimator
#' @description DiP estimator for ATE using paired predictions.

#' DiP (Difference-in-Predictions) Estimator
#'
#' Computes the Difference-in-Predictions (DiP) estimator for the average
#' treatment effect (ATE). This estimator uses both treatment and control
#' predictions for each unlabeled unit, computing the contrast S^(1) - S^(0)
#' at the unit level.
#'
#' @param formula_or_data Either an msd_data object created by \code{\link{msd_data}},
#'   or a formula of the form \code{outcome ~ treatment | pred_treated + pred_control}.
#' @param data If \code{formula_or_data} is a formula, this should be either:
#'   an msd_data object, a combined dataframe, or NULL (if using observed/unobserved).
#' @param observed If using formula with separate dataframes, the observed data.
#' @param unobserved If using formula with separate dataframes, the unobserved data.
#' @param conf_level Confidence level for the confidence interval (default 0.95)
#'
#' @return An msd_result object containing:
#'   \item{estimate}{Point estimate of the ATE}
#'   \item{variance}{Estimated variance}
#'   \item{se}{Standard error}
#'   \item{ci_lower, ci_upper}{Confidence interval bounds}
#'   \item{method}{Name of the estimation method}
#'   \item{lambda}{Tuning parameter (always 1 for DiP)}
#'
#' @details
#' The DiP estimator is:
#' \deqn{\hat{\tau}^{DiP} = \frac{1}{|\mathcal{U}|}\sum_{i \in \mathcal{U}}
#'   (S_i^{(1)} - S_i^{(0)}) + \frac{1}{n_1}\sum_{i \in \mathcal{O}_1}
#'   (Y_i - S_i^{(1)}) - \frac{1}{n_0}\sum_{i \in \mathcal{O}_0}(Y_i - S_i^{(0)})}
#'
#' @note
#' DiP requires BOTH S0 and S1 predictions for ALL units.
#' The key advantage of DiP over GREG is that when S^(1) and S^(0) are positively
#' correlated, the variance of their difference is smaller.
#'
#' @examples
#' # Using msd_data object
#' obs_df <- data.frame(
#'   Y = c(1.2, 1.4, 0.8, 0.6),
#'   S0 = c(1.0, 1.2, 0.7, 0.5),
#'   S1 = c(1.1, 1.3, 0.9, 0.7),
#'   D = c(1, 1, 0, 0)
#' )
#' unobs_df <- data.frame(
#'   S0 = c(1.1, 0.9, 1.0, 0.8),
#'   S1 = c(1.2, 1.0, 1.1, 0.9),
#'   D = c(1, 1, 0, 0)
#' )
#' msd <- msd_data(observed = obs_df, unobserved = unobs_df)
#' result <- msd_dip(msd)
#'
#' # Using formula interface
#' result2 <- msd_dip(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
#'
#' @export
msd_dip <- function(formula_or_data, data = NULL, observed = NULL,
                    unobserved = NULL, conf_level = 0.95) {

  # Resolve flexible input to msd_data
  msd <- resolve_msd_data(formula_or_data, data, observed, unobserved)

  # Check for both predictions
  if (!msd$has_both_predictions) {
    stop("DiP requires both S0 and S1 predictions for all units.")
  }

  if (is.null(msd$unobserved) || msd$m == 0) {
    stop("DiP requires unlabeled data. Use msd_dim() for labeled-only estimation.")
  }

  # Extract data
  obs <- msd$observed
  unobs <- msd$unobserved

  # Split observed by treatment arm
  Y1 <- obs$Y[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]
  n1 <- length(Y1)
  n0 <- length(Y0)

  # Get predictions on observed units (arm-appropriate)
  S1_obs_1 <- obs$S1[obs$D == 1]
  S0_obs_0 <- obs$S0[obs$D == 0]

  # Get both predictions on ALL unlabeled units
  S1_unobs <- unobs$S1
  S0_unobs <- unobs$S0
  m <- length(S1_unobs)

  # Compute DiP estimate
  mean_S_diff_unobs <- mean(S1_unobs - S0_unobs)
  mean_resid_1 <- mean(Y1 - S1_obs_1)
  mean_resid_0 <- mean(Y0 - S0_obs_0)

  estimate <- mean_S_diff_unobs + mean_resid_1 - mean_resid_0

  # Compute variance
  var_S_diff_unobs <- var(S1_unobs - S0_unobs)
  var_resid_1 <- var(Y1 - S1_obs_1)
  var_resid_0 <- var(Y0 - S0_obs_0)

  variance <- var_S_diff_unobs / m + var_resid_1 / n1 + var_resid_0 / n0

  se <- sqrt(variance)
  ci <- compute_ci(estimate, se, conf_level)

  new_msd_result(
    estimate = estimate,
    variance = variance,
    se = se,
    ci_lower = ci["ci_lower"],
    ci_upper = ci["ci_upper"],
    method = "DiP (Difference-in-Predictions, lambda = 1)",
    lambda = 1,
    n1 = n1,
    n0 = n0,
    m = m,
    conf_level = conf_level,
    additional = list(
      cov_S1_S0 = cov(S1_unobs, S0_unobs),
      var_S1 = var(S1_unobs),
      var_S0 = var(S0_unobs)
    )
  )
}

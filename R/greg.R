#' @title GREG Estimator
#' @description Generalized Regression (GREG) calibration estimator for ATE.

#' GREG Estimator
#'
#' Computes the Generalized Regression (GREG) calibration estimator for the
#' average treatment effect (ATE). This estimator corresponds to PPI with
#' tuning parameter lambda = 1.
#'
#' @param formula_or_data Either an msd_data object created by \code{\link{msd_data}},
#'   or a formula of the form \code{outcome ~ treatment | prediction}.
#'   For GREG, the formula specifies which prediction column(s) to use.
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
#'   \item{lambda}{Tuning parameter (always 1 for GREG)}
#'
#' @details
#' The GREG estimator for arm \eqn{d} is:
#' \deqn{\hat{\mu}_d^{GREG} = \bar{S}^{(d)}_{\mathcal{U}_d} +
#'   (\bar{Y}_{\mathcal{O}_d} - \bar{S}^{(d)}_{\mathcal{O}_d})}
#'
#' The ATE estimate is:
#' \deqn{\hat{\tau}^{GREG} = \hat{\mu}_1^{GREG} - \hat{\mu}_0^{GREG}}
#'
#' The variance is:
#' \deqn{\widehat{Var}(\hat{\tau}^{GREG}) =
#'   \sum_{d \in \{0,1\}} \left[\frac{s^2_{S^{(d)}}}{m_d} +
#'   \frac{Var(Y(d) - S^{(d)})}{n_d}\right]}
#'
#' @note
#' GREG requires predictions for each unit's assigned arm:
#' - Treatment arm units need S1
#' - Control arm units need S0
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
#' result <- msd_greg(msd)
#'
#' # Using formula interface
#' result2 <- msd_greg(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
#'
#' @export
msd_greg <- function(formula_or_data, data = NULL, observed = NULL,
                     unobserved = NULL, conf_level = 0.95) {

  # Resolve flexible input to msd_data
  msd <- resolve_msd_data(formula_or_data, data, observed, unobserved)

  # Check for predictions
  if (!msd$has_S0 && !msd$has_S1) {
    stop("GREG requires predictions. Neither S0 nor S1 found in data.")
  }

  if (is.null(msd$unobserved) || msd$m == 0) {
    stop("GREG requires unlabeled data. Use msd_dim() for labeled-only estimation.")
  }

  # Extract data
  obs <- msd$observed
  unobs <- msd$unobserved

  # Split by treatment arm
  Y1 <- obs$Y[obs$D == 1]
  Y0 <- obs$Y[obs$D == 0]
  n1 <- length(Y1)
  n0 <- length(Y0)

  # Get arm-specific predictions on observed
  S1_obs_1 <- obs$S1[obs$D == 1]
  S0_obs_0 <- obs$S0[obs$D == 0]

  # Get arm-specific predictions on unobserved
  S1_unobs_1 <- unobs$S1[unobs$D == 1]
  S0_unobs_0 <- unobs$S0[unobs$D == 0]
  m1 <- length(S1_unobs_1)
  m0 <- length(S0_unobs_0)

  if (m1 == 0 || m0 == 0) {
    stop("GREG requires unlabeled units in both treatment arms.")
  }

  # Compute GREG estimate for each arm
  mu_1 <- mean(S1_unobs_1) + (mean(Y1) - mean(S1_obs_1))
  mu_0 <- mean(S0_unobs_0) + (mean(Y0) - mean(S0_obs_0))

  estimate <- mu_1 - mu_0

  # Compute variance
  var_S1_unobs <- var(S1_unobs_1)
  var_S0_unobs <- var(S0_unobs_0)

  resid_1 <- Y1 - S1_obs_1
  resid_0 <- Y0 - S0_obs_0
  var_resid_1 <- var(resid_1)
  var_resid_0 <- var(resid_0)

  variance <- var_S1_unobs / m1 + var_S0_unobs / m0 +
              var_resid_1 / n1 + var_resid_0 / n0

  se <- sqrt(variance)
  ci <- compute_ci(estimate, se, conf_level)

  new_msd_result(
    estimate = estimate,
    variance = variance,
    se = se,
    ci_lower = ci["ci_lower"],
    ci_upper = ci["ci_upper"],
    method = "GREG (lambda = 1)",
    lambda = 1,
    n1 = n1,
    n0 = n0,
    m = m1 + m0,
    conf_level = conf_level,
    additional = list(m1 = m1, m0 = m0)
  )
}

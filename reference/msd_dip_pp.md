# DiP++ Estimator

Power-tuned Difference-in-Predictions estimator for ATE. DiP++ Estimator

Computes the DiP++ (power-tuned difference-in-predictions) estimator for
the average treatment effect (ATE). This estimator uses paired
predictions S^(1) and S^(0) for each unlabeled unit, with a single
tuning parameter lambda estimated via cross-fitting.

## Usage

``` r
msd_dip_pp(
  formula_or_data,
  data = NULL,
  observed = NULL,
  unobserved = NULL,
  n_folds = 2,
  conf_level = 0.95,
  seed = NULL
)
```

## Arguments

- formula_or_data:

  Either an msd_data object created by
  [`msd_data`](http://klintkanopka.com/mixedsubjects/reference/msd_data.md),
  or a formula of the form
  `outcome ~ treatment | pred_treated + pred_control`.

- data:

  If `formula_or_data` is a formula, this should be either: an msd_data
  object, a combined dataframe, or NULL (if using observed/unobserved).

- observed:

  If using formula with separate dataframes, the observed data.

- unobserved:

  If using formula with separate dataframes, the unobserved data.

- n_folds:

  Number of folds for cross-fitting (default 2)

- conf_level:

  Confidence level for the confidence interval (default 0.95)

- seed:

  Random seed for fold splitting (optional)

## Value

An msd_result object containing:

- estimate:

  Point estimate of the ATE

- variance:

  Estimated variance of the cross-fit estimator (see the Variance
  section)

- se:

  Standard error

- ci_lower, ci_upper:

  Confidence interval bounds

- method:

  Name of the estimation method

- lambda:

  Estimated tuning parameter (single value)

## Details

The DiP++ estimator is: \$\$\hat{\tau}^{DiP++}(\lambda) =
\frac{\lambda}{\|\mathcal{U}\|} \sum\_{i \in \mathcal{U}} (S_i^{(1)} -
S_i^{(0)}) + \frac{1}{n_1}\sum\_{i \in \mathcal{O}\_1}(Y_i - \lambda
S_i^{(1)}) - \frac{1}{n_0}\sum\_{i \in \mathcal{O}\_0}(Y_i - \lambda
S_i^{(0)})\$\$

## Note

DiP++ requires BOTH S0 and S1 predictions for ALL units. For
arm-specific tuning, use
[`msd_dt_dip`](http://klintkanopka.com/mixedsubjects/reference/msd_dt_dip.md).

## Variance

The reported variance is an analytic estimate of the finite-sample
variance of the cross-fit estimator, which averages K fold estimates
that share the unlabeled pool. It has three parts: (i) the variance of
each fold estimate, evaluated at the tuning actually used in that fold
with that fold's own sample moments (these are independent of the
tuning, so estimation of the tuning is accounted for without
approximation); (ii) the covariance between fold estimates induced by
the shared unlabeled pool; and (iii) a second-order cross-fold term that
arises because the tuning used in one fold is estimated on the data that
the other fold is evaluated on. Part (iii) is the product of the
covariances between each tuning estimate and the labeled prediction
means of the fold it was estimated on; it is estimated in closed form
from third central moments, is zero when outcomes and predictions are
jointly symmetric (for example Gaussian), and is of order \\n^{-2}\\.
The components are returned in `additional$variance_components`. The
fixed-tuning plug-in used in version 1.0.0, which omits the variance of
the estimated tuning and part (iii), is returned in
`additional$variance_fixed_lambda`.

## Tuning

The tuning parameter used in each fold is estimated from the other
folds' labeled data only. The unlabeled pool enters the tuning rule only
through its size. In version 1.0.0 the variance of \\S^{(1)} - S^{(0)}\\
in the tuning rule was estimated from the unlabeled pool; it is now
estimated from the labeled tuning units, so that the tuning is
independent of the data it is applied to and the cross-fit estimator is
unbiased in finite samples.

## Examples

``` r
# Using msd_data object
set.seed(123)
n <- 100
obs_df <- data.frame(
  Y = rnorm(n),
  D = rep(c(1, 0), each = n/2)
)
obs_df$Y <- obs_df$Y + 0.3 * obs_df$D
obs_df$S1 <- 0.5 * obs_df$Y + rnorm(n, 0, 0.5)
obs_df$S0 <- 0.5 * obs_df$Y + rnorm(n, 0, 0.5) - 0.1

unobs_df <- data.frame(
  S0 = rnorm(200, 0, 0.5),
  S1 = rnorm(200, 0.2, 0.5),
  D = rep(c(1, 0), each = 100)
)

msd <- msd_data(observed = obs_df, unobserved = unobs_df)
result <- msd_dip_pp(msd)

# Using formula interface
result2 <- msd_dip_pp(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
```

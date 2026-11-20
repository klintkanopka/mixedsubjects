# D-T DiP (Doubly-Tuned Difference-in-Predictions) Estimator

D-T DiP estimator with arm-specific tuning for ATE. D-T DiP
(Doubly-Tuned Difference-in-Predictions) Estimator

Computes the D-T DiP estimator for the average treatment effect (ATE).
This estimator uses paired predictions S^(1) and S^(0) for each
unlabeled unit, with arm-specific tuning parameters (lambda_1 and
lambda_0) estimated via cross-fitting.

## Usage

``` r
msd_dt_dip(
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

  Vector of arm-specific tuning parameters (lambda_1, lambda_0)

## Details

The D-T DiP estimator is: \$\$\hat{\tau}^{D-T DiP} =
\frac{1}{\|\mathcal{U}\|} \sum\_{i \in \mathcal{U}} (\lambda_1
S_i^{(1)} - \lambda_0 S_i^{(0)}) + \frac{1}{n_1}\sum\_{i \in
\mathcal{O}\_1}(Y_i - \lambda_1 S_i^{(1)}) - \frac{1}{n_0}\sum\_{i \in
\mathcal{O}\_0}(Y_i - \lambda_0 S_i^{(0)})\$\$

The tuning parameters \\(\lambda_1, \lambda_0)\\ are chosen jointly to
minimize the total variance. Because the two arms share the unobserved
pool, the variance is jointly quadratic in \\(\lambda_1, \lambda_0)\\
with a cross-term from \\Cov(S^{(1)}, S^{(0)})\\, so the optimum solves
a coupled 2x2 linear system rather than two independent regressions.

The tuning parameters are estimated via cross-fitting:

1.  Split labeled data into K folds

2.  For each fold k, estimate (lambda_1, lambda_0) on opposite folds

3.  Compute the fold-k estimate using estimated lambdas

4.  Average across folds with equal weights

## Note

D-T DiP requires BOTH S0 and S1 predictions for ALL units. This
corresponds to 2 predictions per unlabeled unit.

D-T DiP combines the benefits of:

- DiP: exploiting positive correlation between S^(1) and S^(0)

- D-T: arm-specific tuning for heterogeneous prediction quality

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

The tuning parameters used in each fold are estimated from the other
folds' labeled data only. The unlabeled pool enters the tuning rule only
through its size. In version 1.0.0 the variances of \\S^{(1)}\\ and
\\S^{(0)}\\ and their covariance in the terms weighted by \\1/m\\ were
estimated from the unlabeled pool; they are now estimated from the
labeled tuning units, so that the tunings are independent of the data
they are applied to and the cross-fit estimator is unbiased in finite
samples.

## Examples

``` r
# Create sample data with both predictions
set.seed(123)
n <- 100
obs_df <- data.frame(
  Y = rnorm(n),
  D = rep(c(1, 0), each = n/2)
)
obs_df$Y <- obs_df$Y + 0.3 * obs_df$D
obs_df$S1 <- 0.6 * obs_df$Y + rnorm(n, 0, 0.4)
obs_df$S0 <- 0.4 * obs_df$Y + rnorm(n, 0, 0.5)

unobs_df <- data.frame(
  S0 = rnorm(300, 0, 0.5),
  S1 = rnorm(300, 0.2, 0.4),
  D = rep(c(1, 0), 150)
)
# Add correlation between S0 and S1
unobs_df$S1 <- unobs_df$S1 + 0.5 * unobs_df$S0

msd <- msd_data(observed = obs_df, unobserved = unobs_df)
result <- msd_dt_dip(msd)

# Using formula interface
result2 <- msd_dt_dip(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
```

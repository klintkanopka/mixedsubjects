# D-T (Doubly-Tuned) Estimator

D-T estimator with arm-specific tuning parameters for ATE. D-T
(Doubly-Tuned) Estimator

Computes the Doubly-Tuned (D-T) estimator for the average treatment
effect (ATE). This estimator uses arm-specific tuning parameters
(lambda_1 and lambda_0) estimated via cross-fitting.

## Usage

``` r
msd_dt(
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
  or a formula of the form `outcome ~ treatment | prediction`.

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

The D-T estimator uses arm-specific tuning parameters:
\$\$\hat{\mu}\_d^{D-T}(\lambda_d) = \bar{Y}\_{\mathcal{O}\_d} +
\lambda_d(\bar{S}^{(d)}\_{\mathcal{U}\_d} -
\bar{S}^{(d)}\_{\mathcal{O}\_d})\$\$

Each lambda_d is chosen to minimize the variance in arm d:
\$\$\lambda_d^\* = \frac{Cov(Y(d), S^{(d)}) / n_d}{Var(S^{(d)})(1/m_d +
1/n_d)}\$\$

With cross-fitting, the arm mean is computed on each fold using a tuning
estimated on the other folds' labeled data, and the fold-specific arm
means are combined with weights equal to each fold's share of the arm's
labeled units, \\\pi_d^{(k)} = n_d^{(k)} / n_d\\. The weights are
returned in `additional$fold_weights`. They equal \\1/K\\ whenever an
arm's labeled sample size is divisible by the number of folds.

The tuning parameters are estimated via cross-fitting to avoid bias.

## Note

D-T differs from PPI++ by using separate tuning parameters for each arm,
which can improve efficiency when the prediction quality differs between
treatment and control.

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

## Examples

``` r
# Create sample data
set.seed(123)
n <- 100
obs_df <- data.frame(
  Y = rnorm(n),
  S0 = rnorm(n, 0, 0.5),
  S1 = rnorm(n, 0.2, 0.5),
  D = rep(c(1, 0), each = n/2)
)
obs_df$Y <- obs_df$Y + 0.3 * obs_df$D
obs_df$S1[obs_df$D == 1] <- obs_df$S1[obs_df$D == 1] + 0.5 * obs_df$Y[obs_df$D == 1]
obs_df$S0[obs_df$D == 0] <- obs_df$S0[obs_df$D == 0] + 0.5 * obs_df$Y[obs_df$D == 0]

unobs_df <- data.frame(
  S0 = rnorm(200, 0, 0.5),
  S1 = rnorm(200, 0.2, 0.5),
  D = rep(c(1, 0), each = 100)
)

msd <- msd_data(observed = obs_df, unobserved = unobs_df)
result <- msd_dt(msd)

# Using formula interface
result2 <- msd_dt(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
```

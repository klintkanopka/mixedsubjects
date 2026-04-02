# PPI++ Estimator

Power-tuned Prediction-Powered Inference estimator for ATE. PPI++
Estimator

Computes the PPI++ (power-tuned prediction-powered inference) estimator
for the average treatment effect (ATE). This estimator uses a single
tuning parameter lambda that is estimated via cross-fitting to minimize
variance.

## Usage

``` r
msd_ppi(
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

  Estimated variance (delta-method)

- se:

  Standard error

- ci_lower, ci_upper:

  Confidence interval bounds

- method:

  Name of the estimation method

- lambda:

  Estimated tuning parameter (single value)

## Details

The PPI++ estimator uses a single tuning parameter lambda across both
arms: \$\$\hat{\mu}\_d^{PPI}(\lambda) = \bar{Y}\_{\mathcal{O}\_d} +
\lambda(\bar{S}^{(d)}\_{\mathcal{U}\_d} -
\bar{S}^{(d)}\_{\mathcal{O}\_d})\$\$

The tuning parameter is estimated via cross-fitting:

1.  Split labeled data into K folds

2.  For each fold k, estimate lambda on the opposite folds

3.  Compute the fold-k estimate using the estimated lambda

4.  Average across folds

Lambda is chosen to minimize the combined variance across arms.

## Note

PPI++ uses a single lambda across arms, unlike D-T which uses
arm-specific tuning parameters. For arm-specific tuning, use
[`msd_dt`](http://klintkanopka.com/mixedsubjects/reference/msd_dt.md).

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
result <- msd_ppi(msd)

# Using formula interface
result2 <- msd_ppi(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
```

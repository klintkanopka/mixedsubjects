# Compare variance estimates across methods

Computes and compares variance estimates using both delta-method and
bootstrap for a given estimator.

## Usage

``` r
compare_variance_methods(
  data,
  estimator,
  n_bootstrap = 1000,
  n_folds = 2,
  seed = NULL
)
```

## Arguments

- data:

  An msd_data object

- estimator:

  Character string specifying the estimator

- n_bootstrap:

  Number of bootstrap replications

- n_folds:

  Number of folds for cross-fitting

- seed:

  Random seed

## Value

A data frame comparing variance estimates

## Examples

``` r
obs_df <- data.frame(
  Y = rnorm(100), S0 = rnorm(100), S1 = rnorm(100),
  D = rep(c(1, 0), each = 50)
)
unobs_df <- data.frame(
  S0 = rnorm(200), S1 = rnorm(200), D = rep(c(1, 0), each = 100)
)
msd <- msd_data(observed = obs_df, unobserved = unobs_df)
comparison <- compare_variance_methods(msd, "dt_dip", n_bootstrap = 100, seed = 1)
print(comparison)
#>                Method   Estimate   Variance        SE   CI_Lower  CI_Upper
#> ci_lower Delta-method 0.08423398 0.04932197 0.2220855 -0.3510456 0.5195136
#>             Bootstrap 0.10638804 0.05242130 0.2289570 -0.4321612 0.4863055
```

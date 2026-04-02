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
if (FALSE) { # \dontrun{
msd <- msd_data(observed = obs_df, unobserved = unobs_df)
comparison <- compare_variance_methods(msd, "dt_dip", n_bootstrap = 500)
print(comparison)
} # }
```

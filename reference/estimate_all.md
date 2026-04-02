# Estimate all available estimators

Runs all applicable estimators on the data and returns a summary table.

## Usage

``` r
estimate_all(data, n_folds = 2, conf_level = 0.95)
```

## Arguments

- data:

  An msd_data object

- n_folds:

  Number of folds for cross-fitting (default 2)

- conf_level:

  Confidence level (default 0.95)

## Value

A data frame with estimates from all applicable estimators

## Examples

``` r
if (FALSE) { # \dontrun{
msd <- msd_data(observed = obs_df, unobserved = unobs_df)
all_estimates <- estimate_all(msd)
print(all_estimates)
} # }
```

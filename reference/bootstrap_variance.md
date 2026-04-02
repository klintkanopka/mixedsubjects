# Variance Estimation for Mixed-Subjects Design

Variance estimation methods including fold-respecting bootstrap.
Bootstrap Variance Estimation

Computes bootstrap variance estimates for MSD estimators using a
fold-respecting resampling procedure.

## Usage

``` r
bootstrap_variance(
  data,
  estimator = c("dim", "greg", "ppi", "dt", "dip", "dip_pp", "dt_dip"),
  n_bootstrap = 1000,
  n_folds = 2,
  conf_level = 0.95,
  seed = NULL
)
```

## Arguments

- data:

  An msd_data object created by
  [`msd_data`](http://klintkanopka.com/mixedsubjects/reference/msd_data.md)

- estimator:

  Character string specifying the estimator. One of: "dim", "greg",
  "ppi", "dt", "dip", "dip_pp", "dt_dip"

- n_bootstrap:

  Number of bootstrap replications (default 1000)

- n_folds:

  Number of folds for cross-fitting (default 2, used for tuned
  estimators)

- conf_level:

  Confidence level for confidence intervals (default 0.95)

- seed:

  Random seed for reproducibility (optional)

## Value

A list containing:

- estimate:

  Point estimate from the original data

- variance:

  Bootstrap variance estimate

- se:

  Bootstrap standard error

- ci_lower, ci_upper:

  Bootstrap percentile confidence interval

- bootstrap_estimates:

  Vector of bootstrap estimates

## Details

The fold-respecting bootstrap resamples within each stratum:

1.  Resample observed treatment units with replacement

2.  Resample observed control units with replacement

3.  Resample unlabeled units with replacement

4.  Recompute the estimator on the bootstrap sample

For cross-fit estimators, the fold assignments are regenerated for each
bootstrap replicate to properly account for the cross-fitting variance.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create sample data
obs_df <- data.frame(
  Y = rnorm(100),
  S0 = rnorm(100),
  S1 = rnorm(100),
  D = rep(c(1, 0), each = 50)
)
unobs_df <- data.frame(
  S0 = rnorm(200),
  S1 = rnorm(200),
  D = rep(c(1, 0), each = 100)
)
msd <- msd_data(observed = obs_df, unobserved = unobs_df)

# Bootstrap variance for D-T DiP
boot_result <- bootstrap_variance(msd, "dt_dip", n_bootstrap = 500)
print(boot_result)
} # }
```

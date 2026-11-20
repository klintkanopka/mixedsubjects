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
boot_result <- bootstrap_variance(msd, "dt_dip", n_bootstrap = 100, seed = 1)
print(boot_result)
#> $estimate
#> [1] -0.3205735
#> 
#> $variance
#> [1] 0.0450314
#> 
#> $se
#> [1] 0.212206
#> 
#> $ci_lower
#> [1] -0.7331161
#> 
#> $ci_upper
#> [1] 0.07166747
#> 
#> $conf_level
#> [1] 0.95
#> 
#> $n_bootstrap
#> [1] 100
#> 
#> $bootstrap_estimates
#>   [1] -0.26974974 -0.51555599 -0.16670351 -0.10329068  0.03042996 -0.18644122
#>   [7] -0.52093397  0.09144463 -0.67565865 -0.13329390 -0.43757263  0.06926163
#>  [13] -0.11073444 -0.34703979 -0.07662014  0.02124554 -0.38039215 -0.49318962
#>  [19] -0.47052645 -0.42768691 -0.28051531 -0.73590112 -0.67878794 -0.56955837
#>  [25] -0.27157715  0.05440719 -0.35939225 -0.26780270 -0.29044425 -0.27591813
#>  [31] -0.23662250 -0.08566847 -0.13298016 -0.11387436 -0.16917846 -0.41117010
#>  [37] -0.45648809 -0.34007279 -0.37881793 -0.06751452 -0.13420662 -0.43997985
#>  [43] -0.09967829 -0.27398678 -0.73062341 -0.11938959 -0.77258512 -0.20797005
#>  [49] -0.11464418 -0.27180797 -0.39865524  0.31992314 -0.16762088  0.02028708
#>  [55] -0.65335807 -0.45158252 -0.36783898  0.02388423 -0.44217996 -0.46486651
#>  [61] -0.57766078 -0.73537144 -0.32790729 -0.16072534 -0.37886494 -0.47824605
#>  [67] -0.02275085 -0.51306156 -0.47745925 -0.47971084 -0.38628303 -0.27140255
#>  [73] -0.38621360 -0.56695369 -0.69966970 -0.18620639 -0.47024684 -0.24315918
#>  [79] -0.49307209 -0.29364961 -0.38658358 -0.22795715 -0.29954021 -0.28927797
#>  [85] -0.01873075 -0.18043787 -0.19423600 -0.28603505 -0.30159393 -0.32344633
#>  [91] -0.25286663 -0.52324312 -0.33619066 -0.22688892  0.07384418 -0.26359122
#>  [97] -0.55667910 -0.15548253 -0.18394414 -0.40380436
#> 
#> $method
#> [1] "Bootstrap (dt_dip, B=100)"
#> 
```

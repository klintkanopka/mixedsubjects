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
#> [1] -0.3206185
#> 
#> $variance
#> [1] 0.04494434
#> 
#> $se
#> [1] 0.2120008
#> 
#> $ci_lower
#> [1] -0.733097
#> 
#> $ci_upper
#> [1] 0.07037486
#> 
#> $conf_level
#> [1] 0.95
#> 
#> $n_bootstrap
#> [1] 100
#> 
#> $bootstrap_estimates
#>   [1] -0.26864458 -0.51661552 -0.16763678 -0.10331910  0.03072156 -0.18608884
#>   [7] -0.52090544  0.08423339 -0.67677209 -0.13195015 -0.43241436  0.06780552
#>  [13] -0.11116158 -0.34739121 -0.07396856  0.01931440 -0.37942003 -0.49199017
#>  [19] -0.47144105 -0.42765536 -0.28092503 -0.73515731 -0.67905354 -0.57027202
#>  [25] -0.27056301  0.05520295 -0.35881491 -0.26885930 -0.29093297 -0.27552365
#>  [31] -0.23609528 -0.08552076 -0.13251844 -0.11429088 -0.16955814 -0.41434147
#>  [37] -0.45691441 -0.34026256 -0.37859883 -0.06786513 -0.13173245 -0.43989278
#>  [43] -0.10603069 -0.27364173 -0.73094638 -0.11804863 -0.77143346 -0.20620451
#>  [49] -0.11470766 -0.27140146 -0.39851113  0.32142884 -0.16758444  0.02240419
#>  [55] -0.65452473 -0.45123902 -0.36931356  0.02506665 -0.44119780 -0.46298991
#>  [61] -0.57797188 -0.73504290 -0.32690784 -0.15940339 -0.38062030 -0.47745875
#>  [67] -0.02322449 -0.51049336 -0.47673005 -0.47727368 -0.38570210 -0.27049797
#>  [73] -0.38679571 -0.56706378 -0.70029967 -0.18618822 -0.47133639 -0.24367413
#>  [79] -0.49443146 -0.29371591 -0.38670928 -0.22917800 -0.30047723 -0.28998543
#>  [85] -0.02183866 -0.17692518 -0.19335651 -0.28610907 -0.30189490 -0.32361728
#>  [91] -0.25218618 -0.51720075 -0.33569546 -0.22385134  0.07269950 -0.26246726
#>  [97] -0.55653691 -0.15879304 -0.18290196 -0.40370967
#> 
#> $method
#> [1] "Bootstrap (dt_dip, B=100)"
#> 
```

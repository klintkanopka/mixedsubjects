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
obs_df <- data.frame(
  Y = rnorm(100), S0 = rnorm(100), S1 = rnorm(100),
  D = rep(c(1, 0), each = 50)
)
unobs_df <- data.frame(
  S0 = rnorm(200), S1 = rnorm(200), D = rep(c(1, 0), each = 100)
)
msd <- msd_data(observed = obs_df, unobserved = unobs_df)
all_estimates <- estimate_all(msd)
print(all_estimates)
#> 
#> Mixed-Subjects Design: All Estimators
#> ======================================
#> 
#>                                    Estimator Estimate     SE 95% CI Lower
#>                    Difference-in-Means (DiM)   0.0854 0.2197      -0.3451
#>                            GREG (lambda = 1)   0.2438 0.3104      -0.3646
#>                       PPI++ (cross-fit, K=2)   0.0874 0.2197      -0.3432
#>           D-T (Doubly-Tuned, cross-fit, K=2)   0.0709 0.2197      -0.3597
#>  DiP (Difference-in-Predictions, lambda = 1)   0.3167 0.2953      -0.2621
#>                       DiP++ (cross-fit, K=2)   0.0840 0.2197      -0.3466
#>                     D-T DiP (cross-fit, K=2)   0.0709 0.2195      -0.3594
#>  95% CI Upper
#>        0.5160
#>        0.8522
#>        0.5179
#>        0.5016
#>        0.8955
#>        0.5146
#>        0.5012
#> 
```

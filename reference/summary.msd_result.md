# Summary method for msd_result

Produces a detailed summary of the MSD estimation results, including the
coefficient table with z-statistics and p-values, sample size
information, and interpretation guidance.

## Usage

``` r
# S3 method for class 'msd_result'
summary(object, ...)
```

## Arguments

- object:

  An msd_result object

- ...:

  Additional arguments (ignored)

## Value

A summary.msd_result object (invisibly when printed)

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
result <- msd_dt_dip(msd, seed = 1)
summary(result)
#> 
#> Mixed-Subjects Design: Treatment Effect Estimation
#> ===================================================
#> 
#> Estimator: D-T DiP (cross-fit, K=2) 
#> 
#> Coefficients:
#>     Estimate Std. Error z value Pr(>|z|)  
#> ATE  -0.0704     0.1940  -0.363   0.7166  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> 95%Confidence Interval: [-0.4507, 0.3098]
#> 
#> Tuning Parameters:
#>   lambda_1 (treatment): -0.0902 
#>   lambda_0 (control):   -0.0019 
#> 
#> Sample Sizes:
#>   Observed (human subjects): 100
#>     Treated (n1):  50
#>     Control (n0):  50
#>   Unobserved (predictions only): 200
#>   Total units: 300
#> 
#> Interpretation:
#>   The estimated effect is 0.0704 units, but this is not statistically significant (p = 0.717).
#>   Using 200 unlabeled predictions improved estimation precision.
#> 
```

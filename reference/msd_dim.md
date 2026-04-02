# Difference-in-Means Estimator

Classical difference-in-means estimator for ATE. Difference-in-Means
Estimator

Computes the classical difference-in-means estimator for the average
treatment effect (ATE). This estimator uses only observed (labeled) data
and does not incorporate any predictions.

## Usage

``` r
msd_dim(
  formula_or_data,
  data = NULL,
  observed = NULL,
  unobserved = NULL,
  conf_level = 0.95
)
```

## Arguments

- formula_or_data:

  Either an msd_data object created by
  [`msd_data`](http://klintkanopka.com/mixedsubjects/reference/msd_data.md),
  or a formula of the form `outcome ~ treatment` (predictions not needed
  for DiM).

- data:

  If `formula_or_data` is a formula, this should be either: an msd_data
  object, a combined dataframe, or NULL (if using observed/unobserved).

- observed:

  If using formula with separate dataframes, the observed data.

- unobserved:

  If using formula with separate dataframes, the unobserved data.

- conf_level:

  Confidence level for the confidence interval (default 0.95)

## Value

An msd_result object containing:

- estimate:

  Point estimate of the ATE: mean(Y\|D=1) - mean(Y\|D=0)

- variance:

  Estimated variance: var(Y\|D=1)/n1 + var(Y\|D=0)/n0

- se:

  Standard error

- ci_lower, ci_upper:

  Confidence interval bounds

- method:

  Name of the estimation method

## Details

The difference-in-means estimator is: \$\$\hat{\tau}^{DiM} =
\bar{Y}\_1 - \bar{Y}\_0\$\$

where \\\bar{Y}\_d\\ is the sample mean of outcomes in arm \\d\\.

The variance is estimated as: \$\$\widehat{Var}(\hat{\tau}^{DiM}) =
\frac{s^2\_{Y(1)}}{n_1} + \frac{s^2\_{Y(0)}}{n_0}\$\$

where \\s^2\_{Y(d)}\\ is the sample variance of outcomes in arm \\d\\.

## Examples

``` r
# Using msd_data object (standard interface)
obs_df <- data.frame(
  Y = c(1.2, 1.4, 0.8, 0.6, 1.1, 0.9, 1.3, 0.7),
  S0 = c(1.0, 1.2, 0.7, 0.5, 1.0, 0.8, 1.1, 0.6),
  S1 = c(1.1, 1.3, 0.9, 0.7, 1.1, 0.9, 1.2, 0.8),
  D = c(1, 1, 0, 0, 1, 0, 1, 0)
)
msd <- msd_data(observed = obs_df)

result <- msd_dim(msd)
print(result)
#> 
#> Mixed-Subjects Design Estimation
#> =================================
#> Estimator: Difference-in-Means (DiM) 
#> 
#> Point Estimate:  0.5000 
#> Standard Error:  0.0913 
#> 95% CI:         [0.3211, 0.6789]
#> 
#> Sample Sizes:
#>   Observed:   n_1=4, n_0=4
#> 

# Using formula interface with custom column names
df <- data.frame(
  response = c(1.2, 1.4, 0.8, 0.6),
  treated = c(1, 1, 0, 0)
)
result2 <- msd_dim(response ~ treated, observed = df)
```

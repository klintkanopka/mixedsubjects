# GREG Estimator

Generalized Regression (GREG) calibration estimator for ATE. GREG
Estimator

Computes the Generalized Regression (GREG) calibration estimator for the
average treatment effect (ATE). This estimator corresponds to PPI with
tuning parameter lambda = 1.

## Usage

``` r
msd_greg(
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
  or a formula of the form `outcome ~ treatment | prediction`. For GREG,
  the formula specifies which prediction column(s) to use.

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

  Point estimate of the ATE

- variance:

  Estimated variance

- se:

  Standard error

- ci_lower, ci_upper:

  Confidence interval bounds

- method:

  Name of the estimation method

- lambda:

  Tuning parameter (always 1 for GREG)

## Details

The GREG estimator for arm \\d\\ is: \$\$\hat{\mu}\_d^{GREG} =
\bar{S}^{(d)}\_{\mathcal{U}\_d} + (\bar{Y}\_{\mathcal{O}\_d} -
\bar{S}^{(d)}\_{\mathcal{O}\_d})\$\$

The ATE estimate is: \$\$\hat{\tau}^{GREG} = \hat{\mu}\_1^{GREG} -
\hat{\mu}\_0^{GREG}\$\$

The variance is: \$\$\widehat{Var}(\hat{\tau}^{GREG}) = \sum\_{d \in
\\0,1\\} \left\[\frac{s^2\_{S^{(d)}}}{m_d} + \frac{Var(Y(d) -
S^{(d)})}{n_d}\right\]\$\$

## Note

GREG requires predictions for each unit's assigned arm:

- Treatment arm units need S1

- Control arm units need S0

## Examples

``` r
# Using msd_data object
obs_df <- data.frame(
  Y = c(1.2, 1.4, 0.8, 0.6),
  S0 = c(1.0, 1.2, 0.7, 0.5),
  S1 = c(1.1, 1.3, 0.9, 0.7),
  D = c(1, 1, 0, 0)
)
unobs_df <- data.frame(
  S0 = c(1.1, 0.9, 1.0, 0.8),
  S1 = c(1.2, 1.0, 1.1, 0.9),
  D = c(1, 1, 0, 0)
)
msd <- msd_data(observed = obs_df, unobserved = unobs_df)
result <- msd_greg(msd)

# Using formula interface
result2 <- msd_greg(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
```

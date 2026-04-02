# DiP (Difference-in-Predictions) Estimator

DiP estimator for ATE using paired predictions. DiP
(Difference-in-Predictions) Estimator

Computes the Difference-in-Predictions (DiP) estimator for the average
treatment effect (ATE). This estimator uses both treatment and control
predictions for each unlabeled unit, computing the contrast S^(1) -
S^(0) at the unit level.

## Usage

``` r
msd_dip(
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
  or a formula of the form
  `outcome ~ treatment | pred_treated + pred_control`.

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

  Tuning parameter (always 1 for DiP)

## Details

The DiP estimator is: \$\$\hat{\tau}^{DiP} =
\frac{1}{\|\mathcal{U}\|}\sum\_{i \in \mathcal{U}} (S_i^{(1)} -
S_i^{(0)}) + \frac{1}{n_1}\sum\_{i \in \mathcal{O}\_1} (Y_i -
S_i^{(1)}) - \frac{1}{n_0}\sum\_{i \in \mathcal{O}\_0}(Y_i -
S_i^{(0)})\$\$

## Note

DiP requires BOTH S0 and S1 predictions for ALL units. The key advantage
of DiP over GREG is that when S^(1) and S^(0) are positively correlated,
the variance of their difference is smaller.

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
result <- msd_dip(msd)

# Using formula interface
result2 <- msd_dip(Y ~ D | S1 + S0, observed = obs_df, unobserved = unobs_df)
```

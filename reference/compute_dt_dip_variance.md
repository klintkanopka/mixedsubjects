# Compute D-T DiP variance using delta method

The unobserved component var_U = Var(lambda1*S1 - lambda0*S0) / m is
shared across all folds and is NOT divided by K. The labeled components'
K factors cancel when averaging across folds.

## Usage

``` r
compute_dt_dip_variance(data, lambda1, lambda0, n_folds)
```

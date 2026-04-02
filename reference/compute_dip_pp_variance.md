# Compute DiP++ variance using delta method

The unobserved component var_U = lambda^2 \* Var(S1-S0) / m is shared
across all folds and is NOT divided by K. The labeled components' K
factors cancel when averaging across folds.

## Usage

``` r
compute_dip_pp_variance(data, lambda, n_folds)
```

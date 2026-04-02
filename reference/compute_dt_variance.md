# Compute D-T variance using delta method

See compute_ppi_variance for the derivation. The labeled terms' K
factors cancel when averaging across folds; the unobserved term is
shared across folds and is NOT divided by K.

## Usage

``` r
compute_dt_variance(data, lambda1, lambda0, n_folds)
```

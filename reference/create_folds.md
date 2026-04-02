# Create fold indices for cross-fitting

Create fold indices for cross-fitting

## Usage

``` r
create_folds(n, n_folds = 2, seed = NULL)
```

## Arguments

- n:

  Number of observations

- n_folds:

  Number of folds (default 2)

- seed:

  Random seed for reproducibility (optional)

## Value

A vector of fold assignments (integers 1 to n_folds)

# Validate that a column exists in a dataframe

Validate that a column exists in a dataframe

## Usage

``` r
validate_column(df, col_name, type, required = TRUE)
```

## Arguments

- df:

  Dataframe to check

- col_name:

  Column name to validate

- type:

  Type of column (for error messages): "outcome", "treatment",
  "pred_control", "pred_treated"

- required:

  If TRUE, error when not found; if FALSE, return NULL

## Value

The column name if found, or NULL if not required and not found

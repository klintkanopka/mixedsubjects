# Resolve data from formula or msd_data object

Handles the flexible interface where users can provide either:

1.  An msd_data object directly

2.  A formula + raw dataframe(s)

## Usage

``` r
resolve_msd_data(
  formula_or_data,
  data = NULL,
  observed = NULL,
  unobserved = NULL,
  require_predictions = TRUE
)
```

## Arguments

- formula_or_data:

  Either a formula or an msd_data object

- data:

  If formula provided, this should be the data (msd_data, combined df,
  or NULL)

- observed:

  If formula provided and data is NULL, the observed dataframe

- unobserved:

  If formula provided and data is NULL, the unobserved dataframe

- require_predictions:

  Logical, whether to require prediction columns

## Value

An msd_data object ready for estimation

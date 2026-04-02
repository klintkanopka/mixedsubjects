# Apply formula to create/update msd_data

Takes a parsed formula and either creates a new msd_data object or
updates column references in an existing one.

## Usage

``` r
apply_formula_to_data(
  parsed_formula,
  data = NULL,
  observed = NULL,
  unobserved = NULL
)
```

## Arguments

- parsed_formula:

  Result from parse_msd_formula()

- data:

  Either raw dataframes or an msd_data object

- observed:

  Raw observed dataframe (if data is NULL)

- unobserved:

  Raw unobserved dataframe (if data is NULL)

## Value

An msd_data object with appropriate column mapping

# Parse MSD formula

Parses formulas of the form: outcome ~ treatment \| predictions

## Usage

``` r
parse_msd_formula(formula)
```

## Arguments

- formula:

  A formula object

## Value

A list with outcome, treatment, and prediction variable names

## Details

For GREG-type estimators (single prediction per arm): Y ~ D \| S

For DiP-type estimators (both predictions): Y ~ D \| S1 + S0

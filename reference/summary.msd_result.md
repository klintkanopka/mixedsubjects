# Summary method for msd_result

Produces a detailed summary of the MSD estimation results, including the
coefficient table with z-statistics and p-values, sample size
information, and interpretation guidance.

## Usage

``` r
# S3 method for class 'msd_result'
summary(object, ...)
```

## Arguments

- object:

  An msd_result object

- ...:

  Additional arguments (ignored)

## Value

A summary.msd_result object (invisibly when printed)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- msd_dt_dip(msd)
summary(result)
} # }
```

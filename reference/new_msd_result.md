# Utility Functions for Mixed-Subjects Design

Internal helper functions for the mixedsubjects package.

## Usage

``` r
new_msd_result(
  estimate,
  variance,
  se,
  ci_lower,
  ci_upper,
  method,
  lambda = NULL,
  n1 = NULL,
  n0 = NULL,
  m = NULL,
  conf_level = 0.95,
  additional = list()
)
```

## Arguments

- estimate:

  Point estimate of the ATE

- variance:

  Estimated variance

- se:

  Standard error

- ci_lower:

  Lower bound of confidence interval

- ci_upper:

  Upper bound of confidence interval

- method:

  Name of the estimation method

- lambda:

  Tuning parameter(s), if applicable

- n1:

  Number of treated observed units

- n0:

  Number of control observed units

- m:

  Number of unobserved (unlabeled) units

- conf_level:

  Confidence level used

- additional:

  Additional method-specific information

## Value

An S3 object of class "msd_result"

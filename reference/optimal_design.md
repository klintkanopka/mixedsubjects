# Optimal Design for Mixed-Subjects Experiments

Find optimal budget allocation between human observations and
predictions. Optimal Design Selection

Determines the optimal allocation of budget between human observations
and LLM predictions, and recommends the best estimator for the given
pilot data.

## Usage

``` r
optimal_design(
  pilot_data,
  budget,
  cost_human,
  cost_prediction,
  treatment_prob = 0.5,
  estimators = "all",
  min_observed = 20,
  n_grid = 100
)
```

## Arguments

- pilot_data:

  An msd_data object from a pilot study

- budget:

  Total budget available (in dollars)

- cost_human:

  Cost per human observation (in dollars)

- cost_prediction:

  Cost per LLM prediction (in dollars)

- treatment_prob:

  Probability of treatment assignment (default 0.5)

- estimators:

  Which estimators to consider. Either "all" or a character vector with
  subset of: "dim", "greg", "ppi", "dt", "dip", "dip_pp", "dt_dip"

- min_observed:

  Minimum number of observed units required (default 20)

- n_grid:

  Number of grid points for optimization (default 100)

## Value

An S3 object of class "msd_design" containing:

- optimal_n_obs:

  Recommended number of observed (human) units

- optimal_n_unobs:

  Recommended number of unobserved (prediction) units

- optimal_estimator:

  Recommended estimator

- optimal_variance:

  Expected variance at the optimum

- optimal_se:

  Expected standard error at the optimum

- budget_used:

  Total budget used

- all_results:

  Full grid search results for all estimators

## Details

The function uses grid search to find the optimal (n_O, n_U) allocation
that minimizes expected variance given the budget constraint:

\$\$n_O \times cost\_{human} + n_U \times (k \times cost\_{prediction})
\leq budget\$\$

where k is the number of predictions per unit:

- k = 0 for DiM (no predictions)

- k = 1 for GREG, PPI++, D-T (one prediction per arm)

- k = 2 for DiP, DiP++, D-T DiP (both S^(0) and S^(1))

The expected variance is computed using population moments estimated
from the pilot data.

## Examples

``` r
# Pilot study data
pilot_obs <- data.frame(
  Y = rnorm(50),
  S0 = rnorm(50),
  S1 = rnorm(50),
  D = rep(c(1, 0), each = 25)
)
pilot_unobs <- data.frame(
  S0 = rnorm(100),
  S1 = rnorm(100),
  D = rep(c(1, 0), each = 50)
)
pilot <- msd_data(observed = pilot_obs, unobserved = pilot_unobs)

# Find optimal design with $10,000 budget
design <- optimal_design(
  pilot_data = pilot,
  budget = 10000,
  cost_human = 10,      # $10 per human observation
  cost_prediction = 0.01 # $0.01 per prediction
)
print(design)
#> 
#> Optimal Mixed-Subjects Design
#> ==============================
#> 
#> Recommended Design:
#>   Estimator:         DT 
#>   Observed units:    990 
#>     - Treated:       495 
#>     - Control:       495 
#>   Unobserved units:  10000 
#> 
#> Expected Performance:
#>   Variance:  0.0051 
#>   SE:        0.0712 
#> 
#> Budget:
#>   Total:    $10,000
#>   Used:     $10,000
#>   Human:    $10/observation
#>   Predict:  $0.01/prediction
#> 
#> Comparison Across Estimators:
#>  Estimator n_obs n_unobs     SE
#>        DIM  1000       0 0.0724
#>       GREG   980   20000 0.0977
#>        PPI   990   10000 0.0727
#>         DT   990   10000 0.0712
#>        DIP   980   10000 0.0979
#>     DIP_PP   990    5000 0.0727
#>     DT_DIP   990    5000 0.0712
#> 
```

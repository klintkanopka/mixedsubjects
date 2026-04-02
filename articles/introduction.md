# Introduction to Mixed-Subjects Designs with the mixedsubjects Package

## Overview

The `mixedsubjects` package provides tools for analyzing randomized
experiments that combine traditional human subjects data with
predictions from large language models (LLMs) or other machine learning
algorithms. This approach, which we call a **Mixed-Subjects Design
(MSD)**, can substantially increase the precision of treatment effect
estimates while maintaining valid statistical inference.

### Why Mixed-Subjects Designs?

Running experiments with human subjects is expensive. A typical survey
response might cost several dollars, and detecting small treatment
effects often requires thousands of participants. Meanwhile, LLMs can
generate predictions about how humans might respond to experimental
treatments at a fraction of the cost—often less than a cent per
prediction.

The key insight behind mixed-subjects designs is that **LLM predictions,
even if imperfect, contain information about human behavior**. By
combining a smaller sample of actual human responses (the “gold
standard”) with a larger pool of LLM predictions, we can:

1.  **Reduce variance** in our treatment effect estimates
2.  **Maintain valid inference**—our confidence intervals will have
    correct coverage
3.  **Lower costs** of experimentation

Importantly, this approach **does not assume LLM predictions are
accurate**. Even if predictions are biased, the methods in this package
automatically adjust for this bias using the human subjects data.

### What This Package Does

The `mixedsubjects` package implements seven estimators for the Average
Treatment Effect (ATE) in mixed-subjects designs:

| Estimator   | Description                          | Predictions Needed |
|-------------|--------------------------------------|--------------------|
| **DiM**     | Difference-in-Means                  | None               |
| **GREG**    | Calibration estimator                | 1 per unit         |
| **PPI++**   | Power-tuned, single $\lambda$        | 1 per unit         |
| **D-T**     | Doubly-tuned, arm-specific $\lambda$ | 1 per unit         |
| **DiP**     | Difference-in-Predictions            | 2 per unit         |
| **DiP++**   | Power-tuned DiP                      | 2 per unit         |
| **D-T DiP** | Doubly-tuned DiP                     | 2 per unit         |

The package also helps you **design optimal experiments** by determining
how to allocate your budget between human subjects and LLM predictions.

## Installation

``` r
# Install from github
devtools::install_github('klintkanopka/mixedsubjects')
```

## A Quick Example

Let’s start with a simple example to see the package in action.

### Simulating Experimental Data

Imagine you’re running a survey experiment to test whether a persuasive
message changes attitudes. You have:

- **200 human respondents** who completed the survey (100 treated, 100
  control)
- **1,000 LLM predictions** of how additional respondents would have
  answered

``` r
set.seed(123)

# True treatment effect
true_ate <- 0.3

# Human subjects data (observed)
n_observed <- 200
observed_df <- data.frame(
  # Treatment assignment (balanced)
  D = rep(c(1, 0), each = n_observed / 2)
)

# Generate outcomes: Y(0) ~ N(0, 1), Y(1) ~ N(0.3, 1)
observed_df$Y <- ifelse(
  observed_df$D == 1,
  rnorm(n_observed / 2, mean = true_ate, sd = 1),
  rnorm(n_observed / 2, mean = 0, sd = 1)
)

# LLM predictions (correlated with true outcomes, but with some error)
# S^(1) predicts Y(1), S^(0) predicts Y(0)
observed_df$S1 <- 0.6 * observed_df$Y + rnorm(n_observed, 0, 0.5)
observed_df$S0 <- 0.5 * observed_df$Y + rnorm(n_observed, 0, 0.6)

# Unobserved units (only have predictions, no actual Y)
# Predictions must come from the same pipeline as observed (Assumption C: Random Labeling)
n_unobserved <- 1000
D_unobs <- rep(c(1, 0), each = n_unobserved / 2)
latent_Y <- ifelse(D_unobs == 1,
                   rnorm(n_unobserved / 2, mean = true_ate, sd = 1),
                   rnorm(n_unobserved / 2, mean = 0, sd = 1))
unobserved_df <- data.frame(
  D = D_unobs,
  S1 = 0.6 * latent_Y + rnorm(n_unobserved, 0, 0.5),
  S0 = 0.5 * latent_Y + rnorm(n_unobserved, 0, 0.6)
)
```

### Creating an MSD Data Object

First, we combine our data into an `msd_data` object:

``` r
msd <- msd_data(observed = observed_df, unobserved = unobserved_df)
print(msd)
#> 
#> Mixed-Subjects Design Data
#> ==========================
#> 
#> Sample Sizes:
#>   Observed (labeled):     200 
#>     - Treated (D=1):      100 
#>     - Control (D=0):      100 
#>   Unobserved (unlabeled): 1000 
#> 
#> Predictions Available:
#>   S0 (control arm):   Yes 
#>   S1 (treatment arm): Yes 
#> 
#> Column Mapping (original names):
#>   Observed:  Y=Y, D=D, S0=S0, S1=S1
#>   Unobserved: D=D, S0=S0, S1=S1
#> 
#> Available Estimators:
#>   - DiM:     Yes (no predictions needed)
#>   - GREG:    Yes
#>   - PPI++:   Yes
#>   - D-T:     Yes
#>   - DiP:     Yes (both S0 and S1 available)
#>   - DiP++:   Yes
#>   - D-T DiP: Yes
```

### Estimating the Treatment Effect

Now let’s estimate the ATE using different estimators:

``` r
# Classical difference-in-means (ignores predictions)
result_dim <- msd_dim(msd)

# D-T DiP (uses both predictions with arm-specific tuning)
result_dt_dip <- msd_dt_dip(msd)

# Compare
cat("Difference-in-Means:\n")
#> Difference-in-Means:
cat("  Estimate:", round(result_dim$estimate, 3), "\n")
#>   Estimate: 0.498
cat("  SE:", round(result_dim$se, 3), "\n")
#>   SE: 0.133
cat("  95% CI: [", round(result_dim$ci_lower, 3), ", ",
    round(result_dim$ci_upper, 3), "]\n\n")
#>   95% CI: [ 0.237 ,  0.759 ]

cat("D-T DiP:\n")
#> D-T DiP:
cat("  Estimate:", round(result_dt_dip$estimate, 3), "\n")
#>   Estimate: 0.168
cat("  SE:", round(result_dt_dip$se, 3), "\n")
#>   SE: 0.097
cat("  95% CI: [", round(result_dt_dip$ci_lower, 3), ", ",
    round(result_dt_dip$ci_upper, 3), "]\n")
#>   95% CI: [ -0.022 ,  0.357 ]
```

Notice that the D-T DiP estimator has a **smaller standard error** than
difference-in-means, even though both use the same human subjects data.
This variance reduction comes from leveraging the LLM predictions.

## Understanding the Data Structure

### What Data Do You Need?

For a mixed-subjects design, you need:

1.  **Observed (labeled) data**: Human subjects who were randomized and
    provided actual responses
    - `Y`: The outcome variable
    - `D`: Treatment indicator (1 = treated, 0 = control)
    - `S1` and/or `S0`: LLM predictions (optional for DiM)
2.  **Unobserved (unlabeled) data**: Units for which you only have LLM
    predictions
    - `D`: Treatment indicator
    - `S1` and/or `S0`: LLM predictions

### Two Types of Predictions

The package distinguishes between two prediction setups:

**Setup 1: Arm-specific predictions (1 prediction per unit)**

Each unit only has the prediction for their assigned treatment
condition: - Treated units have $S^{(1)}$ (prediction of their outcome
under treatment) - Control units have $S^{(0)}$ (prediction of their
outcome under control)

This is typical when you ask the LLM to predict each unit’s response to
their actual assigned condition.

**Setup 2: Both predictions (2 predictions per unit)**

Each unit has predictions for *both* conditions: - $S^{(1)}$: Prediction
of outcome if treated - $S^{(0)}$: Prediction of outcome if in control

This requires asking the LLM to predict each unit’s response under both
conditions, which doubles the prediction cost but enables the DiP family
of estimators.

### Creating Your Data Object

The
[`msd_data()`](http://klintkanopka.com/mixedsubjects/reference/msd_data.md)
function accepts data in two formats:

**Option 1: Separate dataframes**

``` r
msd <- msd_data(
  observed = observed_df,    # Must have Y, D, and optionally S0, S1
  unobserved = unobserved_df # Must have D and S0, S1 (no Y)
)
```

**Option 2: Combined dataframe**

If your data is in a single dataframe where unobserved units have
`Y = NA`:

``` r
# Combine into single dataframe
combined_df <- rbind(
  observed_df,
  data.frame(Y = NA, D = unobserved_df$D,
             S0 = unobserved_df$S0, S1 = unobserved_df$S1)
)

# Create msd_data object
msd_combined <- msd_data(data = combined_df)
print(msd_combined)
#> 
#> Mixed-Subjects Design Data
#> ==========================
#> 
#> Sample Sizes:
#>   Observed (labeled):     200 
#>     - Treated (D=1):      100 
#>     - Control (D=0):      100 
#>   Unobserved (unlabeled): 1000 
#> 
#> Predictions Available:
#>   S0 (control arm):   Yes 
#>   S1 (treatment arm): Yes 
#> 
#> Column Mapping (original names):
#>   Observed:  Y=Y, D=D, S0=S0, S1=S1
#>   Unobserved: D=D, S0=S0, S1=S1
#> 
#> Available Estimators:
#>   - DiM:     Yes (no predictions needed)
#>   - GREG:    Yes
#>   - PPI++:   Yes
#>   - D-T:     Yes
#>   - DiP:     Yes (both S0 and S1 available)
#>   - DiP++:   Yes
#>   - D-T DiP: Yes
```

### Flexible Column Names

By default, the package expects columns named `Y` (outcome), `D`
(treatment), `S0` (control prediction), and `S1` (treatment prediction).
If your columns use different names, specify them explicitly:

``` r
# Custom column names
my_data <- data.frame(
  response_var = rnorm(100),
  is_treated = rep(c(1, 0), each = 50),
  claude_pred_ctrl = rnorm(100),
  claude_pred_trt = rnorm(100)
)

msd <- msd_data(
  observed = my_data,
  outcome = "response_var",
  treatment = "is_treated",
  pred_control = "claude_pred_ctrl",
  pred_treated = "claude_pred_trt"
)
```

You can even use different column names for observed and unobserved
data:

``` r
msd <- msd_data(
  observed = obs_df,
  unobserved = unobs_df,
  obs_outcome = "Y",
  obs_treatment = "treated",
  obs_pred_control = "gpt_pred_0",
  obs_pred_treated = "gpt_pred_1",
  unobs_treatment = "D",
  unobs_pred_control = "S0",
  unobs_pred_treated = "S1"
)
```

### Formula Interface

For even more flexibility, you can use a formula interface directly in
the estimator functions. This is similar to the syntax used by `ivreg`
and other regression packages:

``` r
# Formula: outcome ~ treatment | predictions
# Using custom column names directly in the estimator

# For GREG-type estimators (single prediction per arm)
result <- msd_greg(response ~ treatment | pred_1 + pred_0,
                   observed = obs_df, unobserved = unobs_df)

# For DiP-type estimators (both predictions)
result <- msd_dt_dip(Y ~ D | S1 + S0,
                     observed = obs_df, unobserved = unobs_df)
```

The formula has three parts: 1. **Left of `~`**: The outcome variable
name 2. **Between `~` and `|`**: The treatment indicator variable name
3. **After `|`**: The prediction variable name(s)

For DiM (which doesn’t use predictions), you can omit the prediction
part:

``` r
result <- msd_dim(response ~ treated, observed = obs_df)
```

## The Seven Estimators

### 1. Difference-in-Means (DiM)

The classical estimator that ignores LLM predictions entirely:

$${\widehat{\tau}}^{DiM} = {\bar{Y}}_{1} - {\bar{Y}}_{0}$$

Use this as your baseline. It’s unbiased and requires no predictions.

``` r
result <- msd_dim(msd)
print(result)
#> 
#> Mixed-Subjects Design Estimation
#> =================================
#> Estimator: Difference-in-Means (DiM) 
#> 
#> Point Estimate:  0.4980 
#> Standard Error:  0.1330 
#> 95% CI:         [0.2373, 0.7586]
#> 
#> Sample Sizes:
#>   Observed:   n_1=100, n_0=100
```

### 2. GREG (Calibration Estimator)

GREG uses predictions to “calibrate” the estimate:

$${\widehat{\tau}}^{GREG} = \left\lbrack {\bar{S}}_{U}^{(1)} + \left( {\bar{Y}}_{1} - {\bar{S}}_{O_{1}}^{(1)} \right) \right\rbrack - \left\lbrack {\bar{S}}_{U}^{(0)} + \left( {\bar{Y}}_{0} - {\bar{S}}_{O_{0}}^{(0)} \right) \right\rbrack$$

This adjusts the prediction mean by the observed “rectifier” (the
average prediction error in the labeled data).

``` r
result <- msd_greg(msd)
print(result)
#> 
#> Mixed-Subjects Design Estimation
#> =================================
#> Estimator: GREG (lambda = 1) 
#> 
#> Point Estimate:  0.3307 
#> Standard Error:  0.1073 
#> 95% CI:         [0.1203, 0.5410]
#> 
#> Tuning Parameters:
#>   lambda:        1.0000 
#> 
#> Sample Sizes:
#>   Observed:   n_1=100, n_0=100
#>   Unobserved: |U|=1000
```

### 3. PPI++ (Power-Tuned)

PPI++ improves on GREG by estimating an optimal tuning parameter
$\lambda$:

$${\widehat{\mu}}_{d}^{PPI}(\lambda) = {\bar{Y}}_{d} + \lambda\left( {\bar{S}}_{U}^{(d)} - {\bar{S}}_{O}^{(d)} \right)$$

When $\lambda = 1$, this equals GREG. PPI++ estimates $\lambda$ to
minimize variance using cross-fitting.

``` r
result <- msd_ppi(msd, n_folds = 2)
print(result)
#> 
#> Mixed-Subjects Design Estimation
#> =================================
#> Estimator: PPI++ (cross-fit, K=2) 
#> 
#> Point Estimate:  0.3617 
#> Standard Error:  0.1033 
#> 95% CI:         [0.1593, 0.5642]
#> 
#> Tuning Parameters:
#>   lambda:        0.8222 
#> 
#> Sample Sizes:
#>   Observed:   n_1=100, n_0=100
#>   Unobserved: |U|=1000
```

### 4. D-T (Doubly-Tuned)

D-T uses separate tuning parameters for each arm ($\lambda_{1}$ for
treatment, $\lambda_{0}$ for control):

$${\widehat{\tau}}^{D - T} = {\widehat{\mu}}_{1}^{PPI}\left( \lambda_{1} \right) - {\widehat{\mu}}_{0}^{PPI}\left( \lambda_{0} \right)$$

This is useful when prediction quality differs between treatment and
control conditions.

``` r
result <- msd_dt(msd, n_folds = 2)
print(result)
#> 
#> Mixed-Subjects Design Estimation
#> =================================
#> Estimator: D-T (Doubly-Tuned, cross-fit, K=2) 
#> 
#> Point Estimate:  0.3457 
#> Standard Error:  0.1037 
#> 95% CI:         [0.1425, 0.5489]
#> 
#> Tuning Parameters:
#>   lambda_1 (treatment):  0.9039 
#>   lambda_0 (control):    0.7806 
#> 
#> Sample Sizes:
#>   Observed:   n_1=100, n_0=100
#>   Unobserved: |U|=1000
```

### 5. DiP (Difference-in-Predictions)

DiP uses the *contrast* $S^{(1)} - S^{(0)}$ for each unlabeled unit:

$${\widehat{\tau}}^{DiP} = {\overline{S^{(1)} - S^{(0)}}}_{U} + \left( {\bar{Y}}_{1} - {\bar{S}}_{O_{1}}^{(1)} \right) - \left( {\bar{Y}}_{0} - {\bar{S}}_{O_{0}}^{(0)} \right)$$

**Key insight**: When $S^{(1)}$ and $S^{(0)}$ are positively correlated
(which they often are—units that score high under treatment also tend to
score high under control), the variance of their difference is smaller
than the sum of their individual variances.

``` r
result <- msd_dip(msd)
print(result)
#> 
#> Mixed-Subjects Design Estimation
#> =================================
#> Estimator: DiP (Difference-in-Predictions, lambda = 1) 
#> 
#> Point Estimate:  0.1454 
#> Standard Error:  0.0976 
#> 95% CI:         [-0.0459, 0.3367]
#> 
#> Tuning Parameters:
#>   lambda:        1.0000 
#> 
#> Sample Sizes:
#>   Observed:   n_1=100, n_0=100
#>   Unobserved: |U|=1000
```

### 6. DiP++ (Power-Tuned DiP)

DiP++ applies power-tuning to the DiP estimator with a single $\lambda$:

$${\widehat{\tau}}^{DiP + +}(\lambda) = \lambda \cdot {\overline{S^{(1)} - S^{(0)}}}_{U} + \left( {\bar{Y}}_{1} - \lambda{\bar{S}}_{O_{1}}^{(1)} \right) - \left( {\bar{Y}}_{0} - \lambda{\bar{S}}_{O_{0}}^{(0)} \right)$$

``` r
result <- msd_dip_pp(msd, n_folds = 2)
print(result)
#> 
#> Mixed-Subjects Design Estimation
#> =================================
#> Estimator: DiP++ (cross-fit, K=2) 
#> 
#> Point Estimate:  0.1891 
#> Standard Error:  0.0965 
#> 95% CI:         [0.0000, 0.3782]
#> 
#> Tuning Parameters:
#>   lambda:        0.8804 
#> 
#> Sample Sizes:
#>   Observed:   n_1=100, n_0=100
#>   Unobserved: |U|=1000
```

### 7. D-T DiP (Doubly-Tuned DiP)

The most flexible estimator combines DiP with arm-specific tuning:

$${\widehat{\tau}}^{D - T\ DiP} = {\overline{\lambda_{1}S^{(1)} - \lambda_{0}S^{(0)}}}_{U} + \left( {\bar{Y}}_{1} - \lambda_{1}{\bar{S}}_{O_{1}}^{(1)} \right) - \left( {\bar{Y}}_{0} - \lambda_{0}{\bar{S}}_{O_{0}}^{(0)} \right)$$

``` r
result <- msd_dt_dip(msd, n_folds = 2)
print(result)
#> 
#> Mixed-Subjects Design Estimation
#> =================================
#> Estimator: D-T DiP (cross-fit, K=2) 
#> 
#> Point Estimate:  0.1571 
#> Standard Error:  0.0966 
#> 95% CI:         [-0.0323, 0.3465]
#> 
#> Tuning Parameters:
#>   lambda_1 (treatment):  0.9944 
#>   lambda_0 (control):    0.8544 
#> 
#> Sample Sizes:
#>   Observed:   n_1=100, n_0=100
#>   Unobserved: |U|=1000
```

### Comparing All Estimators

Use
[`estimate_all()`](http://klintkanopka.com/mixedsubjects/reference/estimate_all.md)
to run all applicable estimators at once:

``` r
all_results <- estimate_all(msd)
print(all_results)
#> 
#> Mixed-Subjects Design: All Estimators
#> ======================================
#> 
#>                                    Estimator Estimate     SE 95% CI Lower
#>                    Difference-in-Means (DiM)   0.4980 0.1330       0.2373
#>                            GREG (lambda = 1)   0.3307 0.1073       0.1203
#>                       PPI++ (cross-fit, K=2)   0.3664 0.1032       0.1643
#>           D-T (Doubly-Tuned, cross-fit, K=2)   0.3585 0.1034       0.1557
#>  DiP (Difference-in-Predictions, lambda = 1)   0.1454 0.0976      -0.0459
#>                       DiP++ (cross-fit, K=2)   0.1893 0.0965       0.0002
#>                     D-T DiP (cross-fit, K=2)   0.1772 0.0966      -0.0121
#>  95% CI Upper
#>        0.7586
#>        0.5410
#>        0.5686
#>        0.5612
#>        0.3367
#>        0.3785
#>        0.3665
```

## Choosing the Right Estimator

### Decision Tree

Here’s a practical guide for choosing an estimator:

    Do you have LLM predictions?
    ├── No → Use DiM
    └── Yes → Do you have BOTH S^(1) and S^(0) for each unit?
        ├── No (only arm-specific) →
        │   └── Is prediction quality similar across arms?
        │       ├── Yes → Use PPI++
        │       └── No → Use D-T
        └── Yes (both predictions) →
            └── Are S^(1) and S^(0) positively correlated?
                ├── Yes → Use D-T DiP (recommended)
                └── No/Unsure → Use D-T

### When to Use DiP-Type Estimators

The DiP family (DiP, DiP++, D-T DiP) is most beneficial when:

1.  **Predictions are positively correlated**:
    $Cov\left( S^{(1)},S^{(0)} \right) > 0$
2.  **You can afford 2 predictions per unit**: This doubles the
    prediction cost

You can check the correlation in your pilot data:

``` r
cor(unobserved_df$S1, unobserved_df$S0)
#> [1] 0.52053
```

A positive correlation (which is common) means DiP estimators will
likely outperform GREG-type estimators.

## Variance Estimation

### Delta-Method (Default)

All estimators use the delta-method for variance estimation by default.
This provides fast, analytical variance estimates based on the formulas
derived in the accompanying paper.

### Bootstrap (Optional)

For robustness, you can also use the fold-respecting bootstrap:

``` r
# Bootstrap variance for D-T DiP
boot_result <- bootstrap_variance(
  msd,
  estimator = "dt_dip",
  n_bootstrap = 500,
  seed = 42
)

cat("Delta-method SE:", round(result_dt_dip$se, 4), "\n")
#> Delta-method SE: 0.0966
cat("Bootstrap SE:", round(boot_result$se, 4), "\n")
#> Bootstrap SE: 0.0962
```

The bootstrap is computationally more expensive but provides a useful
check on the delta-method variance.

## Optimal Experimental Design

### The Design Problem

Before running your experiment, you need to decide:

1.  **How many human subjects** to recruit ($n_{O}$)
2.  **How many LLM predictions** to generate ($n_{U}$)
3.  **Which estimator** to use

These decisions depend on: - Your **total budget** - The **cost per
human subject** - The **cost per LLM prediction** - The **quality of
predictions** (estimated from pilot data)

### Using `optimal_design()`

The
[`optimal_design()`](http://klintkanopka.com/mixedsubjects/reference/optimal_design.md)
function finds the allocation that minimizes expected variance given
your budget constraints:

``` r
design <- optimal_design(
  pilot_data = msd,          # Pilot study data
  budget = 5000,             # Total budget in dollars
  cost_human = 5,            # $5 per human response
  cost_prediction = 0.01,    # $0.01 per LLM prediction
  treatment_prob = 0.5       # Balanced treatment assignment
)

print(design)
#> 
#> Optimal Mixed-Subjects Design
#> ==============================
#> 
#> Recommended Design:
#>   Estimator:         DT_DIP 
#>   Observed units:    970 
#>     - Treated:       485 
#>     - Control:       485 
#>   Unobserved units:  7500 
#> 
#> Expected Performance:
#>   Variance:  0.0019 
#>   SE:        0.0431 
#> 
#> Budget:
#>   Total:    $5,000
#>   Used:     $5,000
#>   Human:    $5/observation
#>   Predict:  $0.01/prediction
#> 
#> Comparison Across Estimators:
#>  Estimator n_obs n_unobs     SE
#>        DIM  1000       0 0.0595
#>       GREG   951   24500 0.0443
#>        PPI   960   20000 0.0439
#>         DT   960   20000 0.0438
#>        DIP   960   10000 0.0438
#>     DIP_PP   970    7500 0.0434
#>     DT_DIP   970    7500 0.0431
```

### Interpreting the Results

The output tells you:

1.  **Recommended estimator**: Which of the 7 estimators to use
2.  **Sample sizes**: How many observed and unobserved units
3.  **Expected SE**: The anticipated standard error

### Comparing Designs

The output also shows how different estimators perform at their optimal
allocations. This helps you understand the trade-offs:

- **DiM** uses all budget on human subjects (no predictions)
- **GREG/PPI++/D-T** use 1 prediction per unlabeled unit
- **DiP/DiP++/D-T DiP** use 2 predictions per unlabeled unit

## Practical Workflow

Here’s a recommended workflow for using mixed-subjects designs in your
research:

### Step 1: Pilot Study

Run a small pilot study with both human subjects and LLM predictions:

``` r
# Collect pilot data
pilot_observed <- collect_human_responses(n = 100)
pilot_observed$S1 <- get_llm_predictions(pilot_observed, condition = "treatment")
pilot_observed$S0 <- get_llm_predictions(pilot_observed, condition = "control")

pilot_unobserved <- generate_synthetic_units(n = 200)
pilot_unobserved$S1 <- get_llm_predictions(pilot_unobserved, condition = "treatment")
pilot_unobserved$S0 <- get_llm_predictions(pilot_unobserved, condition = "control")

pilot_msd <- msd_data(observed = pilot_observed, unobserved = pilot_unobserved)
```

### Step 2: Estimate Prediction Quality

Examine how well LLM predictions correlate with actual outcomes:

``` r
# Check correlations
pilot_results <- estimate_all(pilot_msd)
print(pilot_results)

# Compare to DiM baseline
dim_se <- pilot_results$SE[pilot_results$Estimator == "Difference-in-Means (DiM)"]
best_se <- min(pilot_results$SE)
variance_reduction <- 1 - (best_se / dim_se)^2
cat("Potential variance reduction:", round(variance_reduction * 100, 1), "%\n")
```

### Step 3: Plan Your Main Study

Use the pilot data to design your main study:

``` r
design <- optimal_design(
  pilot_data = pilot_msd,
  budget = 10000,
  cost_human = 5,
  cost_prediction = 0.01
)
print(design)
```

### Step 4: Run the Main Study

Collect data according to your design:

``` r
# Collect human responses
main_observed <- collect_human_responses(n = design$optimal_n_obs)

# Generate LLM predictions
main_observed$S1 <- get_llm_predictions(main_observed, "treatment")
main_observed$S0 <- get_llm_predictions(main_observed, "control")

main_unobserved <- generate_synthetic_units(n = design$optimal_n_unobs)
main_unobserved$S1 <- get_llm_predictions(main_unobserved, "treatment")
main_unobserved$S0 <- get_llm_predictions(main_unobserved, "control")

main_msd <- msd_data(observed = main_observed, unobserved = main_unobserved)
```

### Step 5: Analyze and Report

Estimate the treatment effect using the recommended estimator:

``` r
# Run the recommended estimator
result <- msd_dt_dip(main_msd)  # Or whichever was recommended
print(result)

# Optionally verify with bootstrap
boot_result <- bootstrap_variance(main_msd, "dt_dip", n_bootstrap = 1000)
```

## Common Questions

### Q: What if my LLM predictions are biased?

**A**: That’s fine! The MSD estimators are designed to be valid even
when predictions are biased. The human subjects data is used to
“rectify” any systematic bias in the predictions. You will still get
unbiased estimates with valid confidence intervals.

### Q: How much do predictions need to correlate with outcomes?

**A**: Any positive correlation helps, but more is better. As a rough
guide: - $\rho > 0.3$: Modest variance reduction - $\rho > 0.5$:
Substantial variance reduction - $\rho > 0.7$: Large variance reduction

You can estimate this correlation from your pilot data.

### Q: Should I use cross-fitting?

**A**: Yes, for the tuned estimators (PPI++, D-T, DiP++, D-T DiP).
Cross-fitting ensures that the tuning parameters are estimated
independently of the data used to compute the estimate, which prevents
bias. The default `n_folds = 2` works well in most cases.

### Q: What if predictions are worse in one treatment arm?

**A**: Use the D-T or D-T DiP estimators. These use arm-specific tuning
parameters ($\lambda_{1}$ and $\lambda_{0}$), which automatically adapt
to different prediction quality across arms.

### Q: How many folds should I use for cross-fitting?

**A**: We recommend `n_folds = 2` for most applications. More folds use
data more efficiently but increase computational cost. With 2 folds,
half the labeled data is used to estimate tuning parameters and half for
estimation in each fold.

### Q: Can I use predictions from multiple LLMs?

**A**: The current package supports one set of predictions. If you have
multiple LLMs, you could: 1. Use the single best-performing LLM 2.
Average predictions across LLMs before using them 3. Use the LLM whose
predictions have the highest correlation with outcomes

## Technical Details

### Assumptions

The MSD estimators require:

1.  **Random sampling**: Both observed and unobserved units are drawn
    from the same population
2.  **Random treatment assignment**: Treatment is randomly assigned in
    both pools
3.  **Independent predictions**: The LLM predictions don’t use the
    observed outcomes
4.  **Finite variances**: Standard regularity conditions

### Variance Formulas

For reference, the variance formulas for key estimators:

**DiM**:
$$Var\left( {\widehat{\tau}}^{DiM} \right) = \frac{\sigma_{Y{(1)}}^{2}}{n_{1}} + \frac{\sigma_{Y{(0)}}^{2}}{n_{0}}$$

**GREG**:
$$Var\left( {\widehat{\tau}}^{GREG} \right) = \sum\limits_{d \in \{ 0,1\}}\left\lbrack \frac{\sigma_{S^{(d)}}^{2}}{m_{d}} + \frac{Var\left( Y(d) - S^{(d)} \right)}{n_{d}} \right\rbrack$$

**DiP**:
$$Var\left( {\widehat{\tau}}^{DiP} \right) = \frac{Var\left( S^{(1)} - S^{(0)} \right)}{m} + \sum\limits_{d \in \{ 0,1\}}\frac{Var\left( Y(d) - S^{(d)} \right)}{n_{d}}$$

The tuned estimators (PPI++, D-T, DiP++, D-T DiP) have more complex
variance formulas that account for the estimation of tuning parameters.

## Citation

If you use this package in your research, please cite:

    van Loon, A. and Kanopka, K. (2025). Causal Inference in Experiments with
    Mixed-Subjects Designs. R package version 0.1.0.

## Session Info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] mixedsubjects_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.5         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.5.3    tools_4.5.3       ragg_1.5.2        evaluate_1.0.5   
#> [21] bslib_0.10.0      yaml_2.3.12       jsonlite_2.0.0    rlang_1.1.7      
#> [25] fs_2.0.1
```

# Data Preparation for Mixed-Subjects Design

Create and validate data objects for MSD estimation. Create an MSD data
object

Creates a validated data object for mixed-subjects design estimation.
Accepts either a single combined dataframe or two separate dataframes
for observed and unobserved units. Column names default to "Y", "D",
"S0", "S1" and can be overridden explicitly.

## Usage

``` r
msd_data(
  data = NULL,
  observed = NULL,
  unobserved = NULL,
  outcome = "Y",
  treatment = "D",
  pred_control = "S0",
  pred_treated = "S1",
  obs_outcome = NULL,
  obs_treatment = NULL,
  obs_pred_control = NULL,
  obs_pred_treated = NULL,
  unobs_treatment = NULL,
  unobs_pred_control = NULL,
  unobs_pred_treated = NULL
)
```

## Arguments

- data:

  A combined dataframe containing both observed and unobserved units.
  Observed units have non-missing Y values; unobserved units have Y =
  NA. If provided, `observed` and `unobserved` should be NULL.

- observed:

  A dataframe of observed (labeled) units with columns for outcome Y,
  predictions S0/S1, and treatment D. If provided, `data` should be
  NULL.

- unobserved:

  A dataframe of unobserved (unlabeled) units with columns for
  predictions S0/S1 and treatment D (no Y column needed).

- outcome:

  Name of the outcome column. Default: "Y".

- treatment:

  Name of the treatment column. Default: "D".

- pred_control:

  Name of the control prediction column. Default: "S0".

- pred_treated:

  Name of the treatment prediction column. Default: "S1".

- obs_outcome:

  Outcome column name for observed data only (overrides `outcome`).

- obs_treatment:

  Treatment column name for observed data only (overrides `treatment`).

- obs_pred_control:

  Control prediction column for observed data only.

- obs_pred_treated:

  Treatment prediction column for observed data only.

- unobs_treatment:

  Treatment column name for unobserved data only.

- unobs_pred_control:

  Control prediction column for unobserved data only.

- unobs_pred_treated:

  Treatment prediction column for unobserved data only.

## Value

An S3 object of class "msd_data" containing:

- observed:

  Dataframe of observed units with standardized columns Y, S0, S1, D

- unobserved:

  Dataframe of unobserved units with standardized columns S0, S1, D (or
  NULL)

- has_S0:

  Logical indicating if S0 predictions are available

- has_S1:

  Logical indicating if S1 predictions are available

- has_both_predictions:

  Logical indicating if both S0 and S1 are available

- n1:

  Number of treated observed units

- n0:

  Number of control observed units

- m:

  Number of unobserved units

- col_mapping:

  List of original column names used

## Details

The function supports flexible column name specification:

**Default column names:** By default, the function expects columns named
"Y" (outcome), "D" (treatment), "S0" (control prediction), and "S1"
(treatment prediction). Override these using the `outcome`, `treatment`,
`pred_control`, and `pred_treated` arguments.

**Global specification:** Use `outcome`, `treatment`, `pred_control`,
`pred_treated` to set column names that apply to both observed and
unobserved dataframes.

**Per-dataframe specification:** Use `obs_*` and `unobs_*` arguments
when column names differ between observed and unobserved dataframes.
These override global settings.

**Two input modes:**

*Mode 1: Single combined dataframe* Provide a single dataframe via the
`data` argument. The function will automatically split it into observed
and unobserved based on whether Y is NA.

*Mode 2: Separate dataframes* Provide two separate dataframes via
`observed` and `unobserved`.

## Examples

``` r
# Default column names (Y, D, S0, S1)
obs_df <- data.frame(
  Y = c(1.2, 0.8, 1.5, 0.9),
  S0 = c(1.0, 0.7, 1.3, 0.8),
  S1 = c(1.1, 0.9, 1.4, 1.0),
  D = c(1, 0, 1, 0)
)
msd <- msd_data(observed = obs_df)

# Custom column names (same in both dataframes)
obs_df2 <- data.frame(
  response = c(1.2, 0.8, 1.5, 0.9),
  pred_ctrl = c(1.0, 0.7, 1.3, 0.8),
  pred_trt = c(1.1, 0.9, 1.4, 1.0),
  treated = c(1, 0, 1, 0)
)
unobs_df2 <- data.frame(
  pred_ctrl = c(1.1, 0.9),
  pred_trt = c(1.2, 1.0),
  treated = c(1, 0)
)
msd2 <- msd_data(
  observed = obs_df2,
  unobserved = unobs_df2,
  outcome = "response",
  treatment = "treated",
  pred_control = "pred_ctrl",
  pred_treated = "pred_trt"
)

# Different column names in observed vs unobserved
obs_df3 <- data.frame(
  outcome = c(1.2, 0.8),
  claude_pred_0 = c(1.0, 0.7),
  claude_pred_1 = c(1.1, 0.9),
  treatment = c(1, 0)
)
unobs_df3 <- data.frame(
  s0_claude = c(1.1, 0.9),
  s1_claude = c(1.2, 1.0),
  D = c(1, 0)
)
msd3 <- msd_data(
  observed = obs_df3,
  unobserved = unobs_df3,
  obs_outcome = "outcome",
  obs_treatment = "treatment",
  obs_pred_control = "claude_pred_0",
  obs_pred_treated = "claude_pred_1",
  unobs_treatment = "D",
  unobs_pred_control = "s0_claude",
  unobs_pred_treated = "s1_claude"
)
```

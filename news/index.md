# Changelog

## mixedsubjects 1.1.0

This release brings the package in line with the finite-sample theory in
the Supplementary Information of the accompanying paper. Standard errors
change for all four cross-fit estimators, and point estimates change
slightly for three of them. Analyses run with 1.0.0 should be rerun.

### Variance of the cross-fit estimators

- New analytic variance for
  [`msd_ppi()`](http://klintkanopka.com/mixedsubjects/reference/msd_ppi.md),
  [`msd_dt()`](http://klintkanopka.com/mixedsubjects/reference/msd_dt.md),
  [`msd_dip_pp()`](http://klintkanopka.com/mixedsubjects/reference/msd_dip_pp.md)
  and
  [`msd_dt_dip()`](http://klintkanopka.com/mixedsubjects/reference/msd_dt_dip.md).
  The variance is now assembled fold by fold: each fold’s fixed-tuning
  variance is evaluated at the tuning used in that fold with that fold’s
  own moments, the shared-unlabeled covariance uses the realized fold
  tunings, and a new second-order cross-fold term accounts for the
  dependence between each fold’s tuning estimate and the labeled data it
  was estimated on. The cross-fold term is estimated in closed form from
  third central moments. Standard errors are slightly larger than in
  1.0.0, most visibly with small labeled samples and skewed outcomes or
  predictions.
- `additional$variance_components` reports the three pieces of the new
  variance (`fold_variances`, `cross_unlabeled`, `cross_fold`), which
  sum to the reported variance. `additional$variance_fixed_lambda`
  reports the 1.0.0 fixed-tuning plug-in for comparison.

### Tuning from labeled data only (changes point estimates)

- [`msd_dip_pp()`](http://klintkanopka.com/mixedsubjects/reference/msd_dip_pp.md)
  and
  [`msd_dt_dip()`](http://klintkanopka.com/mixedsubjects/reference/msd_dt_dip.md)
  now estimate their tuning parameters from labeled data only. In 1.0.0
  the second moments of the predictions that enter the tuning rule with
  weight 1/m were estimated from the unlabeled pool, which made the
  tuning depend on data the estimator is evaluated on. They are now
  estimated from the labeled tuning folds (every labeled unit carries
  both predictions), so the tuning is independent of the unlabeled pool
  and of the fold it is applied to, as the unbiasedness results require.
  Point estimates and tuning parameters from these two functions change
  slightly.

### Fold weights for D-T (changes point estimates when folds are unequal)

- [`msd_dt()`](http://klintkanopka.com/mixedsubjects/reference/msd_dt.md)
  now combines the fold-specific arm means with weights equal to each
  fold’s share of the arm’s labeled units, instead of equal weights. The
  two coincide whenever each arm’s labeled sample size is divisible by
  `n_folds`; otherwise point estimates change slightly. The weights are
  returned in `additional$fold_weights`, and the reported `lambda` is
  the correspondingly weighted average of the fold tunings.

### Unchanged

- [`msd_dim()`](http://klintkanopka.com/mixedsubjects/reference/msd_dim.md),
  [`msd_greg()`](http://klintkanopka.com/mixedsubjects/reference/msd_greg.md)
  and
  [`msd_dip()`](http://klintkanopka.com/mixedsubjects/reference/msd_dip.md)
  are unaffected. Point estimates from
  [`msd_ppi()`](http://klintkanopka.com/mixedsubjects/reference/msd_ppi.md)
  are unchanged.

## mixedsubjects 1.0.0

CRAN release: 2026-07-02

- Initial release.
- Seven estimators for the average treatment effect in mixed-subjects
  designs: Difference-in-Means, GREG, PPI++, Doubly-Tuned (D-T),
  Difference-in-Predictions (DiP), DiP++, and D-T DiP.
- Variance estimation via delta-method or fold-respecting bootstrap.
- Optimal design selection for allocating budget between human
  observations and LLM predictions.

# Compute PPI++ variance using delta method

The cross-fit estimator is theta = (1/K) \* sum_k theta_k, where each
theta_k = L_k + U. L_k is the labeled component (independent across
folds, using n_d/K obs each) and U is the shared unobserved component
(same for all folds). Therefore: Var(theta) = (1/K^2) \* sum_k
Var(L_k) + Var(U) = labeled_var / n_d + lambda^2 \* Var(S_d) / m_d The
labeled terms' K factors cancel; the unobserved term is NOT divided by
K.

## Usage

``` r
compute_ppi_variance(data, lambda, n_folds)
```

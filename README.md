# mixedsubjects

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/mixedsubjects)](https://CRAN.R-project.org/package=mixedsubjects)
[![R-CMD-check](https://github.com/klintkanopka/mixedsubjects/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/klintkanopka/mixedsubjects/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->


`mixedsubjects` is a package for conducting social science experiments using the Mixed-Subjects Design and estimating causal effects. It implements seven estimators for average treatment effect (ATE) estimation in mixed-subjects designs (MSDs), where human subjects data is augmented with predictions from large language models (LLMs). Includes Difference-in-Means, GREG, PPI++, Doubly-Tuned, Difference-in-Predictions (DiP), DiP++, and D-T DiP estimators. Provides point estimates, variance estimation via delta-method or bootstrap, and optimal design selection for budget allocation between human observations and LLM predictions.

## Installation

The most recent CRAN release can be installed using:

```r
install.packages('mixedsubjects')
```

## Development version

Adventurous users can install the current development version using:

```r
devtools::install_github('klintkanopka/mixedsubjects')
```

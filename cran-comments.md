## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

The single NOTE is from the CRAN incoming feasibility check:

* "New submission" — this is a first submission.
* "Possibly misspelled words in DESCRIPTION: DiP, MSDs, PPI" — these are
  correctly spelled domain-specific abbreviations (DiP = Difference-in-
  Predictions, MSDs = Mixed-Subjects Designs, PPI = Prediction-Powered
  Inference), not misspellings.

## Test environments

* local macOS, R 4.5.2
* GitHub Actions (via .github/workflows/R-CMD-check.yaml):
  * macOS-latest (release)
  * windows-latest (release)
  * ubuntu-latest (devel, release, oldrel-1)
* win-builder (devel and release)

## Notes for submission

* Examples run in well under the per-example time limit; the bootstrap-based
  examples use a reduced number of replications.
* Vignettes build with knitr/rmarkdown and require pandoc.

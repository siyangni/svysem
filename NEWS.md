# svysem 0.0.1

* CRAN check fixes. `coef`, `vcov` and `summary` are registered S3 methods, the unused `numDeriv` import is removed, and examples that need `mitools` are guarded.

# svysem 0.0.0.9000

* MLM estimation for `survey.design` and `survey.design2`. The design-based covariance of the sample moments is passed to lavaan as `NACOV`, with the raw data and sampling weights.
* Multiple imputation. A `svyimputationList` is fitted one completed design at a time and pooled with Rubin's rules (`pooled_estimates()`, `pooled_wald()`, `fitMeasures()`). The design degrees of freedom are the complete-data degrees of freedom.

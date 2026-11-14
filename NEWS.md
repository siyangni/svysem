# svysem 0.0.1

* CRAN check fixes. `coef`, `vcov` and `summary` are registered S3 methods, the unused `numDeriv` import is removed, and examples that need `mitools` are guarded.
* Examples fit simulated designs. They do not register methods in the `stats` namespace.
* `pooled_wald()` accepts an unlabeled model. lavaan omits the label column when no parameter is labeled.
* The description spell check accepts Oberski, Satorra, Bentler, and subpopulation.
* R 4.2.0 or newer is required.
* `coef()` and `vcov()` on a pooled object are also S4 methods, so they still run after lavaan turns those generics into S4 methods. Standardized pooling requests lavaan standard errors with `se = TRUE`.
* A `svyrep.design` is refused with `svysem_error_design`. Replicate-weight designs are not supported yet.
* `pooled_wald()` keeps a finite p-value when the small-sample denominator is undefined or negative. One constraint then uses the scalar Barnard-Rubin degrees of freedom.

# svysem 0.0.0.9000

* MLM estimation for `survey.design` and `survey.design2`. The design-based covariance of the sample moments is passed to lavaan as `NACOV`, with the raw data and sampling weights.
* Multiple imputation. A `svyimputationList` is fitted one completed design at a time and pooled with Rubin's rules (`pooled_estimates()`, `pooled_wald()`, `fitMeasures()`). The design degrees of freedom are the complete-data degrees of freedom.

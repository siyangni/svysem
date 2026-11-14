# svysem 0.0.0.9000

* MLM estimation for `survey.design` and `survey.design2`. The design-based covariance of the sample moments is passed to lavaan as `NACOV`, with the raw data and sampling weights.
* Multiple imputation. A `svyimputationList` is fitted one completed design at a time and pooled with Rubin's rules (`pooled_estimates()`, `pooled_wald()`, `fitMeasures()`). The design degrees of freedom are the complete-data degrees of freedom. `coef()` and `vcov()` on the pooled object are S4 methods. Standardized pooling requests lavaan standard errors with `se = TRUE`.
* A `svyrep.design` is refused with `svysem_error_design`. Replicate-weight designs are not supported yet.
* `pooled_wald()` keeps a finite p-value when the small-sample denominator is undefined or negative. One constraint then uses the scalar Barnard-Rubin degrees of freedom.

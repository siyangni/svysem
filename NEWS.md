# svysem 0.3.1

* `parameterEstimates()` and `standardizedSolution()` on a svysem fit include standard errors and confidence intervals. lavaan's own functions omit those columns unless `se = TRUE`, because the object's first class is `svysem` rather than `lavaan`.
* MLR on lavaan 0.7 uses `lavScores()` as returned. That version already multiplies casewise scores by the sampling weight. Multiplying again made the sandwich too large.
* The comparison vignette cites Krewski and Rao (1981) for the stratified jackknife and balanced repeated replication.

# svysem 0.3.0

* Multigroup MLM (`group=`) uses replicate-weight standard errors. A per-group `NACOV` would drop the covariance between groups when primary sampling units overlap the groups. The chi-square for that fit is the maximum-likelihood statistic.
* `predict()` returns posterior-mean factor scores for continuous outcomes, and plausible values when `draws` is positive. `lavaan::lavPredict()` on a weighted fit in lavaan 0.6-21 does not return these scores.
* WLSMV for ordered indicators. Point estimates use lavaan's weighted diagonally weighted least squares matrix. Standard errors use replicate weights (`survey::as.svrepdesign()`, type `"auto"` unless `replicate.type` is set). Continuous WLSMV passes `ordered = FALSE`, which lavaan 0.7 requires when no indicator is ordered.
* A `svyrep.design` is accepted for MLM, WLSMV, and MLR. `variance_method` is then `"replication"`.
* MLR for continuous outcomes, including `missing = "ml"`. The variance is a design-based sandwich from casewise scores, and the Yuan-Bentler statistic uses that sandwich. Listwise deletion keeps incomplete rows in the design with weight 0. Defined-parameter standard errors use the delta method.
* Multiple imputation accepts MLM, WLSMV, and MLR. `parallel` loops imputations; replicate weights inside one imputation stay sequential.
* Three vignettes: a design-based fit, multiply imputed designs, and the comparison with lavaan, lavaan.survey, and Mplus.

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

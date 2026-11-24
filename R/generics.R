#' Design-based structural equation models.
#'
#' `svysem()`, `svycfa()` and `svygrowth()` fit a lavaan model to a
#' [survey::svydesign()] object. Sampling weights, strata and primary sampling
#' units enter the standard errors through the design-based covariance of the
#' sample moments (Gamma), which is passed to lavaan as `NACOV` together with
#' the raw data and the sampling weights. The fit is not built from sample
#' moments. `svysem()` uses the defaults of [lavaan::sem()], `svycfa()` those
#' of [lavaan::cfa()], and `svygrowth()` those of [lavaan::growth()].
#' Arguments in `...` are passed on to that fitter.
#'
#' The returned object extends lavaan, so [coef()], [vcov()],
#' [lavaan::lavTestWald()] and [lavaan::fitMeasures()] report design-based
#' results. [parameterEstimates()][parameterEstimates] and
#' [standardizedSolution()][standardizedSolution] include the
#' design-based standard errors and confidence intervals. lavaan's own
#' functions omit those columns unless `se = TRUE`, because the object's
#' class is `svysem` rather than `lavaan`. [summary()] and the printed
#' display add the design degrees of freedom and the estimator.
#'
#' `"MLM"` is maximum likelihood with a Satorra-Bentler adjustment for
#' complete continuous data. `missing = "ml"` is an error and names MLR as the
#' estimator that handles it. Rows that are incomplete on the model variables
#' are not dropped: their sampling weight is set to 0 inside the design used
#' for Gamma, which is how [subset()] represents a domain, and the primary
#' sampling unit structure is kept. Zero weights are legal. A negative weight,
#' or a design with no positive weight, is an error.
#'
#' lavaan's default `sampling.weights.normalization = "total"` rescales
#' weights so that they sum to the number of rows. Means and covariances of
#' products are invariant to that constant, so Gamma is computed from the
#' design's sampling weights on the same rows. The products use the maximum
#' likelihood divisor (the weighted residual products, with no `n / (n - 1)`
#' factor). [survey::svyvar()] is not used. Means and second-moment products
#' are stacked in one linearization, so the covariance between them is kept;
#' that term is what makes a regression intercept match [survey::svyglm()].
#'
#' Gamma follows lavaan's model order (`ov.order = "force.model"` on the fit
#' that receives `NACOV`). Without that argument lavaan reorders variables to
#' the data order and pairs Gamma with the wrong columns. `fixed.x = TRUE`,
#' the default for [lavaan::sem()], is supported when those columns are the
#' observed second moments, plus the means when a mean structure is estimated.
#' `conditional.x = TRUE` is not supported for MLM or MLR. MLR fits one group
#' of continuous variables. Single-group MLM does too, by linearization.
#' `group=` on MLM refits the weighted model on replicate weights, because
#' lavaan keeps one `NACOV` per group and that drops the covariance between
#' groups. The multigroup chi-square is the maximum-likelihood statistic, not
#' a Satorra-Bentler statistic. `fitMeasures()` does not report a scaled
#' chi-square for that fit: the scaled statistic lavaan adds from the sampling
#' weights ignores primary sampling units, so it is dropped. Clustering is
#' taken from the design; do not
#' pass `cluster=`, which would double-count primary sampling units.
#' `survey.lonely.psu` is read by [survey::svymean()] and is never set here.
#'
#' `"MLR"` is robust maximum likelihood with a Yuan-Bentler adjustment. The
#' variance is a sandwich: lavaan's observed information, and a meat from the
#' design variance of the casewise scores ([survey::svytotal()]). Equality
#' constraints, including a shared label, `==`, and a nonlinear `==`, are
#' included. The sandwich is solved on lavaan's constraint basis and expanded
#' back to the free parameters. `std.lv` and `effect.coding` are accepted.
#' `missing = "ml"` scores the observed indicators; it does not fall back to
#' complete cases. On a `survey.design` the variance method is linearization.
#' On a `svyrep.design` it is replication, because [survey::svytotal()] uses
#' the replicate weights.
#'
#' `"WLSMV"` is lavaan's weighted least squares with the diagonally weighted
#' least squares matrix taken from the weighted sample. Point estimates match
#' that weighted fit. Standard errors and the scaled test use the
#' replicate-weight covariance of lavaan's sample statistics (thresholds,
#' continuous means and variances, and correlations), with the scale from
#' [survey::as.svrepdesign()] and [survey::svrVar()]. `replicate.type = "auto"`
#' is a stratified jackknife when the design has strata and JK1 otherwise.
#' `replicates` is only the bootstrap count. Ordered indicators, `group=`, and
#' `missing = "pairwise"` or `"listwise"` are accepted. lavaan stores one
#' `NACOV` per group, so the covariance of sample statistics between groups is
#' not passed on. The weight matrix Mplus uses for `TYPE = COMPLEX` has not
#' been compared with lavaan's. For WLSMV, MLR, and multigroup MLM, `se` and
#' `test` are set by the design and are not taken from `...`.
#'
#' A `svyrep.design` is accepted for every estimator. `variance_method` is
#' then `"replication"`. `replicates` and `replicate.type` are an error on a
#' design that already has replicate weights.
#'
#' [predict()] returns factor scores for continuous outcomes. See
#' [predict.svysem()].
#'
#' @param model Model syntax, as for [lavaan::sem()], [lavaan::cfa()] or
#'   [lavaan::growth()].
#' @param design A `survey.design`, `survey.design2`, or `svyrep.design`
#'   object, including the result of [subset()] or [survey::as.svrepdesign()].
#' @param estimator Variance estimator. `"MLM"` (Satorra-Bentler for one
#'   group, replicate weights when `group=` is set; complete continuous data),
#'   `"WLSMV"` (ordinal weighted least squares with replicate-weight standard
#'   errors), or `"MLR"` (Yuan-Bentler, including `missing = "ml"`).
#' @param replicates Number of bootstrap replicate weights. Used when
#'   `replicate.type = "bootstrap"` for WLSMV or for MLM with `group=`.
#'   Jackknife and BRR counts are fixed by the primary sampling units.
#' @param replicate.type How [survey::as.svrepdesign()] builds weights for
#'   WLSMV and for MLM with `group=`. One of `"auto"`, `"JK1"`, `"JKn"`,
#'   `"BRR"`, or `"bootstrap"`. `NULL` is `"auto"`.
#' @param parallel `NULL`, or a function of `(X, FUN)` with the same contract
#'   as [lapply()], used to loop WLSMV replicates and multigroup MLM
#'   replicates. On a `svyimputationList`, `parallel` loops imputations
#'   instead, and the replicates inside each imputation stay sequential.
#' @param ... Passed to [lavaan::sem()], [lavaan::cfa()] or [lavaan::growth()].
#'   `data`, `sampling.weights`, sample moments, `NACOV`, and `WLS.V` are
#'   taken from the design. `ov.order` is set internally so Gamma matches the
#'   model.
#' @return An S4 object of class `svysem` that contains the lavaan fit. Slots
#'   are `df` (design degrees of freedom), `estimator` (`"MLM"`, `"WLSMV"`,
#'   or `"MLR"`), `variance_method` (`"linearization"` or `"replication"`)
#'   and `svysem_call`.
#' @references
#' Oberski D (2014). "lavaan.survey: An R Package for Complex Survey Analysis
#' of Structural Equation Models." Journal of Statistical Software, 57(1),
#' 1--27. \doi{10.18637/jss.v057.i01}
#'
#' Yuan KH, Bentler PM (2000). "Three Likelihood-Based Methods for Mean and
#' Covariance Structure Analysis with Nonnormal Missing Data." Sociological
#' Methodology, 30(1), 165--200. \doi{10.1111/0081-1750.00078}
#' @seealso [svysem-class], [predict.svysem()], [simulate_design()],
#'   [lavaan::sem()], [survey::svydesign()], [survey::svyglm()]
#' @examples
#' design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 1)
#' fit <- svysem(
#'   "y ~ x1 + x2", design,
#'   fixed.x = FALSE, meanstructure = TRUE
#' )
#' fit
#' coef(fit)
#' vcov(fit)
#' summary(fit)
#' cfa_fit <- svycfa("f =~ y1 + y2 + y3 + y4", design)
#' growth_fit <- svygrowth(
#'   "i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
#'    s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4",
#'   design
#' )
#' coef(cfa_fit)
#' coef(growth_fit)
#' predict(cfa_fit)[1, , drop = FALSE]
#' mlr_fit <- svycfa("f =~ y1 + y2 + y3 + y4", design, estimator = "MLR")
#' coef(mlr_fit)
#' @export
svysem <- function(model, design, estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_prepare_dispatch(estimator, list(...))
  UseMethod("svysem", design)
}

#' @rdname svysem
#' @export
svycfa <- function(model, design, estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_prepare_dispatch(estimator, list(...))
  UseMethod("svycfa", design)
}

#' @rdname svysem
#' @export
svygrowth <- function(model, design, estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_prepare_dispatch(estimator, list(...))
  UseMethod("svygrowth", design)
}

#' @rdname svysem
#' @export
svysem.survey.design2 <- function(model, design,
                                  estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::sem,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...),
    replicates = replicates,
    replicate.type = replicate.type,
    parallel = parallel
  )
}

#' @rdname svysem
#' @export
svysem.survey.design <- function(model, design,
                                 estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::sem,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...),
    replicates = replicates,
    replicate.type = replicate.type,
    parallel = parallel
  )
}

#' @rdname svysem
#' @export
svysem.svyrep.design <- function(model, design,
                                 estimator = c("MLM", "WLSMV", "MLR"),
                                 replicates = NULL,
                                 replicate.type = NULL,
                                 parallel = NULL, ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::sem,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...),
    replicates = replicates,
    replicate.type = replicate.type,
    parallel = parallel
  )
}

#' @rdname svysem
#' @export
svysem.default <- function(model, design,
                           estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_stop_design(design)
}

#' @rdname svysem
#' @export
svycfa.survey.design2 <- function(model, design,
                                  estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::cfa,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...),
    replicates = replicates,
    replicate.type = replicate.type,
    parallel = parallel
  )
}

#' @rdname svysem
#' @export
svycfa.survey.design <- function(model, design,
                                 estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::cfa,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...),
    replicates = replicates,
    replicate.type = replicate.type,
    parallel = parallel
  )
}

#' @rdname svysem
#' @export
svycfa.svyrep.design <- function(model, design,
                                 estimator = c("MLM", "WLSMV", "MLR"),
                                 replicates = NULL,
                                 replicate.type = NULL,
                                 parallel = NULL, ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::cfa,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...),
    replicates = replicates,
    replicate.type = replicate.type,
    parallel = parallel
  )
}

#' @rdname svysem
#' @export
svycfa.default <- function(model, design,
                           estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_stop_design(design)
}

#' @rdname svysem
#' @export
svygrowth.survey.design2 <- function(model, design,
                                     estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::growth,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...),
    replicates = replicates,
    replicate.type = replicate.type,
    parallel = parallel
  )
}

#' @rdname svysem
#' @export
svygrowth.survey.design <- function(model, design,
                                    estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::growth,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...),
    replicates = replicates,
    replicate.type = replicate.type,
    parallel = parallel
  )
}

#' @rdname svysem
#' @export
svygrowth.svyrep.design <- function(model, design,
                                    estimator = c("MLM", "WLSMV", "MLR"),
                                    replicates = NULL,
                                    replicate.type = NULL,
                                    parallel = NULL, ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::growth,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...),
    replicates = replicates,
    replicate.type = replicate.type,
    parallel = parallel
  )
}

#' @rdname svysem
#' @export
svygrowth.default <- function(model, design,
                              estimator = c("MLM", "WLSMV", "MLR"),
                                  replicates = NULL,
                                  replicate.type = NULL,
                                  parallel = NULL, ...) {
  svysem_stop_design(design)
}

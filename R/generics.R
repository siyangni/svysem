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
#' results. In lavaan 0.6-21, [lavaan::parameterEstimates()] and
#' [lavaan::standardizedSolution()] omit standard errors unless `se = TRUE`
#' when the object's class is not exactly `lavaan`. [summary()] and the
#' printed display add the design degrees of freedom and the estimator.
#'
#' MLM uses complete data. `missing = "ml"` is an error and names MLR as the
#' estimator that will handle it. Rows that are incomplete on the model
#' variables are not dropped: their sampling weight is set to 0 inside the
#' design used for Gamma, which is how [subset()] represents a domain, and the
#' primary sampling unit structure is kept. Zero weights are legal. A negative
#' weight, or a design with no positive weight, is an error.
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
#' `conditional.x = TRUE` is not. One group of continuous variables only.
#' Clustering is taken from the design; do not pass `cluster=`, which would
#' double-count primary sampling units. `survey.lonely.psu` is read by
#' [survey::svymean()] and is never set here.
#'
#' @param model Model syntax, as for [lavaan::sem()], [lavaan::cfa()] or
#'   [lavaan::growth()].
#' @param design A `survey.design` or `survey.design2` object, including the
#'   result of [subset()].
#' @param estimator Variance estimator. Only `"MLM"` is implemented. `"WLSMV"`
#'   and `"MLR"` are reserved and raise `svysem_error_estimator`.
#' @param ... Passed to [lavaan::sem()], [lavaan::cfa()] or [lavaan::growth()].
#'   `data`, `sampling.weights`, sample moments and `NACOV` are taken from the
#'   design. `ov.order` is set internally so Gamma matches the model.
#' @return An S4 object of class `svysem` that contains the lavaan fit. Slots
#'   are `df` (design degrees of freedom), `estimator` (`"MLM"`),
#'   `variance_method` (`"linearization"`) and `svysem_call`.
#' @references Oberski D (2014). "lavaan.survey: An R Package for Complex
#'   Survey Analysis of Structural Equation Models." Journal of Statistical
#'   Software, 57(1), 1--27. \doi{10.18637/jss.v057.i01}
#' @seealso [svysem-class], [simulate_design()], [lavaan::sem()],
#'   [survey::svydesign()], [survey::svyglm()]
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
#' @export
svysem <- function(model, design, estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_prepare_dispatch(estimator, list(...))
  UseMethod("svysem", design)
}

#' @rdname svysem
#' @export
svycfa <- function(model, design, estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_prepare_dispatch(estimator, list(...))
  UseMethod("svycfa", design)
}

#' @rdname svysem
#' @export
svygrowth <- function(model, design, estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_prepare_dispatch(estimator, list(...))
  UseMethod("svygrowth", design)
}

#' @rdname svysem
#' @export
svysem.survey.design2 <- function(model, design,
                                  estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::sem,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...)
  )
}

#' @rdname svysem
#' @export
svysem.survey.design <- function(model, design,
                                 estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::sem,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...)
  )
}

#' @rdname svysem
#' @export
svysem.default <- function(model, design,
                           estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_stop_design(design)
}

#' @rdname svysem
#' @export
svycfa.survey.design2 <- function(model, design,
                                  estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::cfa,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...)
  )
}

#' @rdname svysem
#' @export
svycfa.survey.design <- function(model, design,
                                 estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::cfa,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...)
  )
}

#' @rdname svysem
#' @export
svycfa.default <- function(model, design,
                           estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_stop_design(design)
}

#' @rdname svysem
#' @export
svygrowth.survey.design2 <- function(model, design,
                                     estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::growth,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...)
  )
}

#' @rdname svysem
#' @export
svygrowth.survey.design <- function(model, design,
                                    estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_fit_design(
    model = model,
    design = design,
    fitter = lavaan::growth,
    svysem_call = match.call(),
    estimator = estimator,
    dots = list(...)
  )
}

#' @rdname svysem
#' @export
svygrowth.default <- function(model, design,
                              estimator = c("MLM", "WLSMV", "MLR"), ...) {
  svysem_stop_design(design)
}

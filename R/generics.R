#' Design-based structural equation models.
#'
#' `svysem()`, `svycfa()` and `svygrowth()` fit a lavaan model to a survey
#' design. This skeleton dispatches on the class of `design`. Estimators other
#' than MLM are not implemented yet.
#'
#' @param model Model syntax, as for [lavaan::sem()], [lavaan::cfa()] or
#'   [lavaan::growth()].
#' @param design A design object. `survey.design` and `survey.design2` are the
#'   linearization designs. `svyimputationList` is a set of completed designs.
#' @param estimator Variance estimator. `"MLM"` is implemented. `"WLSMV"` and
#'   `"MLR"` are reserved.
#' @param ... Passed to the lavaan fitter (`sem`, `cfa` or `growth`).
#' @return An object of class `svysem`, or `svysem_mi` for an imputation list.
#' @seealso [lavaan::sem()], [survey::svydesign()]
#' @references Oberski D (2014). "lavaan.survey: An R Package for Complex
#'   Survey Analysis of Structural Equation Models." Journal of Statistical
#'   Software, 57(1), 1--27. \doi{10.18637/jss.v057.i01}
#' @export
svysem <- function(model, design, estimator = c("MLM", "WLSMV", "MLR"), ...) {
  estimator <- match.arg(estimator)
  UseMethod("svysem")
}

#' @rdname svysem
#' @export
svycfa <- function(model, design, estimator = c("MLM", "WLSMV", "MLR"), ...) {
  estimator <- match.arg(estimator)
  UseMethod("svycfa")
}

#' @rdname svysem
#' @export
svygrowth <- function(model, design, estimator = c("MLM", "WLSMV", "MLR"), ...) {
  estimator <- match.arg(estimator)
  UseMethod("svygrowth")
}

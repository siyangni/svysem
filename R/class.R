#' Design-based SEM fit
#'
#' An object of class `svysem` is a [lavaan::lavaan-class] fit whose standard
#' errors and scaled test use the survey design. [coef()], [vcov()],
#' [lavaan::parameterEstimates()], [lavaan::standardizedSolution()],
#' [lavaan::lavTestWald()] and [lavaan::fitMeasures()] therefore report
#' design-based results with no further step.
#'
#' @slot df Design degrees of freedom, [survey::degf()] of the design passed
#'   to [svysem()], [svycfa()] or [svygrowth()].
#' @slot estimator Estimator used for the variance. `"MLM"` in this version.
#' @slot variance_method How the variance was estimated. `"linearization"`
#'   for a `survey.design` or `survey.design2`.
#' @slot svysem_call The call to [svysem()], [svycfa()] or [svygrowth()].
#' @return The class extends `lavaan`. See [svysem()] for the fitting functions.
#' @references Oberski D (2014). "lavaan.survey: An R Package for Complex
#'   Survey Analysis of Structural Equation Models." Journal of Statistical
#'   Software, 57(1), 1--27. \doi{10.18637/jss.v057.i01}
#' @seealso [svysem()], [survey::svydesign()], [lavaan::sem()]
#' @examples
#' design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 1)
#' fit <- svysem(
#'   "y ~ x1 + x2", design,
#'   fixed.x = FALSE, meanstructure = TRUE
#' )
#' fit@df
#' fit@estimator
#' fit@variance_method
#' @importClassesFrom lavaan lavaan
#' @export
setClass(
  "svysem",
  contains = "lavaan",
  slots = c(
    df = "numeric",
    estimator = "character",
    variance_method = "character",
    svysem_call = "call"
  ),
  validity = function(object) {
    if (length(object@df) != 1L || !is.finite(object@df)) {
      return("df must be a single finite design degrees of freedom.")
    }
    if (length(object@estimator) != 1L || object@estimator != "MLM") {
      return("estimator must be \"MLM\".")
    }
    if (length(object@variance_method) != 1L ||
          object@variance_method != "linearization") {
      return("variance_method must be \"linearization\".")
    }
    if (!is.call(object@svysem_call)) {
      return("svysem_call must be a call.")
    }
    TRUE
  }
)

# stats::coef and stats::vcov stay S3 generics unless lavaan is attached.
# The lavaan methods are S4, so call them on the parent object directly.
svysem_parent_method <- function(generic, object, ...) {
  method <- methods::selectMethod(generic, signature = "lavaan")
  method(methods::as(object, "lavaan"), ...)
}

#' @rdname svysem-class
#' @method coef svysem
#' @export
coef.svysem <- function(object, ...) {
  svysem_parent_method("coef", object, ...)
}

#' @rdname svysem-class
#' @method vcov svysem
#' @export
vcov.svysem <- function(object, ...) {
  svysem_parent_method("vcov", object, ...)
}

#' @rdname svysem-class
#' @export
setMethod("coef", "svysem", function(object, ...) {
  svysem_parent_method("coef", object, ...)
})

#' @rdname svysem-class
#' @export
setMethod("vcov", "svysem", function(object, ...) {
  svysem_parent_method("vcov", object, ...)
})

svysem_new <- function(fit, design, svysem_call) {
  methods::new(
    "svysem",
    fit,
    df = as.numeric(survey::degf(design)),
    estimator = "MLM",
    variance_method = "linearization",
    svysem_call = svysem_call
  )
}

#' @rdname svysem-class
#' @param object A `svysem` fit.
#' @export
setMethod("show", "svysem", function(object) {
  cat(
    "svysem fit\n",
    "Estimator: ", object@estimator, "\n",
    "Variance method: ", object@variance_method, "\n",
    "Design df: ", object@df, "\n",
    sep = ""
  )
  methods::show(methods::as(object, "lavaan"))
})

#' @rdname svysem-class
#' @param ... Passed to the lavaan summary.
#' @export
setMethod("summary", "svysem", function(object, ...) {
  ans <- summary(methods::as(object, "lavaan"), ...)
  ans$svysem <- list(
    df = object@df,
    estimator = object@estimator,
    variance_method = object@variance_method
  )
  class(ans) <- c("svysem.summary", class(ans))
  ans
})

#' @rdname svysem-class
#' @param x A summary of a `svysem` fit.
#' @export
print.svysem.summary <- function(x, ...) {
  meta <- x$svysem
  cat(
    "Estimator: ", meta$estimator, "\n",
    "Variance method: ", meta$variance_method, "\n",
    "Design degrees of freedom: ", meta$df, "\n",
    sep = ""
  )
  x$svysem <- NULL
  class(x) <- class(x)[class(x) != "svysem.summary"]
  print(x, ...)
  invisible(x)
}

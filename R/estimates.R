# lavaan drops se, z, pvalue, and the confidence interval when class(object)[1]
# is not "lavaan", unless the caller passes se = TRUE. svysem stays the first
# class so summary() and show() dispatch. These methods pass se = TRUE.

#' Parameter estimates
#'
#' Design-based standard errors and confidence intervals for a [svysem-class]
#' fit.
#'
#' [lavaan::parameterEstimates()] omits standard errors, z statistics,
#' p-values, and confidence intervals when the object's first class is not
#' `"lavaan"`, unless `se = TRUE` is supplied. A svysem fit keeps the class
#' `svysem` so that [summary()] and [show()] can report the design. This
#' method supplies `se = TRUE`, so the default table includes `se`, `z`,
#' `pvalue`, `ci.lower`, and `ci.upper`. Those standard errors match
#' `sqrt(diag(vcov(object)))`. `se = FALSE` returns the point estimates only.
#'
#' `lavaan::parameterEstimates()` still applies lavaan's class check. After
#' `library(lavaan)` is called later than `library(svysem)`, call
#' `svysem::parameterEstimates()`.
#'
#' @param object A fitted model. A `svysem` fit uses the method below. Any
#'   other object is passed to [lavaan::parameterEstimates()].
#' @param se Include standard errors. The default is `TRUE`.
#' @param ... Passed to [lavaan::parameterEstimates()]. `zstat`, `pvalue`, and
#'   `ci` keep lavaan's defaults.
#' @return The data frame from [lavaan::parameterEstimates()].
#' @seealso [standardizedSolution()][standardizedSolution],
#'   [lavaan::parameterEstimates()]
#' @export
parameterEstimates <- function(object, ...) {
  UseMethod("parameterEstimates")
}

#' @rdname parameterEstimates
#' @export
parameterEstimates.svysem <- function(object, se = TRUE, ...) {
  lavaan::parameterEstimates(object, se = se, ...)
}

#' @rdname parameterEstimates
#' @export
parameterEstimates.default <- function(object, ...) {
  lavaan::parameterEstimates(object, ...)
}

#' Standardized solution
#'
#' Design-based standard errors for the standardized estimates of a
#' [svysem-class] fit.
#'
#' [lavaan::standardizedSolution()] uses the same class check as
#' [lavaan::parameterEstimates()]. This method supplies `se = TRUE` for a
#' svysem fit. `se = FALSE` returns the standardized point estimates only.
#' After `library(lavaan)` is called later than `library(svysem)`, call
#' `svysem::standardizedSolution()`.
#'
#' @param object A fitted model. A `svysem` fit uses the method below. Any
#'   other object is passed to [lavaan::standardizedSolution()].
#' @param se Include standard errors. The default is `TRUE`.
#' @param ... Passed to [lavaan::standardizedSolution()].
#' @return The data frame from [lavaan::standardizedSolution()].
#' @seealso [parameterEstimates()][parameterEstimates],
#'   [lavaan::standardizedSolution()]
#' @export
standardizedSolution <- function(object, ...) {
  UseMethod("standardizedSolution")
}

#' @rdname standardizedSolution
#' @export
standardizedSolution.svysem <- function(object, se = TRUE, ...) {
  lavaan::standardizedSolution(object, se = se, ...)
}

#' @rdname standardizedSolution
#' @export
standardizedSolution.default <- function(object, ...) {
  lavaan::standardizedSolution(object, ...)
}

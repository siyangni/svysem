#' Multiply imputed design-based SEM
#'
#' Fits of one model to each completed survey design, pooled afterwards with
#' Rubin's rules. The design degrees of freedom are the complete-data degrees
#' of freedom (`df.com`) unless the user overrides them in [pooled_estimates()]
#' or [pooled_wald()].
#'
#' The object is returned by [svysem()], [svycfa()] or [svygrowth()] when
#' `design` is a `svyimputationList`. Each completed design is fitted
#' separately. `on_fail = "stop"` stops at the first imputation that errors.
#' `on_fail = "omit"` records that imputation in `status` and pools the rest,
#' provided at least two fits remain.
#'
#' @aliases svysem_mi
#' @slot fits Successful fits, in imputation order, with failures removed.
#' @slot df Design degrees of freedom, [survey::degf()] of the completed
#'   designs. This is the default `df.com` when pooling.
#' @slot nobs_used Number of rows with a positive sampling weight in the
#'   first completed design.
#' @slot status One row per imputation, including failures. Columns are
#'   `no.conv`, `no.se`, `no.npd` and `message` (an empty string when the
#'   imputation did not error).
#' @slot estimator Estimator name. Only `"MLM"` is accepted.
#' @slot svysem_call The call that created the object.
#' @seealso [pooled_estimates()], [pooled_wald()], [survey::svydesign()]
#' @references
#' Rubin D (1987). Multiple Imputation for Nonresponse in Surveys.
#' Wiley, New York. \doi{10.1002/9780470316696}
#'
#' Barnard J, Rubin DB (1999). "Small-Sample Degrees of Freedom with
#' Multiple Imputation." Biometrika, 86(4), 948--959.
#' \doi{10.1093/biomet/86.4.948}
#' @examples
#' imp <- list(
#'   data.frame(
#'     y = 1:30, x = rep(1:10, 3), w = 1,
#'     id = 1:30, strata = rep(1:2, each = 15)
#'   ),
#'   data.frame(
#'     y = 1:30 + 0.2, x = rep(1:10, 3), w = 1,
#'     id = 1:30, strata = rep(1:2, each = 15)
#'   )
#' )
#' des <- survey::svydesign(
#'   ids = ~id, strata = ~strata, weights = ~w,
#'   data = mitools::imputationList(imp)
#' )
#' stopifnot(inherits(des, "svyimputationList"))
#' # svysem("y ~ x", des) returns a svysem_mi object once each completed
#' # design can be fitted. Pooling itself only needs coef() and vcov():
#' fit_i <- function(est, v) {
#'   dimnames(v) <- list(names(est), names(est))
#'   structure(list(coefficients = est, variance = v), class = "svysem_eg")
#' }
#' registerS3method(
#'   "coef", "svysem_eg", function(object, ...) object$coefficients,
#'   envir = asNamespace("stats")
#' )
#' registerS3method(
#'   "vcov", "svysem_eg", function(object, ...) object$variance,
#'   envir = asNamespace("stats")
#' )
#' fits <- list(
#'   fit_i(c(yx = 0.4), matrix(0.02)),
#'   fit_i(c(yx = 0.6), matrix(0.03))
#' )
#' pooled_estimates(fits, df.com = survey::degf(des$designs[[1]]))
#' @export
setClass(
  "svysem_mi",
  slots = c(
    fits = "list",
    df = "numeric",
    nobs_used = "numeric",
    status = "data.frame",
    estimator = "character",
    svysem_call = "call"
  )
)

# Non-lavaan objects, including test stubs, count as converged. Standard
# errors are available when vcov() is a finite matrix. A latent covariance
# is not positive definite when an eigenvalue is below -1e-8.
svysem_latent_npd <- function(fit) {
  cv <- tryCatch(lavaan::lavInspect(fit, "cov.lv"), error = function(e) NULL)
  if (is.matrix(cv)) {
    mats <- list(cv)
  } else if (is.list(cv)) {
    mats <- cv
  } else {
    mats <- list()
  }
  for (mat in mats) {
    if (!is.matrix(mat) || !nrow(mat)) {
      next
    }
    ev <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < -1e-8)) {
      return(TRUE)
    }
  }
  FALSE
}

svysem_fit_status <- function(fit) {
  vc <- tryCatch(as.matrix(svysem_vcov(fit)), error = function(e) NULL)
  no_se <- is.null(vc) || !all(is.finite(vc))
  if (inherits(fit, "lavaan")) {
    converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
    no_npd <- isTRUE(svysem_latent_npd(fit))
    no_conv <- !converged
  } else {
    no_conv <- FALSE
    no_npd <- FALSE
  }
  data.frame(
    no.conv = no_conv,
    no.se = no_se,
    no.npd = no_npd,
    message = ""
  )
}

svysem_blank_status <- function(m) {
  data.frame(
    no.conv = rep(FALSE, m),
    no.se = rep(FALSE, m),
    no.npd = rep(FALSE, m),
    message = rep("", m)
  )
}

svysem_impute <- function(model, design, estimator, on_fail, parallel, fitter,
                          svysem_call, ...) {
  estimator <- match.arg(estimator, c("MLM", "WLSMV", "MLR"))
  on_fail <- match.arg(on_fail, c("stop", "omit"))
  if (!identical(estimator, "MLM")) {
    svysem_stop(
      sprintf(
        "Estimator \"%s\" is not implemented. Only \"MLM\" is available.",
        estimator
      ),
      "svysem_error_estimator"
    )
  }
  if (!inherits(design, "svyimputationList") || !is.list(design$designs)) {
    svysem_stop(
      "design must be a svyimputationList with a list of completed designs.",
      "svysem_error_design"
    )
  }
  designs <- design$designs
  m <- length(designs)
  if (m < 2L) {
    svysem_stop(
      "At least two imputations are required for pooling.",
      "svysem_error_pooling"
    )
  }
  dfs <- vapply(designs, function(d) as.numeric(survey::degf(d)), numeric(1))
  if (any(abs(dfs - dfs[[1]]) > 1e-8)) {
    svysem_stop(
      "Completed designs do not share the same degrees of freedom.",
      "svysem_error_design"
    )
  }
  w <- as.numeric(stats::weights(designs[[1]]))
  nobs_used <- sum(w > 0, na.rm = TRUE)
  if (is.null(parallel)) {
    loop <- lapply
  } else if (is.function(parallel)) {
    loop <- parallel
  } else {
    svysem_stop(
      "parallel must be NULL or a function.",
      "svysem_error_design"
    )
  }
  results <- loop(seq_len(m), function(i) {
    tryCatch(
      {
        fit <- fitter(
          model = model, design = designs[[i]], estimator = estimator, ...
        )
        list(ok = TRUE, fit = fit, status = svysem_fit_status(fit))
      },
      error = function(e) {
        list(ok = FALSE, message = conditionMessage(e))
      }
    )
  })
  if (length(results) != m) {
    svysem_stop(
      "parallel did not return one result per imputation.",
      "svysem_error_imputation"
    )
  }
  status <- svysem_blank_status(m)
  fits <- vector("list", m)
  ok <- logical(m)
  for (i in seq_len(m)) {
    res <- results[[i]]
    failed <- inherits(res, "try-error") || !is.list(res) || !isTRUE(res$ok)
    if (failed) {
      msg <- if (inherits(res, "try-error")) {
        as.character(res)[[1]]
      } else if (is.list(res) && !is.null(res$message)) {
        paste(res$message, collapse = " ")
      } else {
        "the imputation did not return a fit"
      }
      ok[[i]] <- FALSE
      status$no.conv[[i]] <- TRUE
      status$no.se[[i]] <- TRUE
      status$no.npd[[i]] <- FALSE
      status$message[[i]] <- msg
    } else {
      ok[[i]] <- TRUE
      fits[[i]] <- res$fit
      status$no.conv[[i]] <- res$status$no.conv
      status$no.se[[i]] <- res$status$no.se
      status$no.npd[[i]] <- res$status$no.npd
      status$message[[i]] <- res$status$message
    }
  }
  if (any(!ok) && identical(on_fail, "stop")) {
    i <- which(!ok)[[1]]
    svysem_stop(
      sprintf("Imputation %d failed: %s", i, status$message[[i]]),
      "svysem_error_imputation"
    )
  }
  fits <- fits[ok]
  if (length(fits) < 2L) {
    svysem_stop(
      "Fewer than two imputations produced a fit, so there is nothing to pool.",
      "svysem_error_pooling"
    )
  }
  methods::new(
    "svysem_mi",
    fits = fits,
    df = dfs[[1]],
    nobs_used = as.numeric(nobs_used),
    status = status,
    estimator = estimator,
    svysem_call = svysem_call
  )
}

#' Fit each completed design
#'
#' [svysem()], [svycfa()] and [svygrowth()] dispatch on a
#' `svyimputationList`. Every completed design is passed to the same generic,
#' which selects [lavaan::sem()], [lavaan::cfa()] or [lavaan::growth()].
#' Only `estimator = "MLM"` is accepted. `parallel` is any function of
#' `(X, FUN)` with the same contract as [lapply()]; `NULL` uses [lapply()].
#' Tests and examples use at most two cores.
#'
#' @param model Model syntax, as for [lavaan::sem()], [lavaan::cfa()] or
#'   [lavaan::growth()].
#' @param design A `svyimputationList`, from
#'   `svydesign(data = mitools::imputationList(...))`.
#' @param estimator Variance estimator. `"MLM"` is the only implemented value.
#' @param ... Passed to the fitter for each completed design.
#' @param on_fail `"stop"` (default) or `"omit"`. `"omit"` drops a failed
#'   imputation, records it in the status table, and continues.
#' @param parallel `NULL` to use [lapply()], or a function called as
#'   `parallel(X, FUN)`. One example is
#'   `function(X, FUN) parallel::mclapply(X, FUN, mc.cores = 2)`.
#' @return A [svysem_mi-class] object. Its `df` slot is the design degrees
#'   of freedom and is the default `df.com` for pooling.
#' @seealso [pooled_estimates()], [pooled_wald()], [survey::svydesign()]
#' @references
#' Rubin D (1987). Multiple Imputation for Nonresponse in Surveys.
#' Wiley, New York. \doi{10.1002/9780470316696}
#'
#' Barnard J, Rubin DB (1999). "Small-Sample Degrees of Freedom with
#' Multiple Imputation." Biometrika, 86(4), 948--959.
#' \doi{10.1093/biomet/86.4.948}
#' @examples
#' imp <- list(
#'   data.frame(
#'     y = 1:30, x = rep(1:10, 3), w = 1,
#'     id = 1:30, strata = rep(1:2, each = 15)
#'   ),
#'   data.frame(
#'     y = 1:30 + 0.2, x = rep(1:10, 3), w = 1,
#'     id = 1:30, strata = rep(1:2, each = 15)
#'   )
#' )
#' des <- survey::svydesign(
#'   ids = ~id, strata = ~strata, weights = ~w,
#'   data = mitools::imputationList(imp)
#' )
#' stopifnot(inherits(des, "svyimputationList"), length(des$designs) == 2L)
#' @export
svysem.svyimputationList <- function(model, design,
                                      estimator = c("MLM", "WLSMV", "MLR"), ...,
                                      on_fail = c("stop", "omit"),
                                      parallel = NULL) {
  svysem_impute(
    model = model, design = design, estimator = estimator,
    on_fail = on_fail, parallel = parallel, fitter = svysem,
    svysem_call = match.call(), ...
  )
}

#' @rdname svysem.svyimputationList
#' @export
svycfa.svyimputationList <- function(model, design,
                                      estimator = c("MLM", "WLSMV", "MLR"), ...,
                                      on_fail = c("stop", "omit"),
                                      parallel = NULL) {
  svysem_impute(
    model = model, design = design, estimator = estimator,
    on_fail = on_fail, parallel = parallel, fitter = svycfa,
    svysem_call = match.call(), ...
  )
}

#' @rdname svysem.svyimputationList
#' @export
svygrowth.svyimputationList <- function(model, design,
                                         estimator = c("MLM", "WLSMV", "MLR"), ...,
                                         on_fail = c("stop", "omit"),
                                         parallel = NULL) {
  svysem_impute(
    model = model, design = design, estimator = estimator,
    on_fail = on_fail, parallel = parallel, fitter = svygrowth,
    svysem_call = match.call(), ...
  )
}

#' @rdname svysem_mi-class
#' @export
setMethod("show", "svysem_mi", function(object) {
  st <- object@status
  cat("svysem_mi (multiple imputation)\n")
  cat("Estimator: ", object@estimator, "\n", sep = "")
  cat(
    "Imputations: ", length(object@fits), " used / ", nrow(st), " supplied\n",
    sep = ""
  )
  cat("Design degrees of freedom: ", format(object@df), "\n", sep = "")
  cat("Observations: ", format(object@nobs_used), "\n", sep = "")
  cat(
    "Not converged: ", sum(st$no.conv),
    "; no SE: ", sum(st$no.se),
    "; not positive definite: ", sum(st$no.npd),
    "\n",
    sep = ""
  )
  invisible(object)
})

# Exported pooling. coef() on a lavaan object has to go through lavaan::coef:
# stats::coef is the S3 generic in this namespace and has no lavaan method.

svysem_coef <- function(fit) {
  if (inherits(fit, "lavaan")) {
    lavaan::coef(fit)
  } else {
    stats::coef(fit)
  }
}

svysem_vcov <- function(fit) {
  if (inherits(fit, "lavaan")) {
    lavaan::vcov(fit)
  } else {
    stats::vcov(fit)
  }
}

svysem_collect <- function(object, df.com, supplied) {
  if (inherits(object, "svysem_mi")) {
    if (!supplied) {
      df.com <- object@df
    }
    fits <- object@fits
  } else if (is.list(object) && !is.data.frame(object)) {
    fits <- object
  } else {
    svysem_stop(
      "object must be a list of fits or a svysem_mi object.",
      "svysem_error_pooling"
    )
  }
  if (length(fits) < 1L) {
    svysem_stop("object contains no fits.", "svysem_error_pooling")
  }
  list(fits = fits, df.com = df.com)
}

svysem_qhat_uhat <- function(fits) {
  coefs <- lapply(fits, function(fit) {
    b <- svysem_coef(fit)
    if (is.null(names(b)) || anyDuplicated(names(b)) || !length(b)) {
      svysem_stop(
        "Coefficients must be a non-empty named vector with unique names.",
        "svysem_error_pooling"
      )
    }
    b
  })
  nms <- names(coefs[[1]])
  if (!all(vapply(coefs, function(b) setequal(names(b), nms), logical(1)))) {
    svysem_stop(
      "Parameter names differ across imputations.",
      "svysem_error_pooling"
    )
  }
  p <- length(nms)
  m <- length(fits)
  Qhat <- matrix(NA_real_, p, m, dimnames = list(nms, NULL))
  Uhat <- array(NA_real_, c(p, p, m), dimnames = list(nms, nms, NULL))
  for (i in seq_len(m)) {
    Qhat[, i] <- as.numeric(coefs[[i]][nms])
    V <- tryCatch(as.matrix(svysem_vcov(fits[[i]])), error = function(e) NULL)
    if (is.null(V)) {
      svysem_stop(
        sprintf("vcov() failed for imputation %d.", i),
        "svysem_error_pooling"
      )
    }
    if (!is.null(rownames(V))) {
      if (!all(nms %in% rownames(V))) {
        svysem_stop(
          "vcov() names do not match coef() names.",
          "svysem_error_pooling"
        )
      }
      V <- V[nms, nms, drop = FALSE]
    } else if (!all(dim(V) == c(p, p))) {
      svysem_stop(
        "vcov() does not have one row per coefficient.",
        "svysem_error_pooling"
      )
    }
    Uhat[, , i] <- V
  }
  list(Qhat = Qhat, Uhat = Uhat)
}

svysem_pool_frame <- function(pool, level, nms) {
  tcrit <- stats::qt(1 - (1 - level) / 2, df = pool$df)
  data.frame(
    est = as.numeric(pool$est),
    se = as.numeric(pool$se),
    t = as.numeric(pool$t),
    df = as.numeric(pool$df),
    p = as.numeric(pool$p),
    ci.lower = as.numeric(pool$est - tcrit * pool$se),
    ci.upper = as.numeric(pool$est + tcrit * pool$se),
    riv = as.numeric(pool$riv),
    lambda = as.numeric(pool$lambda),
    fmi = as.numeric(pool$fmi),
    m = pool$m,
    row.names = nms,
    check.names = FALSE
  )
}

svysem_check_level <- function(level) {
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    svysem_stop(
      "level must be one number strictly between 0 and 1.",
      "svysem_error_pooling"
    )
  }
  level
}

svysem_std_piece <- function(fit) {
  # lavaan 0.6-21 omits se unless se=TRUE when class[1] is not "lavaan".
  ss <- lavaan::standardizedSolution(fit, type = "std.all", se = TRUE)
  needed <- c("lhs", "op", "rhs", "est.std", "se")
  if (!all(needed %in% names(ss))) {
    svysem_stop(
      "standardizedSolution() did not return est.std and se.",
      "svysem_error_pooling"
    )
  }
  key <- paste0(ss$lhs, ss$op, ss$rhs)
  if ("group" %in% names(ss) && length(unique(ss$group)) > 1L) {
    key <- paste0(key, ".g", ss$group)
  }
  if (anyDuplicated(key)) {
    svysem_stop(
      "Standardized parameter names are not unique.",
      "svysem_error_pooling"
    )
  }
  data.frame(
    key = key,
    est = ss$est.std,
    se = ss$se,
    fisher = ss$op == "~~" & ss$lhs != ss$rhs,
    ok = is.finite(ss$est.std) & is.finite(ss$se) & ss$se > 0
  )
}

# Correlations are pooled as z = atanh(r). est, se and the confidence
# limits are transformed back; t, df, p, riv, lambda and fmi stay on z.
svysem_pool_standardized <- function(fits, level, df.com) {
  if (!all(vapply(fits, inherits, logical(1), what = "lavaan"))) {
    svysem_stop(
      "standardized = TRUE requires every fit to be a lavaan object.",
      "svysem_error_pooling"
    )
  }
  pieces <- lapply(fits, svysem_std_piece)
  all_keys <- lapply(pieces, `[[`, "key")
  if (!all(vapply(all_keys, function(k) identical(k, all_keys[[1]]), logical(1)))) {
    svysem_stop(
      "Standardized parameters differ across imputations.",
      "svysem_error_pooling"
    )
  }
  keys <- Reduce(intersect, lapply(pieces, function(piece) piece$key[piece$ok]))
  if (!length(keys)) {
    svysem_stop(
      "No standardized estimate has a positive standard error.",
      "svysem_error_pooling"
    )
  }
  first <- pieces[[1]]
  fisher <- unname(first$fisher[match(keys, first$key)])
  p <- length(keys)
  m <- length(fits)
  Qhat <- matrix(NA_real_, p, m, dimnames = list(keys, NULL))
  Uhat <- array(0, c(p, p, m), dimnames = list(keys, keys, NULL))
  for (i in seq_len(m)) {
    piece <- pieces[[i]]
    idx <- match(keys, piece$key)
    est <- piece$est[idx]
    se <- piece$se[idx]
    if (any(fisher & abs(est) >= 1)) {
      svysem_stop(
        "A standardized correlation is outside (-1, 1) and cannot be Fisher-transformed.",
        "svysem_error_pooling"
      )
    }
    z <- est
    se_z <- se
    z[fisher] <- atanh(est[fisher])
    se_z[fisher] <- se[fisher] / (1 - est[fisher]^2)
    Qhat[, i] <- z
    diag(Uhat[, , i]) <- se_z^2
  }
  pool <- rubin_pool(Qhat, Uhat, df.com = df.com)
  out <- svysem_pool_frame(pool, level, keys)
  if (any(fisher)) {
    z_bar <- out$est[fisher]
    se_z <- out$se[fisher]
    tcrit <- stats::qt(1 - (1 - level) / 2, df = out$df[fisher])
    r <- tanh(z_bar)
    out$est[fisher] <- r
    out$se[fisher] <- se_z * (1 - r^2)
    out$ci.lower[fisher] <- tanh(z_bar - tcrit * se_z)
    out$ci.upper[fisher] <- tanh(z_bar + tcrit * se_z)
  }
  out
}

#' Pooled parameter estimates
#'
#' Combine point estimates and covariance matrices from a list of completed-data
#' fits with Rubin's rules. Degrees of freedom use the Barnard-Rubin adjustment
#' when `df.com` is finite. For a [svysem_mi-class] object, `df.com` defaults to
#' the design degrees of freedom.
#'
#' `t` is the t statistic, the pooled estimate divided by its standard error,
#' not the total variance. With `standardized = TRUE`, each fit must be a
#' lavaan object. Standard errors on that table are requested with
#' `se = TRUE`, because lavaan 0.6-21 omits them when the fit's class is not
#' exactly `lavaan`. Loadings, regressions and variances are pooled on the
#' `std.all` scale. Correlations (`~~` with different left- and right-hand
#' sides) are pooled on the Fisher z scale, `z = atanh(r)` with
#' `se_z = se_r / (1 - r^2)`. Reported `est`, `se` and confidence limits are
#' transformed back (`est = tanh(z)`, `se = se_z * (1 - est^2)`, limits
#' `tanh(z +/- t * se_z)`). For those rows, `t`, `df`, `p`, `riv`, `lambda`
#' and `fmi` stay on the z scale.
#'
#' @param object A list of objects with [stats::coef()] and [stats::vcov()],
#'   or a [svysem_mi-class] object.
#' @param standardized If `TRUE`, pool `std.all` standardized estimates.
#'   Correlations use the Fisher z scale.
#' @param level Confidence level for `ci.lower` and `ci.upper`.
#' @param df.com Complete-data degrees of freedom. For a `svysem_mi` object
#'   the default is the design degrees of freedom. `NULL` or `Inf` uses the
#'   large-sample Rubin degrees of freedom.
#' @return A data frame with one row per parameter and columns `est`, `se`,
#'   `t`, `df`, `p`, `ci.lower`, `ci.upper`, `riv`, `lambda`, `fmi` and `m`.
#' @seealso [pooled_wald()], [svysem()]
#' @references
#' Rubin D (1987). Multiple Imputation for Nonresponse in Surveys.
#' Wiley, New York. \doi{10.1002/9780470316696}
#'
#' Barnard J, Rubin DB (1999). "Small-Sample Degrees of Freedom with
#' Multiple Imputation." Biometrika, 86(4), 948--959.
#' \doi{10.1093/biomet/86.4.948}
#' @examples
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
#'   fit_i(c(b1 = 0.30, b2 = 0.10), diag(2) * 0.02),
#'   fit_i(c(b1 = 0.50, b2 = 0.20), diag(2) * 0.03),
#'   fit_i(c(b1 = 0.40, b2 = 0.15), diag(2) * 0.02),
#'   fit_i(c(b1 = 0.35, b2 = 0.12), diag(2) * 0.025)
#' )
#' pooled_estimates(fits, df.com = 80)
#' @export
#' @importFrom lavaan fitMeasures
pooled_estimates <- function(object, standardized = FALSE, level = 0.95,
                             df.com = NULL) {
  level <- svysem_check_level(level)
  if (!isTRUE(standardized) && !identical(standardized, FALSE)) {
    svysem_stop(
      "standardized must be TRUE or FALSE.",
      "svysem_error_pooling"
    )
  }
  picked <- svysem_collect(object, df.com, supplied = !missing(df.com))
  if (isTRUE(standardized)) {
    return(svysem_pool_standardized(picked$fits, level, picked$df.com))
  }
  qu <- svysem_qhat_uhat(picked$fits)
  pool <- rubin_pool(qu$Qhat, qu$Uhat, df.com = picked$df.com)
  svysem_pool_frame(pool, level, rownames(qu$Qhat))
}

#' Wald test of pooled constraints
#'
#' Combine a set of constraints with the D1 rule. A constraint is an R
#' expression in the free-parameter names that is tested against zero.
#' `a == b` is the scalar `a - b`. Labels and `:=` names from a lavaan fit
#' are rewritten onto the free parameters and differentiated with the delta
#' method inside each imputation. Inequality constraints are not supported.
#'
#' For a [svysem_mi-class] object, `df.com` defaults to the design degrees
#' of freedom. The small-sample denominator follows mitml when that value is
#' positive. When it is undefined (`k * (m - 1) <= 4`) or negative, one
#' constraint uses the scalar Barnard-Rubin degrees of freedom and several
#' constraints use a positive harmonic mean of the large-sample denominator
#' and the observed-data degrees of freedom.
#'
#' @param object A list of objects with [stats::coef()] and [stats::vcov()],
#'   or a [svysem_mi-class] object.
#' @param constraints Character vector of constraints. Lavaan names that are
#'   not syntactic are written in backticks, for example `` "`y~x` - `y~z`" ``.
#' @param df.com Complete-data degrees of freedom. For a `svysem_mi` object
#'   the default is the design degrees of freedom. `NULL` uses the
#'   large-sample D1 denominator.
#' @return A one-row data frame with columns `F`, `df1`, `df2`, `pvalue`
#'   and `riv`.
#' @seealso [pooled_estimates()], [svysem()]
#' @references
#' Rubin D (1987). Multiple Imputation for Nonresponse in Surveys.
#' Wiley, New York. \doi{10.1002/9780470316696}
#'
#' Barnard J, Rubin DB (1999). "Small-Sample Degrees of Freedom with
#' Multiple Imputation." Biometrika, 86(4), 948--959.
#' \doi{10.1093/biomet/86.4.948}
#' @examples
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
#'   fit_i(c(b1 = 0.30, b2 = 0.10), diag(2) * 0.02),
#'   fit_i(c(b1 = 0.50, b2 = 0.20), diag(2) * 0.03),
#'   fit_i(c(b1 = 0.40, b2 = 0.15), diag(2) * 0.02),
#'   fit_i(c(b1 = 0.35, b2 = 0.12), diag(2) * 0.025),
#'   fit_i(c(b1 = 0.45, b2 = 0.18), diag(2) * 0.02),
#'   fit_i(c(b1 = 0.38, b2 = 0.11), diag(2) * 0.022)
#' )
#' pooled_wald(fits, constraints = "b1 - b2", df.com = 80)
#' @export
pooled_wald <- function(object, constraints, df.com = NULL) {
  picked <- svysem_collect(object, df.com, supplied = !missing(df.com))
  qu <- svysem_qhat_uhat(picked$fits)
  cu <- svysem_constraint_moments(picked$fits, constraints, qu)
  d1 <- rubin_d1(cu$Qhat, cu$Uhat, df.com = picked$df.com)
  data.frame(
    F = d1$F, df1 = d1$df1, df2 = d1$df2,
    pvalue = d1$pvalue, riv = d1$riv
  )
}

#' @rdname svysem_mi-class
#' @export
coef.svysem_mi <- function(object, ...) {
  pe <- pooled_estimates(object)
  stats::setNames(pe$est, rownames(pe))
}

#' @rdname svysem_mi-class
#' @export
setMethod("coef", "svysem_mi", function(object, ...) {
  coef.svysem_mi(object, ...)
})

#' @rdname svysem_mi-class
#' @export
vcov.svysem_mi <- function(object, ...) {
  qu <- svysem_qhat_uhat(object@fits)
  rubin_pool(qu$Qhat, qu$Uhat, df.com = NULL)$total
}

#' @rdname svysem_mi-class
#' @export
setMethod("vcov", "svysem_mi", function(object, ...) {
  vcov.svysem_mi(object, ...)
})

#' @rdname svysem_mi-class
#' @export
nobs.svysem_mi <- function(object, ...) {
  object@nobs_used
}

#' @rdname svysem_mi-class
#' @export
summary.svysem_mi <- function(object, ...) {
  st <- object@status
  cat("Pooled estimates (Rubin's rules)\n")
  cat(
    "Complete-data degrees of freedom (design df): ",
    format(object@df), "\n",
    sep = ""
  )
  cat("Imputations: ", nrow(st), "\n", sep = "")
  cat("  no.conv: ", sum(st$no.conv), "\n", sep = "")
  cat("  no.se: ", sum(st$no.se), "\n", sep = "")
  cat("  no.npd: ", sum(st$no.npd), "\n", sep = "")
  cat("  failed: ", sum(nzchar(st$message)), "\n", sep = "")
  est <- pooled_estimates(object)
  print(est)
  invisible(est)
}

svysem_fit_measure_table <- function(object, fm.args) {
  wanted <- c(
    "chisq.scaled", "df.scaled",
    "cfi.robust", "tli.robust", "rmsea.robust", "srmr"
  )
  fms <- lapply(object@fits, function(fit) {
    if (!inherits(fit, "lavaan")) {
      svysem_stop(
        "fitMeasures() requires lavaan fits.",
        "svysem_error_pooling"
      )
    }
    lavaan::fitMeasures(fit, fit.measures = wanted, fm.args = fm.args)
  })
  if (!"chisq.scaled" %in% names(fms[[1]]) || !"df.scaled" %in% names(fms[[1]])) {
    svysem_stop(
      "The fits do not report a scaled chi-square.",
      "svysem_error_pooling"
    )
  }
  dfs <- vapply(fms, function(fm) unname(fm[["df.scaled"]]), numeric(1))
  if (any(abs(dfs - dfs[[1]]) > 1e-8)) {
    svysem_stop(
      "Chi-square degrees of freedom differ across imputations.",
      "svysem_error_pooling"
    )
  }
  w <- vapply(fms, function(fm) unname(fm[["chisq.scaled"]]), numeric(1))
  d2 <- rubin_d2(w, DF = dfs[[1]], asymptotic = TRUE)
  out <- c(
    chisq.scaled = unname(d2$chisq),
    df.scaled = unname(d2$df),
    pvalue.scaled = unname(d2$pvalue),
    ariv = unname(d2$ariv),
    fmi = unname(d2$fmi)
  )
  for (nm in c("cfi.robust", "tli.robust", "rmsea.robust", "srmr")) {
    if (!nm %in% names(fms[[1]])) {
      next
    }
    vals <- vapply(fms, function(fm) unname(fm[[nm]]), numeric(1))
    out[[paste0(nm, ".mean")]] <- mean(vals)
    out[[paste0(nm, ".min")]] <- min(vals)
    out[[paste0(nm, ".max")]] <- max(vals)
  }
  out
}

#' @param object A `svysem_mi` object.
#' @param fit.measures `"all"` or a character vector of measure names to keep.
#' @param baseline.model Accepted for compatibility with lavaan's generic.
#'   Not used: each imputation already carries its own baseline.
#' @param h1.model Accepted for compatibility with lavaan's generic. Not used.
#' @param fm.args Passed to [lavaan::fitMeasures()] on each imputation.
#' @param output `"vector"` (default), or `"matrix"` / `"data.frame"` for a
#'   one-row data frame.
#' @param ... Further arguments, ignored.
#' @rdname svysem_mi-class
#' @export
setMethod(
  "fitMeasures",
  signature(object = "svysem_mi"),
  function(object, fit.measures = "all", baseline.model = NULL, h1.model = NULL,
           fm.args = list(standard.test = "default", scaled.test = "default",
                          rmsea.ci.level = 0.9, rmsea.close.h0 = 0.05,
                          rmsea.notclose.h0 = 0.08, robust = TRUE,
                          cat.check.pd = TRUE),
           output = "vector", ...) {
    out <- svysem_fit_measure_table(object, fm.args)
    if (!identical(fit.measures, "all")) {
      out <- out[intersect(fit.measures, names(out))]
    }
    if (identical(output, "matrix") || identical(output, "data.frame")) {
      return(as.data.frame(as.list(out), check.names = FALSE))
    }
    out
  }
)

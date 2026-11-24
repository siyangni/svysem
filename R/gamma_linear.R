# ML second-moment names in lavaan's vech order. Products are not scaled by
# n / (n - 1); svyvar() is the wrong estimand for this Gamma.
svysem_product_names <- function(ov) {
  p <- length(ov)
  row_idx <- lavaan::lav_matrix_vech_row_idx(p)
  col_idx <- lavaan::lav_matrix_vech_col_idx(p)
  list(
    row = row_idx,
    col = col_idx,
    names = paste0(ov[col_idx], "~~", ov[row_idx])
  )
}

svysem_canonical_names <- function(ov) {
  products <- svysem_product_names(ov)
  means <- paste0(ov, "~1")
  list(full = c(means, products$names), cov = products$names, index = products)
}

# Casewise means (raw y) and centered products. Non-finite rows, which already
# have weight 0, are filled with 0 so they do not poison the PSU totals.
svysem_moment_matrix <- function(design, ov) {
  observed <- as.matrix(design$variables[, ov, drop = FALSE])
  storage.mode(observed) <- "double"
  weights <- as.numeric(weights(design, "sampling"))
  bad <- !stats::complete.cases(observed) | !is.finite(weights) | weights == 0
  observed[!is.finite(observed)] <- 0

  filled <- design
  variables <- data.frame(filled$variables, check.names = FALSE, stringsAsFactors = FALSE)
  variables[, ov] <- observed
  variables[bad, ov] <- 0
  filled$variables <- variables
  center <- stats::coef(survey::svymean(stats::reformulate(ov), filled))
  center <- center[ov]

  centered <- sweep(observed, 2L, center, "-")
  centered[bad, ] <- 0
  index <- svysem_product_names(ov)
  products <- centered[, index$row, drop = FALSE] * centered[, index$col, drop = FALSE]
  means <- observed
  means[bad, ] <- 0
  moments <- cbind(means, products)
  colnames(moments) <- c(paste0(ov, "~1"), index$names)
  moments
}

svysem_plain_matrix <- function(x) {
  x <- as.matrix(x)
  x <- (x + t(x)) / 2
  storage.mode(x) <- "double"
  k <- nrow(x)
  attributes(x) <- list(dim = c(k, k))
  x
}

svysem_single_nobs <- function(fit) {
  ngroups <- lavaan::lavInspect(fit, "ngroups")
  if (length(ngroups) != 1L || ngroups != 1L) {
    svysem_stop("Multiple groups are not supported for MLM.", "svysem_error_data")
  }
  nobs <- lavaan::lavInspect(fit, "nobs")
  if (is.list(nobs)) {
    nobs <- unlist(nobs, use.names = FALSE)
  }
  nobs <- as.numeric(nobs)
  if (length(nobs) != 1L || !is.finite(nobs) || nobs <= 0) {
    svysem_stop(
      "MLM could not read a single positive sample size from the lavaan fit.",
      "svysem_error_data"
    )
  }
  nobs
}

svysem_gamma_matrix <- function(fit) {
  gamma <- lavaan::lavInspect(fit, "gamma")
  if (is.list(gamma)) {
    if (length(gamma) != 1L) {
      svysem_stop("Multiple groups are not supported for MLM.", "svysem_error_data")
    }
    gamma <- gamma[[1L]]
  }
  gamma
}

# NACOV is nobs times the design variance of the stacked moments lavaan names
# in gamma. Means and products are estimated jointly; a block-diagonal Gamma
# drops the intercept's covariance with the second moments.
svysem_nacov <- function(prelim, design) {
  ov <- lavaan::lavNames(prelim, "ov")
  gamma <- svysem_gamma_matrix(prelim)
  gnames <- colnames(gamma)
  canonical <- svysem_canonical_names(ov)
  if (is.null(gnames) ||
        !(identical(gnames, canonical$full) || identical(gnames, canonical$cov))) {
    unexpected <- if (is.null(gnames)) {
      character(0)
    } else {
      setdiff(gnames, canonical$full)
    }
    detail <- if (length(unexpected)) {
      paste0(" Unrecognized moment(s): ", paste(unexpected, collapse = ", "), ".")
    } else {
      ""
    }
    svysem_stop(
      paste0(
        "MLM does not support this sample-moment layout.", detail,
        " fixed.x = TRUE is supported when gamma holds the observed second ",
        "moments, and the means when a mean structure is estimated. ",
        "conditional.x = TRUE is not supported."
      ),
      "svysem_error_data"
    )
  }

  moments <- svysem_moment_matrix(design, ov)[, gnames, drop = FALSE]
  safe <- sprintf(".svysem_m%04d", seq_len(ncol(moments)))
  extra <- as.data.frame(moments, stringsAsFactors = FALSE)
  names(extra) <- safe
  design_moments <- design
  variables <- data.frame(design$variables, check.names = FALSE, stringsAsFactors = FALSE)
  variables[, safe] <- extra
  design_moments$variables <- variables

  estimated <- survey::svymean(stats::reformulate(safe), design_moments)
  variance <- svysem_plain_matrix(attr(estimated, "var"))
  if (nrow(variance) != length(gnames) || any(!is.finite(variance))) {
    svysem_stop(
      "The design variance of the sample moments is not a finite Gamma.",
      "svysem_error_data"
    )
  }
  svysem_single_nobs(prelim) * variance
}

svysem_lavaan_args <- function(model, data, weight, dots, estimator = "MLM",
                               nacov = NULL, wls_v = NULL, extra = NULL) {
  blocked <- c(
    "data", "sampling.weights", "estimator", "NACOV", "ov.order",
    "sample.cov", "sample.mean", "sample.th", "sample.nobs", "WLS.V"
  )
  forbidden <- intersect(
    names(dots),
    c(
      "data", "sampling.weights", "sample.cov", "sample.mean", "sample.nobs",
      "NACOV", "WLS.V"
    )
  )
  if (length(forbidden)) {
    svysem_stop(
      paste0(
        "Arguments ", paste(forbidden, collapse = ", "),
        " are taken from the survey design and cannot be supplied to svysem()."
      ),
      "svysem_error_data"
    )
  }
  user <- dots[setdiff(names(dots), blocked)]
  # NACOV makes lavaan switch ov.order to "data" unless force.model is set,
  # which would pair Gamma with the wrong variables. The preliminary fit stays
  # in model order and does not receive NACOV.
  args <- c(
    list(model = model, data = data),
    user,
    list(sampling.weights = weight, estimator = estimator)
  )
  if (is.null(nacov)) {
    args$ov.order <- "model"
  } else {
    args$NACOV <- nacov
    # lavaan turns ov.order "force.model" into "model" when NACOV is set,
    # and otherwise reorders statistics to the data. A supplied NACOV is also
    # the DWLS weight matrix unless WLS.V is passed alongside it.
    args$ov.order <- "force.model"
  }
  if (!is.null(wls_v)) {
    args$WLS.V <- wls_v
  }
  if (!is.null(extra)) {
    args <- args[!names(args) %in% names(extra)]
    args <- c(args, extra)
  }
  svysem_allow_continuous_wlsmv(args, model, data)
}

# lavaan 0.7 stops on WLSMV when ordered= was omitted, including for
# continuous data. ordered = FALSE keeps that continuous path. Ordered
# factors are left unset so lavaan still treats them as categorical.
svysem_allow_continuous_wlsmv <- function(args, model, data) {
  if (!identical(args$estimator, "WLSMV") || !is.null(args$ordered) ||
        !is.data.frame(data)) {
    return(args)
  }
  observed <- intersect(svysem_observed_names(model), names(data))
  categorical <- length(observed) > 0L && any(vapply(
    data[observed], inherits, logical(1), what = "ordered"
  ))
  if (!categorical) {
    args$ordered <- FALSE
  }
  args
}

svysem_analysis_frame <- function(design, ov, zero_incomplete = TRUE) {
  svysem_check_weights(design)
  # FIML keeps incomplete rows at their sampling weight. Listwise sets those
  # weights to 0 and leaves the row in the design.
  if (isTRUE(zero_incomplete)) {
    na_rows <- !stats::complete.cases(design$variables[, ov, drop = FALSE])
    design <- svysem_zero_incomplete(design, na_rows)
  }
  weight <- svysem_weight_name(c(names(design$variables), ov))
  data <- data.frame(design$variables, check.names = FALSE, stringsAsFactors = FALSE)
  data[[weight]] <- as.numeric(weights(design, "sampling"))
  list(design = design, data = data, weight = weight)
}

svysem_fit_mlm <- function(model, design, fitter, dots, replicates = NULL,
                           replicate.type = NULL, parallel = NULL) {
  ov <- svysem_observed_names(model)
  if (!length(ov)) {
    svysem_stop(
      "The model has no observed variables to pass to lavaan.",
      "svysem_error_variables"
    )
  }
  svysem_check_variables(ov, names(design$variables))
  if (!is.null(dots$group)) {
    return(svysem_fit_mlm_grouped(
      model, design, fitter, dots, ov, replicates, replicate.type, parallel
    ))
  }
  model_cols <- design$variables[, ov, drop = FALSE]
  svysem_check_arguments(dots, model_cols)
  non_numeric <- ov[!vapply(model_cols, is.numeric, logical(1))]
  if (length(non_numeric)) {
    svysem_stop(
      paste0(
        "Model variables must be numeric for MLM: ",
        paste(non_numeric, collapse = ", "), "."
      ),
      "svysem_error_data"
    )
  }
  prepared <- svysem_analysis_frame(design, ov)
  preliminary <- do.call(
    fitter,
    svysem_lavaan_args(model, prepared$data, prepared$weight, dots)
  )
  nacov <- svysem_nacov(preliminary, prepared$design)
  final <- do.call(
    fitter,
    svysem_lavaan_args(model, prepared$data, prepared$weight, dots, nacov = nacov)
  )
  list(fit = final, design = design)
}

svysem_fit_design <- function(model, design, fitter, svysem_call, estimator, dots,
                               replicates = NULL, replicate.type = NULL,
                               parallel = NULL) {
  estimator <- svysem_check_estimator(estimator)
  svysem_check_missing(dots$missing, estimator)
  if (!inherits(design, "survey.design") && !inherits(design, "svyrep.design")) {
    svysem_stop_design(design)
  }
  grouped_mlm <- identical(estimator, "MLM") && !is.null(dots$group)
  svysem_check_replicate_request(
    estimator, design, replicates, replicate.type, parallel,
    allow_replicates = grouped_mlm
  )
  # svymean and svytotal on a svyrep.design already return replicate variances.
  # WLSMV builds that replicate design itself when the input does not have one.
  # Multigroup MLM does too: lavaan cannot take the between-group block of NACOV.
  replication <- inherits(design, "svyrep.design") ||
    identical(estimator, "WLSMV") || grouped_mlm
  variance_method <- if (replication) "replication" else "linearization"
  fitted <- if (identical(estimator, "WLSMV")) {
    svysem_fit_wlsmv(
      model, design, fitter, dots, replicates, replicate.type, parallel
    )
  } else if (identical(estimator, "MLR")) {
    svysem_fit_mlr(model, design, fitter, dots)
  } else {
    svysem_fit_mlm(
      model, design, fitter, dots, replicates, replicate.type, parallel
    )
  }
  svysem_new(
    fitted$fit, design, svysem_call,
    estimator = estimator, variance_method = variance_method
  )
}

# Point estimates are weighted maximum likelihood. se = "standard" on the
# full sample gives lavaan a covariance slot to overwrite. Replicates skip it.
svysem_mlm_point <- function(fitter, model, data, weight, dots, standard_errors) {
  args <- svysem_lavaan_args(model, data, weight, dots, estimator = "ML")
  args$se <- if (isTRUE(standard_errors)) "standard" else "none"
  args$test <- "standard"
  do.call(fitter, args)
}

svysem_one_coef_replicate <- function(fitter, model, prepared, dots, weights, target) {
  data_replicate <- prepared$data
  data_replicate[[prepared$weight]] <- as.numeric(weights)
  fit_replicate <- svysem_mlm_point(
    fitter, model, data_replicate, prepared$weight, dots, FALSE
  )
  estimate <- lavaan::coef(fit_replicate)
  converged <- isTRUE(all(lavaan::lavInspect(fit_replicate, "converged")))
  if (!converged || !identical(names(estimate), names(target)) ||
        any(!is.finite(estimate))) {
    return(list(
      ok = FALSE,
      message = "the replicate did not return the full-sample coefficients"
    ))
  }
  list(ok = TRUE, coef = estimate)
}

# Sampling weights make lavaan append a Yuan-Bentler statistic even when
# test = "standard". That statistic uses the weights only, not the primary
# sampling units. Multigroup MLM keeps the maximum-likelihood chi-square.
svysem_ml_test_only <- function(test) {
  if (!is.list(test) || !length(test)) {
    return(test)
  }
  keep <- vapply(test, function(element) {
    is.list(element) && identical(element$test, "standard")
  }, logical(1))
  if (!any(keep)) {
    return(test)
  }
  test[keep]
}

svysem_strip_scaled_test <- function(fit) {
  fit@test <- svysem_ml_test_only(fit@test)
  if (is.list(fit@baseline$test)) {
    fit@baseline$test <- svysem_ml_test_only(fit@baseline$test)
  }
  if (is.character(fit@Options$test)) {
    fit@Options$test <- "standard"
  }
  fit
}

# lavaan stores one NACOV per group, so the covariance between groups is
# dropped. Primary sampling units that contain both groups make a
# between-group contrast miss that term. The coefficients are refit on each
# replicate weight instead. The chi-square stays the maximum-likelihood
# statistic: there is no Satorra-Bentler adjustment on this route.
svysem_replicate_coefficients <- function(fitter, model, prepared, dots,
                                          rep_design, parallel) {
  full <- svysem_mlm_point(
    fitter, model, prepared$data, prepared$weight, dots, TRUE
  )
  if (!isTRUE(all(lavaan::lavInspect(full, "converged")))) {
    svysem_stop(
      "The multigroup MLM model did not converge.",
      "svysem_error_data"
    )
  }
  target <- lavaan::coef(full)
  replicated <- weights(rep_design, "analysis")
  if (!is.matrix(replicated) || nrow(replicated) != nrow(prepared$data)) {
    svysem_stop(
      "Replicate weights could not be aligned with the analysis rows.",
      "svysem_error_design"
    )
  }
  loop <- svysem_resolve_parallel(parallel)
  n_replicates <- ncol(replicated)
  results <- loop(seq_len(n_replicates), function(replicate) {
    tryCatch(
      svysem_one_coef_replicate(
        fitter, model, prepared, dots, replicated[, replicate], target
      ),
      error = function(e) list(ok = FALSE, message = conditionMessage(e))
    )
  })
  if (length(results) != n_replicates) {
    svysem_stop(
      "parallel did not return one result per replicate.",
      "svysem_error_design"
    )
  }
  statistics <- matrix(NA_real_, n_replicates, length(target))
  for (replicate in seq_len(n_replicates)) {
    result <- results[[replicate]]
    failed <- inherits(result, "try-error") || !is.list(result) || !isTRUE(result$ok)
    if (failed) {
      message <- if (is.list(result) && !is.null(result$message)) {
        paste(result$message, collapse = " ")
      } else {
        "the replicate did not return coefficients"
      }
      svysem_stop(
        sprintf("Replicate %d failed: %s", replicate, message),
        "svysem_error_data"
      )
    }
    statistics[replicate, ] <- result$coef[names(target)]
  }
  list(
    fit = svysem_strip_scaled_test(full),
    statistics = statistics,
    coef = target
  )
}

svysem_coef_variance <- function(collected, rep_design) {
  variance <- tryCatch(
    survey::svrVar(
      thetas = collected$statistics,
      scale = rep_design$scale,
      rscales = rep_design$rscales,
      mse = isTRUE(rep_design$mse),
      coef = collected$coef
    ),
    error = function(e) e
  )
  if (inherits(variance, "error")) {
    svysem_stop(
      paste0("Replicate variance failed: ", conditionMessage(variance)),
      "svysem_error_data"
    )
  }
  variance <- as.matrix(variance)
  if (nrow(variance) != length(collected$coef) || any(!is.finite(variance))) {
    svysem_stop(
      "The replicate covariance of the coefficients is not finite.",
      "svysem_error_data"
    )
  }
  dimnames(variance) <- list(names(collected$coef), names(collected$coef))
  free_names <- colnames(lavaan::lavInspect(collected$fit, "vcov"))
  if (is.null(free_names) || !all(free_names %in% colnames(variance))) {
    svysem_stop(
      "Replicate coefficients do not match the free parameters.",
      "svysem_error_data"
    )
  }
  variance[free_names, free_names, drop = FALSE]
}

svysem_fit_mlm_grouped <- function(model, design, fitter, dots, ov, replicates,
                                   replicate.type, parallel) {
  svysem_check_group_column(dots$group, names(design$variables))
  model_cols <- design$variables[, ov, drop = FALSE]
  svysem_check_arguments(dots, model_cols, allow_group = TRUE)
  non_numeric <- ov[!vapply(model_cols, is.numeric, logical(1))]
  if (length(non_numeric)) {
    svysem_stop(
      paste0(
        "Model variables must be numeric for MLM: ",
        paste(non_numeric, collapse = ", "), "."
      ),
      "svysem_error_data"
    )
  }
  dots$se <- NULL
  dots$test <- NULL
  prepared <- svysem_analysis_frame(design, ov)
  rep_design <- svysem_make_replicates(prepared$design, replicates, replicate.type)
  collected <- svysem_replicate_coefficients(
    fitter, model, prepared, dots, rep_design, parallel
  )
  covariance <- svysem_coef_variance(collected, rep_design)
  list(
    fit = svysem_inject_covariance(collected$fit, covariance),
    design = design
  )
}

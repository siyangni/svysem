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

svysem_lavaan_args <- function(model, data, weight, dots, nacov = NULL) {
  blocked <- c(
    "data", "sampling.weights", "estimator", "NACOV", "ov.order",
    "sample.cov", "sample.mean", "sample.th", "sample.nobs"
  )
  forbidden <- intersect(
    names(dots),
    c("data", "sampling.weights", "sample.cov", "sample.mean", "sample.nobs", "NACOV")
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
    list(sampling.weights = weight, estimator = "MLM")
  )
  if (is.null(nacov)) {
    args$ov.order <- "model"
  } else {
    args$NACOV <- nacov
    args$ov.order <- "force.model"
  }
  args
}

svysem_analysis_frame <- function(design, ov) {
  svysem_check_weights(design)
  na_rows <- !stats::complete.cases(design$variables[, ov, drop = FALSE])
  design <- svysem_zero_incomplete(design, na_rows)
  weight <- svysem_weight_name(c(names(design$variables), ov))
  data <- data.frame(design$variables, check.names = FALSE, stringsAsFactors = FALSE)
  data[[weight]] <- as.numeric(weights(design, "sampling"))
  list(design = design, data = data, weight = weight)
}

svysem_fit_mlm <- function(model, design, fitter, dots) {
  ov <- svysem_observed_names(model)
  if (!length(ov)) {
    svysem_stop(
      "The model has no observed variables to pass to lavaan.",
      "svysem_error_variables"
    )
  }
  svysem_check_variables(ov, names(design$variables))
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

svysem_fit_design <- function(model, design, fitter, svysem_call, estimator, dots) {
  svysem_check_estimator(estimator)
  svysem_check_missing(dots$missing)
  if (!inherits(design, "survey.design")) {
    svysem_stop_design(design)
  }
  fitted <- svysem_fit_mlm(model, design, fitter, dots)
  svysem_new(fitted$fit, design, svysem_call)
}

# WLSMV point estimates use lavaan's weighted DWLS matrix. The design enters
# only through the replicate-weight covariance of lavaan's own sample
# statistics (thresholds, continuous means and variances, and correlations).
# survey::svrVar supplies the jackknife or bootstrap scale; it is not
# reimplemented here. lavaan builds that DWLS matrix from the diagonal of
# NACOV, so the final fit is given lavaan's own WLS.V and the design NACOV
# replaces only the variance. Mplus TYPE=COMPLEX uses its own weight matrix,
# and that choice has not been compared with this one.

svysem_replicate_types <- c("auto", "JK1", "JKn", "BRR", "bootstrap")

svysem_nobs_vector <- function(fit) {
  nobs <- lavaan::lavInspect(fit, "nobs")
  if (is.list(nobs)) {
    nobs <- unlist(nobs, use.names = FALSE)
  }
  nobs <- as.numeric(nobs)
  if (!length(nobs) || any(!is.finite(nobs)) || any(nobs <= 0)) {
    svysem_stop(
      "The lavaan fit did not report a positive sample size.",
      "svysem_error_data"
    )
  }
  nobs
}

svysem_resolve_parallel <- function(parallel) {
  if (is.null(parallel)) {
    return(lapply)
  }
  if (!is.function(parallel)) {
    svysem_stop(
      "parallel must be NULL or a function.",
      "svysem_error_design"
    )
  }
  parallel
}

svysem_make_replicates <- function(design, replicates, replicate.type) {
  if (inherits(design, "svyrep.design")) {
    return(design)
  }
  type <- if (is.null(replicate.type)) "auto" else replicate.type
  if (length(type) != 1L || !is.character(type) || anyNA(type) ||
        !type %in% svysem_replicate_types) {
    svysem_stop(
      paste0(
        "replicate.type must be one of ",
        paste(sprintf("\"%s\"", svysem_replicate_types), collapse = ", "),
        "."
      ),
      "svysem_error_design"
    )
  }
  # Jackknife and BRR counts are fixed by the primary sampling units.
  # replicates is only the bootstrap sample size.
  args <- list(design = design, type = type, compress = FALSE)
  if (identical(type, "bootstrap")) {
    if (is.null(replicates)) {
      svysem_stop(
        "bootstrap requires replicates, the number of replicate weights.",
        "svysem_error_design"
      )
    }
    if (length(replicates) != 1L || !is.numeric(replicates) ||
          !is.finite(replicates) || replicates < 2 ||
          abs(replicates - round(replicates)) > 1e-8) {
      svysem_stop(
        "replicates must be a single integer of 2 or more.",
        "svysem_error_design"
      )
    }
    args$replicates <- as.integer(round(replicates))
  } else if (!is.null(replicates)) {
    svysem_stop(
      paste0(
        "replicates is the number of bootstrap replicates. ",
        "Jackknife and BRR replicate counts are fixed by the primary sampling units."
      ),
      "svysem_error_design"
    )
  }
  rep_design <- tryCatch(
    do.call(survey::as.svrepdesign, args),
    error = function(e) e
  )
  if (inherits(rep_design, "error")) {
    svysem_stop(
      paste0(
        "Replicate weights could not be built: ",
        conditionMessage(rep_design)
      ),
      "svysem_error_design"
    )
  }
  rep_design
}

svysem_check_group_column <- function(group, data_names) {
  if (is.null(group)) {
    return(invisible(NULL))
  }
  if (length(group) != 1L || !is.character(group) || anyNA(group) ||
        !group %in% data_names) {
    svysem_stop(
      "group must name one column in the design.",
      "svysem_error_variables"
    )
  }
  invisible(NULL)
}

svysem_lookup_stat <- function(name, samp) {
  if (grepl("|", name, fixed = TRUE)) {
    threshold <- samp$th
    if (is.null(threshold) || is.null(names(threshold)) ||
          !name %in% names(threshold)) {
      return(NA_real_)
    }
    return(unname(threshold[[name]]))
  }
  if (endsWith(name, "~1")) {
    key <- substr(name, 1L, nchar(name) - 2L)
    mean <- samp$mean
    if (is.null(names(mean)) || !key %in% names(mean)) {
      return(NA_real_)
    }
    return(unname(mean[[key]]))
  }
  if (grepl("~~", name, fixed = TRUE)) {
    parts <- strsplit(name, "~~", fixed = TRUE)[[1L]]
    covariance <- samp$cov
    if (length(parts) != 2L || is.null(dimnames(covariance)[[1L]]) ||
          !parts[[1L]] %in% rownames(covariance) ||
          !parts[[2L]] %in% colnames(covariance)) {
      return(NA_real_)
    }
    return(unname(covariance[parts[[1L]], parts[[2L]]]))
  }
  NA_real_
}

svysem_stat_vector <- function(gamma, samp) {
  if (!is.matrix(gamma) || is.null(colnames(gamma))) {
    svysem_stop(
      "WLSMV could not read the names of lavaan's sample statistics.",
      "svysem_error_data"
    )
  }
  statistic_names <- colnames(gamma)
  out <- vapply(
    statistic_names, svysem_lookup_stat, numeric(1), samp = samp,
    USE.NAMES = TRUE
  )
  if (any(!is.finite(out))) {
    bad <- statistic_names[!is.finite(out)]
    svysem_stop(
      paste0(
        "WLSMV does not support this sample-statistic layout. ",
        "Unrecognized statistic(s): ", paste(bad, collapse = ", "), "."
      ),
      "svysem_error_data"
    )
  }
  out
}

svysem_sample_vectors <- function(fit) {
  ngroups <- lavaan::lavInspect(fit, "ngroups")
  if (length(ngroups) != 1L || !is.finite(ngroups) || ngroups < 1) {
    svysem_stop(
      "WLSMV could not read the number of groups.",
      "svysem_error_data"
    )
  }
  ngroups <- as.integer(ngroups)
  gamma <- lavaan::lavInspect(fit, "gamma")
  samp <- lavaan::lavInspect(fit, "sampstat")
  if (is.matrix(gamma)) {
    gamma <- list(gamma)
  }
  if (ngroups == 1L) {
    if (is.null(samp$cov) && is.list(samp) && length(samp) == 1L) {
      samp <- samp[[1L]]
    }
    return(list(svysem_stat_vector(gamma[[1L]], samp)))
  }
  if (!is.list(gamma) || !is.list(samp) || length(gamma) != ngroups ||
        length(samp) != ngroups) {
    svysem_stop(
      "WLSMV sample statistics and Gamma do not have one block per group.",
      "svysem_error_data"
    )
  }
  vectors <- lapply(seq_len(ngroups), function(group) {
    svysem_stat_vector(gamma[[group]], samp[[group]])
  })
  stats::setNames(vectors, names(samp))
}

# do.fit = FALSE still returns sampstat and the Gamma column names. Those
# names are the order of the replicate statistic. se and test are cleared so
# a user setting cannot turn the replicate fit into a different estimator.
svysem_wlsmv_sample <- function(fitter, model, data, weight, dots) {
  args <- svysem_lavaan_args(
    model, data, weight, dots, estimator = "WLSMV"
  )
  args$do.fit <- FALSE
  args$se <- "none"
  args$test <- "none"
  # One warning on the returned fit is enough. Replicate sample statistics
  # repeat lavaan's continuous-data DWLS note once per weight.
  withCallingHandlers(
    do.call(fitter, args),
    warning = function(w) {
      if (grepl("not recommended for continuous data", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

svysem_replicate_statistics <- function(fitter, model, prepared, dots, rep_design,
                                        parallel) {
  full <- svysem_wlsmv_sample(
    fitter, model, prepared$data, prepared$weight, dots
  )
  full_vectors <- svysem_sample_vectors(full)
  full_stat <- unlist(full_vectors, use.names = FALSE)
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
      {
        data_replicate <- prepared$data
        data_replicate[[prepared$weight]] <- as.numeric(replicated[, replicate])
        fit_replicate <- svysem_wlsmv_sample(
          fitter, model, data_replicate, prepared$weight, dots
        )
        list(
          ok = TRUE,
          stat = unlist(svysem_sample_vectors(fit_replicate), use.names = FALSE)
        )
      },
      error = function(e) list(ok = FALSE, message = conditionMessage(e))
    )
  })
  if (length(results) != n_replicates) {
    svysem_stop(
      "parallel did not return one result per replicate.",
      "svysem_error_design"
    )
  }
  statistics <- matrix(NA_real_, n_replicates, length(full_stat))
  for (replicate in seq_len(n_replicates)) {
    result <- results[[replicate]]
    failed <- inherits(result, "try-error") || !is.list(result) || !isTRUE(result$ok)
    if (failed) {
      message <- if (is.list(result) && !is.null(result$message)) {
        paste(result$message, collapse = " ")
      } else {
        "the replicate did not return sample statistics"
      }
      svysem_stop(
        sprintf("Replicate %d failed: %s", replicate, message),
        "svysem_error_data"
      )
    }
    if (length(result$stat) != length(full_stat) || any(!is.finite(result$stat))) {
      svysem_stop(
        sprintf(
          "Replicate %d did not return a finite sample-statistic vector.",
          replicate
        ),
        "svysem_error_data"
      )
    }
    statistics[replicate, ] <- result$stat
  }
  list(
    full = full,
    wls_v = svysem_wls_matrix(full),
    vectors = full_vectors,
    statistics = statistics,
    coef = full_stat
  )
}

# lavaan uses diag(NACOV) as the DWLS weight matrix when NACOV is supplied.
# The matrix from the weighted sample, before the design NACOV is attached,
# keeps the point estimates on lavaan's weighted solution.
svysem_wls_matrix <- function(fit) {
  weight <- lavaan::lavInspect(fit, "WLS.V")
  if (is.null(weight)) {
    svysem_stop(
      "lavaan did not return a DWLS weight matrix.",
      "svysem_error_data"
    )
  }
  weight
}

# lavaan takes one NACOV per group, so the between-group block of the
# replicate covariance is dropped.
svysem_nacov_from_replicates <- function(collected, rep_design) {
  nobs <- svysem_nobs_vector(collected$full)
  lengths <- vapply(collected$vectors, length, integer(1))
  if (length(nobs) != length(lengths)) {
    svysem_stop(
      "Group sample sizes do not match the sample-statistic blocks.",
      "svysem_error_data"
    )
  }
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
      "The replicate covariance of the sample statistics is not finite.",
      "svysem_error_data"
    )
  }
  if (length(lengths) == 1L) {
    return(nobs[[1L]] * svysem_plain_matrix(variance))
  }
  nacov <- vector("list", length(lengths))
  offset <- 0L
  for (group in seq_along(lengths)) {
    index <- offset + seq_len(lengths[[group]])
    block <- variance[index, index, drop = FALSE]
    nacov[[group]] <- nobs[[group]] * svysem_plain_matrix(block)
    offset <- offset + lengths[[group]]
  }
  stats::setNames(nacov, names(collected$vectors))
}

svysem_fit_wlsmv <- function(model, design, fitter, dots, replicates,
                             replicate.type, parallel) {
  ov <- svysem_observed_names(model)
  if (!length(ov)) {
    svysem_stop(
      "The model has no observed variables to pass to lavaan.",
      "svysem_error_variables"
    )
  }
  svysem_check_variables(ov, names(design$variables))
  svysem_check_group_column(dots$group, names(design$variables))
  svysem_check_arguments(
    dots, design$variables[, ov, drop = FALSE], estimator = "WLSMV"
  )
  # se and test are replaced so the returned standard errors are the
  # replicate-weight ones. The DWLS weight matrix stays lavaan's default.
  dots$se <- NULL
  dots$test <- NULL
  dots$do.fit <- NULL
  missing_value <- dots$missing
  listwise <- is.null(missing_value) || (
    length(missing_value) == 1L && !is.na(missing_value) &&
      tolower(missing_value) %in% c("listwise", "default")
  )
  prepared <- svysem_analysis_frame(
    design, ov, zero_incomplete = listwise
  )
  rep_design <- svysem_make_replicates(prepared$design, replicates, replicate.type)
  collected <- svysem_replicate_statistics(
    fitter, model, prepared, dots, rep_design, parallel
  )
  nacov <- svysem_nacov_from_replicates(collected, rep_design)
  final <- do.call(
    fitter,
    svysem_lavaan_args(
      model, prepared$data, prepared$weight, dots,
      estimator = "WLSMV", nacov = nacov, wls_v = collected$wls_v
    )
  )
  list(fit = final, design = design)
}

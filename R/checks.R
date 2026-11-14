# Classed conditions. The specific class is first, then the parent svysem_error.
svysem_stop <- function(message, class) {
  # call is NULL so the condition has no call. stop() on a condition
  # rejects a further call. argument, which would only add a warning.
  stop(errorCondition(
    message,
    class = c(class, "svysem_error"),
    call = NULL
  ))
}

svysem_check_estimator <- function(estimator) {
  supplied <- estimator
  matched <- tryCatch(
    match.arg(supplied, c("MLM", "WLSMV", "MLR")),
    error = function(e) NA_character_
  )
  if (!identical(matched, "MLM")) {
    label <- if (is.character(supplied) && length(supplied) == 1L) {
      supplied
    } else {
      paste(as.character(supplied), collapse = ", ")
    }
    svysem_stop(
      sprintf(
        "estimator \"%s\" is not implemented. Only \"MLM\" is available.",
        label
      ),
      "svysem_error_estimator"
    )
  }
  "MLM"
}

# missing = "ml" (and the FIML aliases lavaan accepts) is reserved for MLR.
# Any other value except listwise / default is a data error.
svysem_check_missing <- function(missing_value) {
  if (is.null(missing_value)) {
    return(invisible(NULL))
  }
  if (length(missing_value) != 1L || !is.character(missing_value) ||
        anyNA(missing_value)) {
    svysem_stop(
      "missing must be a single string. MLM accepts \"listwise\" or the default.",
      "svysem_error_data"
    )
  }
  value <- tolower(missing_value)
  if (value %in% c("ml", "fiml", "ml.x")) {
    svysem_stop(
      paste0(
        "missing = \"", missing_value, "\" is not available for MLM. ",
        "MLM is complete data only. The MLR estimator will handle missing = \"ml\"."
      ),
      "svysem_error_missing"
    )
  }
  if (!value %in% c("listwise", "default")) {
    svysem_stop(
      paste0(
        "missing = \"", missing_value, "\" is not supported for MLM. ",
        "Use listwise deletion, which is the default. Incomplete rows stay in ",
        "the design with weight 0 so the primary sampling units are kept."
      ),
      "svysem_error_data"
    )
  }
  invisible(NULL)
}

svysem_stop_design <- function(design) {
  svysem_stop(
    paste0(
      "design must be a survey.design or survey.design2 object. ",
      "The class supplied was: ", paste(class(design), collapse = ", "), "."
    ),
    "svysem_error_design"
  )
}

svysem_prepare_dispatch <- function(estimator, dots) {
  svysem_check_estimator(estimator)
  svysem_check_missing(dots$missing)
  invisible(NULL)
}

svysem_observed_names <- function(model) {
  lavaan::lavNames(lavaan::lavaanify(model), "ov")
}

svysem_ordered_requested <- function(ordered) {
  if (is.null(ordered)) {
    return(FALSE)
  }
  if (is.logical(ordered) && length(ordered) == 1L && !is.na(ordered) &&
        !isTRUE(ordered)) {
    return(FALSE)
  }
  if (is.character(ordered) && length(ordered) == 0L) {
    return(FALSE)
  }
  TRUE
}

# group=, cluster=, ordered indicators, and multiple groups stop before lavaan
# sees them. cluster= would double-count the primary sampling units.
svysem_check_arguments <- function(dots, variables) {
  if (svysem_ordered_requested(dots$ordered)) {
    svysem_stop(
      "Ordered indicators are not supported for MLM.",
      "svysem_error_data"
    )
  }
  ordered_cols <- character(0)
  if (length(variables) && !is.null(variables)) {
    ordered_cols <- names(variables)[vapply(variables, is.ordered, logical(1))]
  }
  if (length(ordered_cols)) {
    svysem_stop(
      paste0(
        "Ordered indicators are not supported for MLM: ",
        paste(ordered_cols, collapse = ", "), "."
      ),
      "svysem_error_data"
    )
  }
  if (!is.null(dots$group)) {
    svysem_stop(
      paste0(
        "A group= argument is not supported. MLM fits a single group. ",
        "Stratification comes from the design."
      ),
      "svysem_error_data"
    )
  }
  if (!is.null(dots$cluster)) {
    svysem_stop(
      paste0(
        "A cluster= argument is not supported. Clustering comes from the ",
        "design; passing cluster would double-count primary sampling units."
      ),
      "svysem_error_data"
    )
  }
  ngroups <- dots$ngroups
  if (!is.null(ngroups) && is.numeric(ngroups) && length(ngroups) == 1L &&
        !is.na(ngroups) && ngroups > 1) {
    svysem_stop(
      "Multiple groups are not supported for MLM.",
      "svysem_error_data"
    )
  }
  invisible(NULL)
}

svysem_check_variables <- function(ov, data_names) {
  if (is.null(ov)) {
    return(invisible(NULL))
  }
  missing_vars <- ov[!ov %in% data_names]
  if (length(missing_vars)) {
    svysem_stop(
      paste0(
        "Model variable(s) not found in the design: ",
        paste(missing_vars, collapse = ", "), "."
      ),
      "svysem_error_variables"
    )
  }
  invisible(NULL)
}

svysem_check_weights <- function(design) {
  weights <- as.numeric(weights(design, "sampling"))
  if (any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
    svysem_stop(
      paste0(
        "Sampling weights must be non-negative and finite, and at least one ",
        "weight must be positive. Zero weights are kept."
      ),
      "svysem_error_weights"
    )
  }
  weights
}

# Listwise deletion as a domain: incomplete rows stay, with sampling weight 0,
# so stratum and PSU counts are not rebuilt by dropping rows.
svysem_zero_incomplete <- function(design, rows) {
  if (!any(rows)) {
    return(design)
  }
  prob <- design$prob * 1
  if (is.matrix(prob)) {
    prob[rows, ] <- Inf
  } else {
    prob[rows] <- Inf
  }
  design$prob <- prob
  design
}

svysem_weight_name <- function(reserved) {
  name <- ".svysem_weight"
  while (name %in% reserved) {
    name <- paste0(name, "x")
  }
  name
}

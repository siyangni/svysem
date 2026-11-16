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
  if (is.na(matched)) {
    label <- if (is.character(supplied) && length(supplied) == 1L) {
      supplied
    } else {
      paste(as.character(supplied), collapse = ", ")
    }
    svysem_stop(
      sprintf(
        paste0(
          "estimator \"%s\" is not implemented. ",
          "Available estimators are \"MLM\", \"WLSMV\" and \"MLR\"."
        ),
        label
      ),
      "svysem_error_estimator"
    )
  }
  matched
}

# missing = "ml" (and the FIML aliases lavaan accepts) is MLR only.
# MLM messages below stay fixed: they are the public error text.
svysem_check_missing <- function(missing_value, estimator = "MLM") {
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
  if (identical(estimator, "MLM")) {
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
    return(invisible(NULL))
  }
  if (identical(estimator, "MLR")) {
    if (value %in% c("ml", "fiml", "ml.x", "listwise", "default")) {
      return(invisible(NULL))
    }
    svysem_stop(
      paste0(
        "missing = \"", missing_value, "\" is not supported for MLR. ",
        "Use listwise deletion, which is the default, or missing = \"ml\"."
      ),
      "svysem_error_data"
    )
  }
  # WLSMV: listwise (the default) or pairwise. Full information is MLR.
  if (value %in% c("ml", "fiml", "ml.x")) {
    svysem_stop(
      paste0(
        "missing = \"", missing_value, "\" is not available for WLSMV. ",
        "Full-information missing data is handled by estimator = \"MLR\"."
      ),
      "svysem_error_missing"
    )
  }
  if (!value %in% c("listwise", "default", "pairwise")) {
    svysem_stop(
      paste0(
        "missing = \"", missing_value, "\" is not supported for WLSMV. ",
        "Use \"listwise\" or \"pairwise\"."
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
  matched <- svysem_check_estimator(estimator)
  svysem_check_missing(dots$missing, matched)
  invisible(NULL)
}

# replicates and replicate.type build a svyrep.design for WLSMV and for
# multigroup MLM. An object that is already a svyrep.design keeps its own
# replicate weights. parallel on single-group MLM or on MLR would otherwise
# fall through to lavaan.
svysem_check_replicate_request <- function(estimator, design, replicates,
                                           replicate.type, parallel,
                                           allow_replicates = FALSE) {
  supplied <- !is.null(replicates) || !is.null(replicate.type)
  if (inherits(design, "svyrep.design") && supplied) {
    svysem_stop(
      paste0(
        "replicates and replicate.type apply when svysem builds replicate ",
        "weights. This design already has them."
      ),
      "svysem_error_design"
    )
  }
  allowed <- identical(estimator, "WLSMV") || isTRUE(allow_replicates)
  if (!allowed && (supplied || !is.null(parallel))) {
    svysem_stop(
      paste0(
        "replicates, replicate.type, and parallel are used by WLSMV, ",
        "and by MLM when group= is set. MLR and single-group MLM take the ",
        "variance from the design as it stands."
      ),
      "svysem_error_design"
    )
  }
  invisible(NULL)
}

svysem_observed_names <- function(model) {
  # lavaanify() rejects c(label1, label2)*x, which is how a multigroup
  # model and measEq.syntax name group-specific parameters. The modifier
  # is not an observed variable. The model passed to lavaan is unchanged.
  syntax <- model
  if (is.character(syntax) && length(syntax) == 1L) {
    syntax <- gsub("c\\([^()]*\\)\\*", "", syntax, perl = TRUE)
  }
  lavaan::lavNames(lavaan::lavaanify(syntax), "ov")
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
# sees them, except WLSMV, which uses those sample statistics.
# cluster= would double-count the primary sampling units for every estimator.
# The MLM strings are the public error text and stay as written.
svysem_check_arguments <- function(dots, variables, estimator = "MLM",
                                  allow_group = FALSE) {
  if (!identical(estimator, "WLSMV")) {
    if (svysem_ordered_requested(dots$ordered)) {
      svysem_stop(
        if (identical(estimator, "MLM")) {
          "Ordered indicators are not supported for MLM."
        } else {
          "Ordered indicators are not supported for MLR."
        },
        "svysem_error_data"
      )
    }
    ordered_cols <- character(0)
    if (length(variables) && !is.null(variables)) {
      ordered_cols <- names(variables)[vapply(variables, is.ordered, logical(1))]
    }
    if (length(ordered_cols)) {
      svysem_stop(
        if (identical(estimator, "MLM")) {
          paste0(
            "Ordered indicators are not supported for MLM: ",
            paste(ordered_cols, collapse = ", "), "."
          )
        } else {
          paste0(
            "Ordered indicators are not supported for MLR: ",
            paste(ordered_cols, collapse = ", "), "."
          )
        },
        "svysem_error_data"
      )
    }
    cx <- dots$conditional.x
    if (!is.null(cx) && !isFALSE(cx)) {
      svysem_stop(
        if (identical(estimator, "MLM")) {
          paste0(
            "MLM does not support this sample-moment layout. ",
            "Unrecognized moment(s): y~x1. ",
            "fixed.x = TRUE is supported when gamma holds the observed second ",
            "moments, and the means when a mean structure is estimated. ",
            "conditional.x = TRUE is not supported."
          )
        } else {
          "MLR does not support conditional.x = TRUE."
        },
        "svysem_error_data"
      )
    }
    if (!is.null(dots$group) && !isTRUE(allow_group)) {
      svysem_stop(
        if (identical(estimator, "MLM")) {
          paste0(
            "A group= argument is not supported. MLM fits a single group. ",
            "Stratification comes from the design."
          )
        } else {
          paste0(
            "A group= argument is not supported. MLR fits a single group. ",
            "Stratification comes from the design."
          )
        },
        "svysem_error_data"
      )
    }
    ngroups <- dots$ngroups
    if (!is.null(ngroups) && is.numeric(ngroups) && length(ngroups) == 1L &&
          !is.na(ngroups) && ngroups > 1) {
      svysem_stop(
        if (identical(estimator, "MLM")) {
          "ngroups= is not supported. Name the groups with group=."
        } else {
          "Multiple groups are not supported for MLR."
        },
        "svysem_error_data"
      )
    }
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
# A svyrep.design stores sampling weights in pweights and has no prob slot.
svysem_zero_incomplete <- function(design, rows) {
  if (!any(rows)) {
    return(design)
  }
  if (inherits(design, "svyrep.design")) {
    design$pweights[rows] <- 0
    replicated <- as.matrix(design$repweights)
    replicated[rows, ] <- 0
    # as.matrix expands a compressed replicate matrix. Keep it expanded and
    # leave combined.weights as it was: analysis weights follow that flag.
    design$repweights <- structure(replicated, class = "repweights")
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

# Turn lavaan constraint syntax into per-imputation values and covariances.
# Labels and := names are rewritten onto free parameters before numericDeriv,
# which is the same derivative mitml::testConstraints uses. Not exported.

svysem_expr_symbols <- function(expr) {
  if (is.symbol(expr)) {
    return(as.character(expr))
  }
  if (!is.call(expr)) {
    return(character())
  }
  parts <- lapply(as.list(expr)[-1], svysem_expr_symbols)
  if (!length(parts)) {
    return(character())
  }
  unique(unlist(parts, use.names = FALSE))
}

svysem_has_inequality <- function(expr) {
  if (!is.call(expr)) {
    return(FALSE)
  }
  head <- as.character(expr[[1]])
  if (head %in% c(">", "<", ">=", "<=", "!=")) {
    return(TRUE)
  }
  any(vapply(as.list(expr)[-1], svysem_has_inequality, logical(1)))
}

svysem_rewrite_equal <- function(expr) {
  if (!is.call(expr)) {
    return(expr)
  }
  args <- lapply(as.list(expr)[-1], svysem_rewrite_equal)
  head <- expr[[1]]
  if (identical(head, as.symbol("=="))) {
    head <- as.symbol("-")
  }
  as.call(c(list(head), args))
}

# map values are either coef names (aliases) or := right-hand sides (text).
svysem_replace_symbols <- function(expr, map, parse_text = FALSE) {
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    if (!nm %in% names(map)) {
      return(expr)
    }
    if (parse_text) {
      return(parse(text = map[[nm]], keep.source = FALSE)[[1]])
    }
    return(as.symbol(map[[nm]]))
  }
  if (!is.call(expr)) {
    return(expr)
  }
  args <- lapply(
    as.list(expr)[-1], svysem_replace_symbols,
    map = map, parse_text = parse_text
  )
  as.call(c(list(expr[[1]]), args))
}

svysem_substitute_definitions <- function(expr, definitions) {
  if (!length(definitions)) {
    return(expr)
  }
  seen <- character()
  repeat {
    used <- intersect(svysem_expr_symbols(expr), names(definitions))
    if (!length(used)) {
      break
    }
    if (any(used %in% seen) || length(seen) > length(definitions)) {
      svysem_stop(
        "Defined parameters refer to each other in a cycle.",
        "svysem_error_constraint"
      )
    }
    seen <- c(seen, used)
    expr <- svysem_replace_symbols(
      expr, definitions[used], parse_text = TRUE
    )
  }
  expr
}

# Free-parameter symbol for each label and each lhs-op-rhs name.
svysem_lavaan_maps <- function(fit, coef_names) {
  pe <- lavaan::parameterEstimates(fit)
  needed <- c("lhs", "op", "rhs", "label")
  if (!all(needed %in% names(pe))) {
    svysem_stop(
      "parameterEstimates() is missing lhs, op, rhs, or label.",
      "svysem_error_constraint"
    )
  }
  aliases <- character()
  definitions <- character()
  for (i in seq_len(nrow(pe))) {
    op <- pe$op[[i]]
    lhs <- pe$lhs[[i]]
    rhs <- pe$rhs[[i]]
    label <- pe$label[[i]]
    if (!is.character(label) || !nzchar(label)) {
      label <- ""
    }
    if (identical(op, ":=")) {
      definitions[[lhs]] <- rhs
      next
    }
    key <- paste0(lhs, op, rhs)
    if (key %in% coef_names) {
      if (nzchar(label) && !identical(label, key)) {
        aliases[[label]] <- key
      }
      next
    }
    if (nzchar(label) && label %in% coef_names) {
      aliases[[key]] <- label
    }
  }
  list(aliases = aliases, definitions = definitions)
}

svysem_prepare_constraint <- function(text, fit, coef_names) {
  if (!is.character(text) || length(text) != 1L || !nzchar(text) || is.na(text)) {
    svysem_stop(
      "Each constraint must be one non-empty string.",
      "svysem_error_constraint"
    )
  }
  expr <- tryCatch(
    parse(text = text, keep.source = FALSE)[[1]],
    error = function(e) {
      svysem_stop(
        sprintf("Could not parse constraint %s.", dQuote(text)),
        "svysem_error_constraint"
      )
    }
  )
  if (svysem_has_inequality(expr)) {
    svysem_stop(
      sprintf("Inequality constraints are not supported (%s).", text),
      "svysem_error_constraint"
    )
  }
  aliases <- character()
  definitions <- character()
  if (inherits(fit, "lavaan")) {
    maps <- svysem_lavaan_maps(fit, coef_names)
    aliases <- maps$aliases
    definitions <- maps$definitions
  }
  expr <- svysem_rewrite_equal(expr)
  expr <- svysem_substitute_definitions(expr, definitions)
  if (length(aliases)) {
    expr <- svysem_replace_symbols(expr, aliases, parse_text = FALSE)
  }
  unknown <- setdiff(svysem_expr_symbols(expr), coef_names)
  if (length(unknown)) {
    svysem_stop(
      sprintf(
        "Unknown name in constraint: %s.",
        paste(unknown, collapse = ", ")
      ),
      "svysem_error_constraint"
    )
  }
  expr
}

# Delta method inside each imputation. Qhat/Uhat are the free-parameter moments.
svysem_constraint_moments <- function(fits, constraints, qu) {
  if (!is.character(constraints) || !length(constraints) || anyNA(constraints)) {
    svysem_stop(
      "constraints must be a character vector of expressions.",
      "svysem_error_constraint"
    )
  }
  nms <- rownames(qu$Qhat)
  m <- ncol(qu$Qhat)
  k <- length(constraints)
  Q <- matrix(NA_real_, k, m)
  U <- array(NA_real_, c(k, k, m))
  for (i in seq_len(m)) {
    theta <- qu$Qhat[, i]
    names(theta) <- nms
    exprs <- lapply(constraints, svysem_prepare_constraint, fit = fits[[i]], coef_names = nms)
    env <- list2env(as.list(theta), parent = baseenv())
    jac <- matrix(NA_real_, k, length(nms))
    vals <- numeric(k)
    for (j in seq_len(k)) {
      for (nm in nms) {
        env[[nm]] <- unname(theta[[nm]])
      }
      nd <- stats::numericDeriv(exprs[[j]], nms, env)
      value <- as.numeric(nd)
      if (length(value) != 1L || !is.finite(value)) {
        svysem_stop(
          "Each constraint must evaluate to one finite number.",
          "svysem_error_constraint"
        )
      }
      gradient <- as.numeric(attr(nd, "gradient"))
      if (length(gradient) != length(nms) || !all(is.finite(gradient))) {
        svysem_stop(
          "Could not differentiate a constraint with respect to the free parameters.",
          "svysem_error_constraint"
        )
      }
      vals[[j]] <- value
      jac[j, ] <- gradient
    }
    Q[, i] <- vals
    U[, , i] <- jac %*% qu$Uhat[, , i] %*% t(jac)
  }
  list(Qhat = Q, Uhat = U)
}

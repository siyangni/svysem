# MLR variance is a design-based sandwich. The bread is lavaan's observed
# information. The meat is the sampling variance of the total score, divided
# by nobs squared, which is the scale on which that information is stored.
# lavaan 0.6-21 has no public slot for a user-supplied meat, so the covariance,
# the parameter-table standard errors, and the Yuan-Bentler test are written
# into the fitted object. Complete-data scores come from lavScores(). Weighted
# full-information scores are the numerical derivative of the casewise
# log-likelihood, because lavScores() rejects that combination.

svysem_is_fiml <- function(missing_value) {
  if (is.null(missing_value) || length(missing_value) != 1L ||
        !is.character(missing_value) || anyNA(missing_value)) {
    return(FALSE)
  }
  tolower(missing_value) %in% c("ml", "fiml", "ml.x")
}

svysem_information <- function(fit) {
  info <- lavaan::lavInspect(fit, "information.observed")
  if (is.list(info) && !is.matrix(info)) {
    info <- info[[1L]]
  }
  info <- as.matrix(info)
  # stats::coef is the S3 generic in this namespace and has no lavaan method.
  parameter_names <- names(lavaan::coef(fit))
  # Repeated names are equality constraints (a shared label). Indexing by
  # those names keeps the first column twice and makes the information singular.
  if (!is.null(colnames(info)) && length(parameter_names) &&
        !anyDuplicated(parameter_names) && !anyDuplicated(colnames(info)) &&
        all(parameter_names %in% colnames(info))) {
    info <- info[parameter_names, parameter_names, drop = FALSE]
  }
  if (nrow(info) != length(parameter_names) || any(!is.finite(info))) {
    svysem_stop(
      "The observed information matrix does not match the free parameters.",
      "svysem_error_data"
    )
  }
  dimnames(info) <- list(parameter_names, parameter_names)
  info
}

svysem_align_rows <- function(fit, scores, n_rows) {
  if (nrow(scores) == n_rows) {
    return(scores)
  }
  index <- fit@Data@case.idx
  if (is.list(index)) {
    index <- index[[1L]]
  }
  index <- as.integer(index)
  if (length(index) != nrow(scores) || any(index < 1L | index > n_rows)) {
    svysem_stop(
      "Casewise scores could not be aligned with the design rows.",
      "svysem_error_data"
    )
  }
  aligned <- matrix(
    0, n_rows, ncol(scores),
    dimnames = list(NULL, colnames(scores))
  )
  aligned[index, ] <- scores
  aligned
}

svysem_lavaan_weight_vector <- function(fit, n_rows) {
  weights_lavaan <- as.numeric(fit@Data@weights[[1L]])
  if (length(weights_lavaan) == n_rows) {
    return(weights_lavaan)
  }
  index <- fit@Data@case.idx
  if (is.list(index)) {
    index <- index[[1L]]
  }
  index <- as.integer(index)
  if (length(weights_lavaan) != length(index)) {
    svysem_stop(
      "lavaan sampling weights could not be aligned with the design.",
      "svysem_error_data"
    )
  }
  full <- numeric(n_rows)
  full[index] <- weights_lavaan
  full
}

# loglik.casewise reads @implied, not the GLIST inside @Model, so both are
# updated before the derivative is taken. The step is absolute, not a relative
# fraction of a parameter that may be near zero.
svysem_loglik_at <- function(fit, theta) {
  fit@Model <- lavaan::lav_model_set_parameters(fit@Model, theta)
  fit@implied <- lavaan::lav_model_implied(fit@Model)
  lavaan::lavInspect(fit, "loglik.casewise")
}

svysem_restore_theta <- function(fit, theta) {
  fit@Model <- lavaan::lav_model_set_parameters(fit@Model, theta)
  fit@implied <- lavaan::lav_model_implied(fit@Model)
  invisible(fit)
}

svysem_numerical_scores <- function(fit) {
  theta <- lavaan::lav_model_get_parameters(fit@Model)
  # The derivative moves @Model. Put the MLE back before the information
  # matrix is read from the same object.
  on.exit(svysem_restore_theta(fit, theta), add = TRUE)
  base <- lavaan::lavInspect(fit, "loglik.casewise")
  if (any(!is.finite(base))) {
    svysem_stop(
      paste0(
        "FIML casewise log-likelihoods are not finite. ",
        "The fit will not fall back to complete cases."
      ),
      "svysem_error_data"
    )
  }
  scores <- matrix(NA_real_, length(base), length(theta))
  for (column in seq_along(theta)) {
    # 1e-6 keeps a central difference within 1e-6 of lavScores on complete
    # data. A larger step was about 1e-5 off that analytic score.
    step <- 1e-6 * max(1, abs(theta[[column]]))
    above <- theta
    below <- theta
    above[[column]] <- theta[[column]] + step
    below[[column]] <- theta[[column]] - step
    scores[, column] <- (
      svysem_loglik_at(fit, above) - svysem_loglik_at(fit, below)
    ) / (2 * step)
  }
  if (any(!is.finite(scores))) {
    svysem_stop(
      paste0(
        "FIML scores are not finite. ",
        "The fit will not fall back to complete cases."
      ),
      "svysem_error_data"
    )
  }
  scores
}

svysem_check_score_balance <- function(weighted) {
  total <- colSums(weighted)
  scale <- pmax(sqrt(colSums(weighted^2)), 1)
  if (any(!is.finite(total)) || any(abs(total) > 1e-3 * scale)) {
    svysem_stop(
      paste0(
        "FIML scores are unstable, so the design variance was not computed. ",
        "The fit does not fall back to complete cases."
      ),
      "svysem_error_data"
    )
  }
  invisible(NULL)
}

# lavaan 0.6 returns unweighted scores. lavaan 0.7 multiplies them by the
# sampling weight, so the column sums equal the weighted gradient. After
# listwise deletion lavScores() errors ("non-conformable arrays") even though
# the casewise log-likelihood is finite. That derivative is already multiplied
# by lavaan's sampling weight.
svysem_lavscores_weighted <- function() {
  utils::packageVersion("lavaan") >= "0.7"
}

svysem_casewise_scores <- function(fit, numerical) {
  if (isTRUE(numerical)) {
    return(list(scores = svysem_numerical_scores(fit), weighted = TRUE))
  }
  # remove.duplicated = FALSE keeps one column per coef() entry. lavaan 0.7
  # otherwise drops the repeated equality label and the meat no longer
  # matches the observed information. The dotted argument is the 0.6 name;
  # 0.7 accepts it as an alias of remove_duplicated.
  raw <- tryCatch(
    lavaan::lavScores(
      fit,
      remove.duplicated = FALSE,
      remove.empty.cases = FALSE
    ),
    error = function(e) e
  )
  usable <- !inherits(raw, "error") && is.matrix(raw) &&
    nrow(raw) > 0L && ncol(raw) > 0L
  if (usable) {
    return(list(scores = raw, weighted = svysem_lavscores_weighted()))
  }
  if (!inherits(raw, "error")) {
    svysem_stop(
      "MLR could not compute casewise scores: lavScores returned no scores",
      "svysem_error_data"
    )
  }
  list(scores = svysem_numerical_scores(fit), weighted = TRUE)
}

# Returns the per-row contribution divided by the original sampling weight.
# svytotal of that matrix estimates the total of lavaan's weighted scores.
# Rows lavaan dropped, and rows with a zero sampling weight, contribute 0.
svysem_adjusted_scores <- function(fit, design, numerical) {
  n_rows <- nrow(design$variables)
  parameter_names <- names(lavaan::coef(fit))
  scored <- svysem_casewise_scores(fit, numerical)
  raw <- scored$scores
  if (ncol(raw) != length(parameter_names)) {
    svysem_stop(
      "Casewise scores do not have one column per free parameter.",
      "svysem_error_data"
    )
  }
  if (!is.null(colnames(raw)) && !anyDuplicated(parameter_names) &&
        !anyDuplicated(colnames(raw)) &&
        all(parameter_names %in% colnames(raw))) {
    raw <- raw[, parameter_names, drop = FALSE]
  }
  colnames(raw) <- parameter_names
  raw[!is.finite(raw)] <- 0
  aligned <- svysem_align_rows(fit, raw, n_rows)
  weights_lavaan <- svysem_lavaan_weight_vector(fit, n_rows)
  weights_design <- as.numeric(weights(design, "sampling"))
  weighted <- if (isTRUE(scored$weighted)) {
    aligned
  } else {
    aligned * weights_lavaan
  }
  if (numerical) {
    svysem_check_score_balance(weighted)
  }
  adjusted <- weighted
  positive <- is.finite(weights_design) & weights_design > 0
  adjusted[positive, ] <- weighted[positive, , drop = FALSE] / weights_design[positive]
  adjusted[!positive, ] <- 0
  adjusted[!is.finite(adjusted)] <- 0
  adjusted
}

svysem_total_meat <- function(design, adjusted, nobs) {
  safe <- sprintf(".svysem_s%04d", seq_len(ncol(adjusted)))
  design_scores <- design
  variables <- data.frame(
    design$variables, check.names = FALSE, stringsAsFactors = FALSE
  )
  variables[, safe] <- adjusted
  design_scores$variables <- variables
  estimated <- tryCatch(
    survey::svytotal(stats::reformulate(safe), design_scores),
    error = function(e) e
  )
  if (inherits(estimated, "error")) {
    svysem_stop(
      paste0(
        "The design variance of the scores failed: ",
        conditionMessage(estimated)
      ),
      "svysem_error_data"
    )
  }
  variance <- as.matrix(stats::vcov(estimated))
  if (nrow(variance) != ncol(adjusted) || any(!is.finite(variance))) {
    svysem_stop(
      "The design variance of the scores is not finite.",
      "svysem_error_data"
    )
  }
  variance / nobs^2
}

svysem_sandwich_pieces <- function(fit, design, numerical) {
  nobs <- svysem_nobs_vector(fit)
  if (length(nobs) != 1L) {
    svysem_stop("MLR fits a single group.", "svysem_error_data")
  }
  adjusted <- svysem_adjusted_scores(fit, design, numerical)
  info <- svysem_information(fit)
  if (ncol(adjusted) != ncol(info)) {
    svysem_stop(
      "Casewise scores do not match the observed information matrix.",
      "svysem_error_data"
    )
  }
  meat <- svysem_total_meat(design, adjusted, nobs[[1L]])
  dimnames(meat) <- dimnames(info)
  list(H = info, meat = meat, nobs = nobs[[1L]])
}

svysem_defined_se <- function(fit, covariance) {
  theta <- lavaan::lav_model_get_parameters(fit@Model, type = "free")
  jacobian <- tryCatch(
    lavaan::lav_func_jacobian_complex(fit@Model@def.function, theta),
    error = function(e) e
  )
  if (inherits(jacobian, "error")) {
    jacobian <- tryCatch(
      lavaan::lav_func_jacobian_simple(fit@Model@def.function, theta),
      error = function(e) e
    )
  }
  if (inherits(jacobian, "error")) {
    svysem_stop(
      paste0(
        "Standard errors for defined parameters failed: ",
        conditionMessage(jacobian)
      ),
      "svysem_error_data"
    )
  }
  if (isTRUE(fit@Model@ceq.simple.only)) {
    jacobian <- jacobian %*% t(fit@Model@ceq.simple.K)
  }
  defined <- diag(jacobian %*% covariance %*% t(jacobian))
  defined[defined < 0] <- NA_real_
  sqrt(defined)
}

# parameterEstimates() reads ParTable$se. vcov() and lavTestWald() read @vcov.
# Defined parameters are not rows of @vcov; their standard errors are the
# delta method applied to the free-parameter covariance.
svysem_inject_covariance <- function(fit, covariance) {
  free_names <- colnames(covariance)
  covariance <- svysem_plain_matrix(covariance)
  if (is.null(free_names) || length(free_names) != ncol(covariance) ||
        any(!is.finite(covariance))) {
    svysem_stop(
      "The coefficient covariance matrix is not finite.",
      "svysem_error_data"
    )
  }
  dimnames(covariance) <- list(free_names, free_names)
  partable <- fit@ParTable
  standard_error <- rep(NA_real_, length(partable$lhs))
  free_variance <- diag(covariance)
  free_variance[free_variance < 0] <- NA_real_
  free_se <- sqrt(free_variance)
  free <- partable$free
  standard_error[free == 0L] <- 0
  positive <- free > 0L
  if (any(positive)) {
    standard_error[positive] <- free_se[free[positive]]
  }
  defined <- which(partable$op == ":=")
  if (length(defined)) {
    standard_error[defined] <- svysem_defined_se(fit, covariance)
  }
  fit@vcov$vcov <- covariance
  fit@ParTable$se <- standard_error
  if (length(fit@Fit@se) == length(standard_error)) {
    fit@Fit@se <- standard_error
  }
  if (length(fit@vcov$se) == length(standard_error)) {
    fit@vcov$se <- standard_error
  } else if (length(fit@vcov$se) == ncol(covariance)) {
    fit@vcov$se <- free_se
  }
  fit
}

# lavaan stores equality constraints in two places. Linear constraints,
# including a shared label, have a basis in eq.constraints.K. Nonlinear
# constraints have only con.jac; its null space is the same kind of basis.
# The sandwich is solved on that reduced space and expanded back to coef().
svysem_null_basis <- function(jacobian, n_parameters) {
  if (!is.matrix(jacobian)) {
    return(NULL)
  }
  if (ncol(jacobian) != n_parameters && nrow(jacobian) == n_parameters) {
    jacobian <- t(jacobian)
  }
  if (ncol(jacobian) != n_parameters || nrow(jacobian) < 1L ||
        any(!is.finite(jacobian))) {
    return(NULL)
  }
  decomposed <- svd(jacobian, nv = n_parameters)
  if (!length(decomposed$d) || !is.finite(decomposed$d[[1L]]) ||
        decomposed$d[[1L]] < 1e-8) {
    return(NULL)
  }
  rank <- sum(decomposed$d > 1e-8 * decomposed$d[[1L]])
  if (rank < 1L || rank >= n_parameters) {
    return(NULL)
  }
  decomposed$v[, seq.int(rank + 1L, n_parameters), drop = FALSE]
}

svysem_constraint_basis <- function(fit, n_parameters) {
  stored <- fit@Model@eq.constraints.K
  if (is.matrix(stored) && nrow(stored) == n_parameters &&
        ncol(stored) > 0L && ncol(stored) < n_parameters &&
        all(is.finite(stored))) {
    return(stored)
  }
  svysem_null_basis(fit@Model@con.jac, n_parameters)
}

svysem_project_pieces <- function(fit, pieces) {
  basis <- svysem_constraint_basis(fit, ncol(pieces$H))
  if (is.null(basis)) {
    return(list(
      H = pieces$H, meat = pieces$meat, nobs = pieces$nobs, expand = NULL
    ))
  }
  list(
    H = t(basis) %*% pieces$H %*% basis,
    meat = t(basis) %*% pieces$meat %*% basis,
    nobs = pieces$nobs,
    expand = basis
  )
}

svysem_solve_bread <- function(information) {
  solved <- tryCatch(solve(information), error = function(e) e)
  if (inherits(solved, "error") || any(!is.finite(solved))) {
    svysem_stop(
      paste0(
        "The observed information matrix is singular, so MLR standard errors ",
        "are not available."
      ),
      "svysem_error_data"
    )
  }
  solved
}

svysem_inject_sandwich <- function(fit, pieces) {
  projected <- svysem_project_pieces(fit, pieces)
  solved <- svysem_solve_bread(projected$H)
  covariance <- solved %*% projected$meat %*% solved
  if (!is.null(projected$expand)) {
    covariance <- projected$expand %*% covariance %*% t(projected$expand)
  }
  dimnames(covariance) <- dimnames(pieces$H)
  if (any(!is.finite(covariance))) {
    svysem_stop(
      "The MLR covariance matrix is not finite.",
      "svysem_error_data"
    )
  }
  svysem_inject_covariance(fit, covariance)
}

svysem_trace_product <- function(fit, pieces) {
  # tr(solve(H) %*% (n * meat)) on the constraint-reduced matrices. The
  # product is an elementwise product with the transpose.
  projected <- svysem_project_pieces(fit, pieces)
  solved <- svysem_solve_bread(projected$H)
  first_order <- projected$nobs * projected$meat
  sum(solved * t(first_order))
}

svysem_refit_partable <- function(partable, data, weight, missing) {
  args <- list(
    model = partable,
    data = data,
    sampling.weights = weight,
    estimator = "ML",
    information = "observed",
    se = "none",
    test = "standard"
  )
  if (!is.null(missing)) {
    args$missing <- missing
  }
  fit <- tryCatch(do.call(lavaan::lavaan, args), error = function(e) e)
  if (inherits(fit, "error")) {
    svysem_stop(
      paste0(
        "The Yuan-Bentler reference model could not be fit: ",
        conditionMessage(fit)
      ),
      "svysem_error_data"
    )
  }
  fit
}

svysem_set_yuan_bentler <- function(test, df, stat, trace_h1, trace_h0, ndat, npar) {
  if (length(df) != 1L || !is.finite(df) || df <= 0) {
    return(test)
  }
  ugamma <- trace_h1 - trace_h0
  scaling <- ugamma / df
  if (!is.finite(scaling) || scaling <= 0 || !is.finite(stat) ||
        !is.finite(trace_h1) || !is.finite(trace_h0) ||
        !is.finite(ndat) || ndat <= 0 || !is.finite(npar) || npar <= 0) {
    svysem_stop(
      paste0(
        "The Yuan-Bentler scaling factor is not positive. ",
        "The independence-model scaling was not kept."
      ),
      "svysem_error_data"
    )
  }
  scaled <- stat / scaling
  index <- which(vapply(test, function(element) {
    is.list(element) && identical(element$test, "yuan.bentler.mplus")
  }, logical(1)))
  if (!length(index)) {
    test <- c(test, list(yuan.bentler.mplus = list(
      test = "yuan.bentler.mplus",
      refdistr = "chisq",
      label = "Yuan-Bentler correction (Mplus variant)"
    )))
    index <- length(test)
  }
  element <- test[[index[[1L]]]]
  element$stat <- scaled
  element$stat.group <- scaled
  element$df <- as.integer(df)
  element$pvalue <- stats::pchisq(scaled, df, lower.tail = FALSE)
  element$scaling.factor <- scaling
  element$scaling.factor.h1 <- trace_h1 / ndat
  element$scaling.factor.h0 <- trace_h0 / npar
  element$trace.UGamma <- ugamma
  element$scaled.test.stat <- stat
  element$scaled.test <- "standard"
  attr(element$stat, "h1") <- trace_h1
  attr(element$stat, "h0") <- trace_h0
  attr(element$stat, "h1.ndat") <- ndat
  test[[index[[1L]]]] <- element
  test
}

svysem_inject_yuan_bentler <- function(fit, design, data, weight, numerical, pieces) {
  if (!is.list(fit@test) || !length(fit@test) || is.null(fit@test[[1L]]$stat)) {
    svysem_stop(
      "MLR did not return a chi-square statistic to scale.",
      "svysem_error_data"
    )
  }
  df <- fit@test[[1L]]$df
  stat <- fit@test[[1L]]$stat
  if (length(df) != 1L || !is.finite(df) || df <= 0) {
    return(fit)
  }
  partable <- tryCatch(
    lavaan::lav_partable_unrestricted(fit),
    error = function(e) e
  )
  if (inherits(partable, "error")) {
    svysem_stop(
      paste0("The saturated model could not be built: ", conditionMessage(partable)),
      "svysem_error_data"
    )
  }
  missing <- fit@Options$missing
  saturated <- svysem_refit_partable(partable, data, weight, missing)
  saturated_pieces <- svysem_sandwich_pieces(saturated, design, numerical)
  trace_h1 <- svysem_trace_product(saturated, saturated_pieces)
  trace_h0 <- svysem_trace_product(fit, pieces)
  ndat <- ncol(saturated_pieces$H)
  npar <- lavaan::lavInspect(fit, "npar")
  fit@test <- svysem_set_yuan_bentler(
    fit@test, df, stat, trace_h1, trace_h0, ndat, npar
  )
  baseline <- fit@baseline$partable
  if (is.null(baseline)) {
    return(fit)
  }
  baseline_fit <- svysem_refit_partable(baseline, data, weight, missing)
  baseline_df <- baseline_fit@test[[1L]]$df
  baseline_stat <- baseline_fit@test[[1L]]$stat
  if (length(baseline_df) != 1L || !is.finite(baseline_df) || baseline_df <= 0) {
    return(fit)
  }
  baseline_pieces <- svysem_sandwich_pieces(baseline_fit, design, numerical)
  trace_baseline <- svysem_trace_product(baseline_fit, baseline_pieces)
  fit@baseline$test <- svysem_set_yuan_bentler(
    fit@baseline$test,
    baseline_df,
    baseline_stat,
    trace_h1,
    trace_baseline,
    ndat,
    lavaan::lavInspect(baseline_fit, "npar")
  )
  fit
}

svysem_fit_mlr <- function(model, design, fitter, dots) {
  ov <- svysem_observed_names(model)
  if (!length(ov)) {
    svysem_stop(
      "The model has no observed variables to pass to lavaan.",
      "svysem_error_variables"
    )
  }
  svysem_check_variables(ov, names(design$variables))
  model_columns <- design$variables[, ov, drop = FALSE]
  svysem_check_arguments(dots, model_columns, estimator = "MLR")
  non_numeric <- ov[!vapply(model_columns, is.numeric, logical(1))]
  if (length(non_numeric)) {
    svysem_stop(
      paste0(
        "Model variables must be numeric for MLR: ",
        paste(non_numeric, collapse = ", "), "."
      ),
      "svysem_error_data"
    )
  }
  # information = "observed" is the bread. se and test are lavaan's MLR
  # defaults; the design-based values replace them after the fit.
  numerical <- svysem_is_fiml(dots$missing)
  dots$se <- NULL
  dots$test <- NULL
  dots$do.fit <- NULL
  dots$information <- NULL
  prepared <- svysem_analysis_frame(
    design, ov, zero_incomplete = !numerical
  )
  fit <- do.call(
    fitter,
    svysem_lavaan_args(
      model, prepared$data, prepared$weight, dots,
      estimator = "MLR",
      extra = list(information = "observed")
    )
  )
  pieces <- svysem_sandwich_pieces(fit, prepared$design, numerical)
  fit <- svysem_inject_sandwich(fit, pieces)
  fit <- svysem_inject_yuan_bentler(
    fit, prepared$design, prepared$data, prepared$weight, numerical, pieces
  )
  list(fit = fit, design = design)
}

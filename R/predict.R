# Posterior mean of the latent variables given the observed continuous
# indicators: alpha + Psi Lambda' Sigma^{-1} (y - mu). Sigma, Lambda, Psi, and
# the means are the model-implied moments from the weighted estimates.
# Plausible values are draws from the corresponding posterior covariance
# Psi - Psi Lambda' Sigma^{-1} Lambda Psi. lavPredict() is not used: on a
# weighted fit in lavaan 0.6-21 it does not return this posterior mean.

svysem_unwrap <- function(x) {
  if (is.list(x) && !is.matrix(x) && length(x) == 1L) {
    return(x[[1L]])
  }
  x
}

svysem_as_named <- function(x, target_names, what) {
  values <- as.numeric(x)
  value_names <- names(x)
  if (is.null(value_names) && length(values) == length(target_names)) {
    value_names <- target_names
  }
  if (is.null(value_names) || !all(target_names %in% value_names)) {
    svysem_stop(
      paste0("Factor-score ", what, " could not be aligned with the model."),
      "svysem_error_data"
    )
  }
  stats::setNames(values[match(target_names, value_names)], target_names)
}

svysem_as_square <- function(x, target_names, what) {
  x <- as.matrix(svysem_unwrap(x))
  if (is.null(rownames(x)) && nrow(x) == length(target_names) &&
        ncol(x) == length(target_names)) {
    dimnames(x) <- list(target_names, target_names)
  }
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
        !all(target_names %in% rownames(x)) ||
        !all(target_names %in% colnames(x))) {
    svysem_stop(
      paste0("Factor-score ", what, " could not be aligned with the model."),
      "svysem_error_data"
    )
  }
  x[target_names, target_names, drop = FALSE]
}

svysem_factor_moments <- function(object) {
  estimates <- lavaan::lavInspect(object, "est")
  lambda <- svysem_unwrap(estimates$lambda)
  if (!is.matrix(lambda) || !nrow(lambda) || !ncol(lambda)) {
    svysem_stop(
      "The model has no latent variables to score.",
      "svysem_error_data"
    )
  }
  indicators <- rownames(lambda)
  latents <- colnames(lambda)
  if (is.null(indicators) || is.null(latents)) {
    svysem_stop(
      "Factor loadings are missing row or column names.",
      "svysem_error_data"
    )
  }
  psi <- svysem_as_square(lavaan::lavInspect(object, "cov.lv"), latents, "latent covariance")
  sigma <- svysem_as_square(lavaan::lavInspect(object, "cov.ov"), indicators, "observed covariance")
  alpha <- svysem_unwrap(lavaan::lavInspect(object, "mean.lv"))
  if (length(alpha) == 0L) {
    alpha <- stats::setNames(rep(0, length(latents)), latents)
  }
  alpha <- svysem_as_named(alpha, latents, "latent means")
  mu <- svysem_unwrap(lavaan::lavInspect(object, "mean.ov"))
  if (length(mu) == 0L) {
    # No mean structure: lavaan centers factor scores at the sample means.
    # Those means are not in lavInspect(fit, "mean.ov").
    mu <- object@SampleStats@mean[[1L]]
    observed <- lavaan::lavNames(object, "ov")
    if (is.null(names(mu)) && length(mu) == length(observed)) {
      names(mu) <- observed
    }
  }
  mu <- svysem_as_named(mu, indicators, "observed means")
  list(
    lambda = lambda[, latents, drop = FALSE],
    psi = psi,
    sigma = sigma,
    alpha = alpha,
    mu = mu,
    indicators = indicators,
    latents = latents
  )
}

svysem_check_draws <- function(draws) {
  if (length(draws) != 1L || !is.numeric(draws) || !is.finite(draws) ||
        draws < 0 || abs(draws - round(draws)) > 1e-8) {
    svysem_stop(
      "draws must be a single non-negative integer.",
      "svysem_error_data"
    )
  }
  as.integer(round(draws))
}

svysem_check_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed)) {
    svysem_stop(
      "seed must be a single finite number.",
      "svysem_error_data"
    )
  }
  set.seed(seed)
  invisible(NULL)
}

svysem_score_matrix <- function(object, newdata, indicators) {
  if (is.null(newdata)) {
    observed <- object@Data@X[[1L]]
    observed_names <- object@Data@ov.names[[1L]]
    if (is.null(colnames(observed))) {
      colnames(observed) <- observed_names
    }
    data <- as.data.frame(observed, stringsAsFactors = FALSE)
  } else {
    data <- as.data.frame(newdata, stringsAsFactors = FALSE, check.names = FALSE)
  }
  missing_names <- indicators[!indicators %in% names(data)]
  if (length(missing_names)) {
    svysem_stop(
      paste0(
        "Model variable(s) not found in the data: ",
        paste(missing_names, collapse = ", "), "."
      ),
      "svysem_error_variables"
    )
  }
  selected <- data[, indicators, drop = FALSE]
  if (!all(vapply(selected, is.numeric, logical(1)))) {
    svysem_stop(
      "Factor scores require numeric indicators.",
      "svysem_error_data"
    )
  }
  scores <- as.matrix(selected)
  storage.mode(scores) <- "double"
  rownames(scores) <- rownames(data)
  scores
}

svysem_posterior_map <- function(lambda, psi, sigma) {
  solved <- tryCatch(solve(sigma), error = function(e) NULL)
  if (is.null(solved) || any(!is.finite(solved))) {
    svysem_stop(
      "The model-implied covariance of the observed indicators is singular.",
      "svysem_error_data"
    )
  }
  gain <- psi %*% t(lambda) %*% solved
  variance <- psi - gain %*% lambda %*% psi
  variance <- (variance + t(variance)) / 2
  list(gain = gain, variance = variance)
}

svysem_clip_variance <- function(variance) {
  eig <- eigen(variance, symmetric = TRUE)
  if (any(eig$values < -1e-6)) {
    svysem_stop(
      "A factor-score posterior covariance is not positive semidefinite.",
      "svysem_error_data"
    )
  }
  eig$values[eig$values < 0] <- 0
  list(vectors = eig$vectors, values = eig$values)
}

svysem_draw_pattern <- function(mean, transform, n_draws) {
  n_rows <- nrow(mean)
  n_latents <- ncol(mean)
  draws <- vector("list", n_draws)
  for (draw in seq_len(n_draws)) {
    noise <- matrix(stats::rnorm(n_rows * n_latents), n_rows, n_latents)
    draws[[draw]] <- mean + noise %*% transform
  }
  draws
}

svysem_score_patterns <- function(observed, moments, n_draws) {
  finite <- is.finite(observed)
  pattern <- apply(finite, 1L, function(row) paste(as.integer(row), collapse = ""))
  n_rows <- nrow(observed)
  n_latents <- length(moments$latents)
  posterior <- matrix(
    NA_real_, n_rows, n_latents,
    dimnames = list(rownames(observed), moments$latents)
  )
  draw_rows <- vector("list", n_rows)
  for (key in unique(pattern)) {
    rows <- which(pattern == key)
    seen <- finite[rows[[1L]], ]
    if (!any(seen)) {
      next
    }
    mapped <- svysem_posterior_map(
      moments$lambda[seen, , drop = FALSE],
      moments$psi,
      moments$sigma[seen, seen, drop = FALSE]
    )
    deviation <- sweep(
      observed[rows, seen, drop = FALSE],
      2L, moments$mu[seen], "-"
    )
    mean <- matrix(moments$alpha, length(rows), n_latents, byrow = TRUE)
    mean <- mean + deviation %*% t(mapped$gain)
    posterior[rows, ] <- mean
    if (n_draws > 0L) {
      clipped <- svysem_clip_variance(mapped$variance)
      transform <- t(clipped$vectors %*% diag(sqrt(clipped$values), n_latents, n_latents))
      draw_rows[[rows[[1L]]]] <- list(
        rows = rows,
        draws = svysem_draw_pattern(mean, transform, n_draws)
      )
    }
  }
  if (n_draws == 0L) {
    return(posterior)
  }
  draws <- lapply(seq_len(n_draws), function(draw) {
    out <- matrix(
      NA_real_, n_rows, n_latents,
      dimnames = list(rownames(observed), moments$latents)
    )
    for (piece in draw_rows) {
      if (is.null(piece)) {
        next
      }
      out[piece$rows, ] <- piece$draws[[draw]]
    }
    out
  })
  draws
}

svysem_predict_one <- function(object, newdata, draws, seed) {
  ordered <- lavaan::lavNames(object, "ov.ord")
  if (length(ordered)) {
    svysem_stop(
      "Factor scores are for continuous outcomes. Ordered indicators are not scored.",
      "svysem_error_data"
    )
  }
  ngroups <- lavaan::lavInspect(object, "ngroups")
  if (length(ngroups) != 1L || ngroups != 1L) {
    svysem_stop(
      "Factor scores for more than one group are not implemented.",
      "svysem_error_data"
    )
  }
  draws <- svysem_check_draws(draws)
  svysem_check_seed(seed)
  moments <- svysem_factor_moments(object)
  observed <- svysem_score_matrix(object, newdata, moments$indicators)
  svysem_score_patterns(observed, moments, draws)
}

#' Factor scores and plausible values
#'
#' Posterior means of the latent variables in a [svysem-class] fit, using the
#' model-implied moments and the survey-weighted estimates. `draws = 0`
#' returns one row per case. A positive `draws` returns that many
#' plausible-value matrices, drawn from the posterior covariance of the same
#' regression. Cases are scored from the indicators they actually have, so a
#' full-information fit can leave some indicators missing.
#'
#' [lavaan::lavPredict()] on a weighted fit in lavaan 0.6-21 does not return
#' these values. The comparison in the examples is the unweighted formula,
#' which is what `predict()` implements.
#'
#' For a [svysem_mi-class] object the posterior means are averaged across
#' imputations. Plausible values are not averaged: set `per.imputation = TRUE`
#' to get one set of matrices per imputation. `newdata` is scored under each
#' imputation's estimates. Without `newdata`, each imputation is scored on the
#' data lavaan stored for that fit, and those rows are averaged in imputation
#' order, so the completed designs need the same row order.
#'
#' @param object A `svysem` or `svysem_mi` fit with continuous outcomes and one
#'   group.
#' @param newdata A data frame with the indicators used in the measurement
#'   model. `NULL` scores the cases stored in the fit.
#' @param draws Number of plausible-value data sets. `0` returns the posterior
#'   means.
#' @param seed Optional seed, passed to [set.seed()] before the draws. The
#'   previous random-generator state is not restored. `NULL` leaves it alone.
#' @param per.imputation If `TRUE`, return one result per imputation instead of
#'   the average posterior mean. Required when `draws` is positive.
#' @param ... Unused. Present so the method accepts the lavaan arguments it
#'   does not use.
#' @return For `draws = 0`, a numeric matrix with one row per case and one
#'   column per latent variable. For `draws > 0`, a list of that many
#'   matrices. A multiple-imputation fit returns the average of the posterior
#'   means, or a list of per-imputation results when `per.imputation = TRUE`.
#' @seealso [svysem()], [lavaan::lavPredict()]
#' @examples
#' design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 1)
#' fit <- svycfa("f =~ y1 + y2 + y3 + y4", design)
#' predict(fit)[1, , drop = FALSE]
#' @method predict svysem
#' @export
predict.svysem <- function(object, newdata = NULL, draws = 0, seed = NULL, ...) {
  svysem_predict_one(object, newdata = newdata, draws = draws, seed = seed)
}

#' @rdname predict.svysem
#' @export
setMethod("predict", "svysem", function(object, newdata = NULL, draws = 0,
                                        seed = NULL, ...) {
  svysem_predict_one(object, newdata = newdata, draws = draws, seed = seed)
})

#' @rdname predict.svysem
#' @method predict svysem_mi
#' @export
predict.svysem_mi <- function(object, newdata = NULL, draws = 0, seed = NULL,
                              per.imputation = FALSE, ...) {
  draws <- svysem_check_draws(draws)
  if (draws > 0L && !isTRUE(per.imputation)) {
    svysem_stop(
      paste0(
        "Plausible values are not averaged across imputations. ",
        "Set per.imputation = TRUE."
      ),
      "svysem_error_data"
    )
  }
  # One seed for the whole set, so imputation order is the only extra randomness.
  svysem_check_seed(seed)
  scores <- lapply(object@fits, function(fit) {
    svysem_predict_one(fit, newdata = newdata, draws = draws, seed = NULL)
  })
  if (isTRUE(per.imputation) || draws > 0L) {
    return(scores)
  }
  combined <- scores[[1L]]
  if (length(scores) == 1L) {
    return(combined)
  }
  for (imputation in scores[-1L]) {
    if (!identical(dim(imputation), dim(combined))) {
      svysem_stop(
        "Imputations did not return factor scores of the same dimension.",
        "svysem_error_data"
      )
    }
    combined <- combined + imputation
  }
  combined / length(scores)
}

#' @rdname predict.svysem
#' @export
setMethod("predict", "svysem_mi", function(object, newdata = NULL, draws = 0,
                                           seed = NULL, per.imputation = FALSE, ...) {
  predict.svysem_mi(
    object,
    newdata = newdata,
    draws = draws,
    seed = seed,
    per.imputation = per.imputation,
    ...
  )
})

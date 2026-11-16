library(lavaan)

svysem <- function(...) suppressWarnings(svysem::svysem(...))
svycfa <- function(...) suppressWarnings(svysem::svycfa(...))
svygrowth <- function(...) suppressWarnings(svysem::svygrowth(...))

cfa_model <- "f =~ y1 + y2 + y3 + y4"

posterior_mean <- function(fit) {
  estimates <- lavInspect(fit, "est")
  lambda <- estimates$lambda
  if (is.list(lambda) && !is.matrix(lambda)) {
    lambda <- lambda[[1L]]
  }
  psi <- lavInspect(fit, "cov.lv")
  sigma <- lavInspect(fit, "cov.ov")
  if (is.list(psi) && !is.matrix(psi)) {
    psi <- psi[[1L]]
  }
  if (is.list(sigma) && !is.matrix(sigma)) {
    sigma <- sigma[[1L]]
  }
  psi <- psi[colnames(lambda), colnames(lambda), drop = FALSE]
  sigma <- sigma[rownames(lambda), rownames(lambda), drop = FALSE]
  mu <- lavInspect(fit, "mean.ov")
  if (is.list(mu)) {
    mu <- mu[[1L]]
  }
  if (length(mu) == 0L) {
    mu <- fit@SampleStats@mean[[1L]]
    names(mu) <- lavNames(fit, "ov")
  }
  mu <- mu[rownames(lambda)]
  alpha <- lavInspect(fit, "mean.lv")
  if (is.list(alpha)) {
    alpha <- alpha[[1L]]
  }
  if (length(alpha) == 0L) {
    alpha <- rep(0, ncol(lambda))
  }
  alpha <- as.numeric(alpha)
  observed <- fit@Data@X[[1L]]
  if (is.null(colnames(observed))) {
    colnames(observed) <- lavNames(fit, "ov")
  }
  observed <- observed[, rownames(lambda), drop = FALSE]
  gain <- psi %*% t(lambda) %*% solve(sigma)
  deviation <- sweep(observed, 2L, mu, "-")
  mean <- matrix(alpha, nrow(observed), length(alpha), byrow = TRUE)
  mean <- mean + deviation %*% t(gain)
  colnames(mean) <- colnames(lambda)
  variance <- psi - gain %*% lambda %*% psi
  variance <- (variance + t(variance)) / 2
  list(mean = mean, variance = variance)
}

test_that("posterior means follow the regression formula, and lavPredict does not", {
  design <- simulate_design(n_psu = 8, n_per = 6, n_strata = 2, seed = 11)
  fit <- svycfa(cfa_model, design)
  scores <- predict(fit)
  manual <- posterior_mean(fit)
  expect_equal(unname(scores), unname(manual$mean), tolerance = 1e-8)
  expect_identical(colnames(scores), "f")

  lavaan_scores <- lavPredict(fit)
  gap <- max(abs(as.numeric(scores) - as.numeric(lavaan_scores)))
  if (packageVersion("lavaan") < "0.7") {
    expect_gt(gap, 1e-3)
  } else {
    expect_lt(gap, 1e-6)
  }

  rewritten <- fit
  rewritten@Data@weights[[1L]] <- rep(1, length(rewritten@Data@weights[[1L]]))
  expect_equal(
    as.numeric(lavPredict(rewritten)),
    as.numeric(scores),
    tolerance = 1e-8
  )
})

test_that("newdata is scored in its own row order, including partially missing rows", {
  design <- simulate_design(n_psu = 6, n_per = 5, n_strata = 2, seed = 12)
  fit <- svycfa(cfa_model, design)
  stored <- predict(fit)
  frame <- design$variables
  expect_equal(predict(fit, newdata = frame), stored, tolerance = 1e-8)

  partial <- frame
  partial$y2[1L] <- NA
  partial[2L, c("y1", "y2", "y3", "y4")] <- NA
  scored <- predict(fit, newdata = partial)
  expect_equal(scored[-c(1L, 2L), , drop = FALSE], stored[-c(1L, 2L), , drop = FALSE])
  expect_false(isTRUE(all.equal(scored[1L, ], stored[1L, ])))
  expect_true(is.finite(scored[1L, ]))
  expect_true(all(is.na(scored[2L, ])))

  expect_error(
    predict(fit, newdata = frame[, c("y1", "y2"), drop = FALSE]),
    class = "svysem_error_variables"
  )
})

test_that("plausible values stay near the posterior mean and covariance", {
  design <- simulate_design(n_psu = 6, n_per = 5, n_strata = 2, seed = 13)
  fit <- svycfa(cfa_model, design)
  mean <- predict(fit)
  manual <- posterior_mean(fit)
  draws <- predict(fit, draws = 400, seed = 21)
  expect_length(draws, 400)
  expect_true(is.matrix(predict(fit, draws = 0)))
  case <- vapply(draws, function(draw) draw[1L, 1L], numeric(1))
  expect_lt(abs(mean(case) - mean[1L, 1L]), 0.2)
  expect_lt(abs(stats::var(case) - manual$variance[1L, 1L]), 0.25)
  expect_equal(predict(fit, draws = 5, seed = 4), predict(fit, draws = 5, seed = 4))
  expect_error(predict(fit, draws = -1), class = "svysem_error_data")
  expect_error(predict(fit, seed = "a"), class = "svysem_error_data")
})

test_that("a full-information fit scores the indicators a case actually has", {
  design <- simulate_design(n_psu = 8, n_per = 5, n_strata = 2, seed = 14)
  design$variables$y2[1:5] <- NA
  design$variables$y4[6:8] <- NA
  fit <- svycfa(cfa_model, design, estimator = "MLR", missing = "ml")
  scores <- predict(fit)
  expect_identical(nrow(scores), nrow(design$variables))
  expect_true(all(is.finite(scores[1:8, ])))
  expect_false(isTRUE(all.equal(scores[1L, ], scores[6L, ])))
})

test_that("ordered indicators, several groups, and a model with no factors are refused", {
  design <- simulate_design(n_psu = 6, n_per = 8, n_strata = 2, nest = TRUE, seed = 15)
  ordered_data <- design$variables
  ordered_data$o1 <- ordered(ifelse(ordered_data$y1 > 0, 1, 0))
  ordered_data$o2 <- ordered(ifelse(ordered_data$y2 > 0, 1, 0))
  ordered_data$o3 <- ordered(ifelse(ordered_data$y3 > 0, 1, 0))
  ordered_design <- survey::svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, nest = TRUE, data = ordered_data
  )
  ordered_fit <- svycfa(
    "f =~ o1 + o2 + o3", ordered_design, estimator = "WLSMV"
  )
  expect_error(predict(ordered_fit), class = "svysem_error_data")

  grouped <- design$variables
  grouped$g <- ifelse(grouped$strat == 1, "a", "b")
  grouped_design <- survey::svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, nest = TRUE, data = grouped
  )
  grouped_fit <- svycfa(
    cfa_model, grouped_design, estimator = "WLSMV", group = "g"
  )
  expect_error(predict(grouped_fit), class = "svysem_error_data")

  regression <- svysem(
    "y ~ x1 + x2", design, fixed.x = FALSE, meanstructure = TRUE
  )
  expect_error(predict(regression), class = "svysem_error_data")
})

test_that("multiple imputation averages posterior means and does not average draws", {
  skip_if_not_installed("mitools")
  frames <- lapply(1:2, function(i) {
    simulate_design(n_psu = 4, n_per = 5, n_strata = 2, nest = TRUE, seed = 30L + i)$variables
  })
  design <- survey::svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, nest = TRUE,
    data = mitools::imputationList(frames)
  )
  fit <- svycfa(cfa_model, design)
  per_imputation <- predict(fit, per.imputation = TRUE)
  expect_length(per_imputation, 2)
  expect_equal(predict(fit), (per_imputation[[1L]] + per_imputation[[2L]]) / 2)
  expect_error(predict(fit, draws = 2), class = "svysem_error_data")
  drawn <- predict(fit, draws = 2, seed = 5, per.imputation = TRUE)
  expect_length(drawn, 2)
  expect_length(drawn[[1L]], 2)
})

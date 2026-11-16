library(lavaan)

svysem <- function(...) suppressWarnings(svysem::svysem(...))
svycfa <- function(...) suppressWarnings(svysem::svycfa(...))
svygrowth <- function(...) suppressWarnings(svysem::svygrowth(...))

cfa_model <- "f =~ y1 + y2 + y3 + y4"

iid_frame <- function(n, seed) {
  simulate_design(n_psu = n, n_per = 1L, n_strata = 0L, seed = seed)$variables
}

iid_design <- function(data) {
  survey::svydesign(ids = ~0, weights = ~w, data = data)
}

free_se <- function(fit) {
  estimates <- parameterEstimates(fit, se = TRUE)
  keep <- estimates$free > 0L
  list(
    table = estimates$se[keep],
    covariance = sqrt(diag(vcov(fit)))[estimates$free[keep]]
  )
}

test_that("independent MLR standard errors and the Yuan-Bentler statistic match lavaan", {
  data <- iid_frame(250L, 2L)
  design <- iid_design(data)
  fit <- svycfa(cfa_model, design, estimator = "MLR", se = "none", test = "standard")
  lavaan_fit <- cfa(cfa_model, data = data, estimator = "MLR", sampling.weights = "w")
  expect_s4_class(fit, "svysem")
  expect_identical(fit@estimator, "MLR")
  expect_identical(fit@variance_method, "linearization")
  expect_rel(coef(fit), coef(lavaan_fit), 1e-6)
  expect_rel(sqrt(diag(vcov(fit))), sqrt(diag(vcov(lavaan_fit))), 0.01)
  expect_true(isSymmetric(vcov(fit)))
  expect_identical(colnames(vcov(fit)), names(coef(fit)))

  standard_errors <- free_se(fit)
  expect_equal(
    unname(standard_errors$table),
    unname(standard_errors$covariance),
    tolerance = 1e-8
  )

  measures <- fitMeasures(fit, c("chisq.scaled", "chisq.scaling.factor", "cfi.scaled"))
  lavaan_measures <- fitMeasures(
    lavaan_fit, c("chisq.scaled", "chisq.scaling.factor", "cfi.scaled")
  )
  expect_true(all(is.finite(measures)))
  expect_rel(measures, lavaan_measures, 0.05)
  expect_true(any(vapply(fit@test, function(element) {
    identical(element$test, "yuan.bentler.mplus")
  }, logical(1))))

  standardized <- standardizedSolution(fit, se = TRUE)
  expect_true(any(is.finite(standardized$se)))
})

test_that("a defined parameter uses the delta method and lavTestWald reads it", {
  design <- simulate_design(n_psu = 8, n_per = 6, n_strata = 2, seed = 8)
  model <- "f =~ y1 + y2 + y3 + y4\ny2 ~~ v2*y2\ny3 ~~ v3*y3\nind := v2 + v3"
  fit <- svycfa(model, design, estimator = "MLR")
  covariance <- vcov(fit)
  delta <- sqrt(
    covariance["v2", "v2"] + covariance["v3", "v3"] + 2 * covariance["v2", "v3"]
  )
  estimates <- parameterEstimates(fit, se = TRUE)
  expect_equal(estimates$se[estimates$lhs == "ind"], unname(delta), tolerance = 1e-8)
  wald <- lavTestWald(fit, "ind == 0")
  expect_true(is.finite(wald$stat))
  expect_true(is.finite(wald$p.value))
})

test_that("full-information MLR keeps every row and returns a finite sandwich", {
  data <- iid_frame(180L, 9L)
  data$y2[seq(1L, 180L, by = 7L)] <- NA
  data$y4[seq(3L, 180L, by = 11L)] <- NA
  design <- iid_design(data)
  fit <- svycfa(cfa_model, design, estimator = "MLR", missing = "ml")
  lavaan_fit <- cfa(
    cfa_model, data = data, estimator = "MLR",
    sampling.weights = "w", missing = "ml"
  )
  expect_equal(as.numeric(lavInspect(fit, "nobs")), nrow(data))
  expect_gt(as.numeric(lavInspect(fit, "nobs")), sum(stats::complete.cases(data[c("y1", "y2", "y3", "y4")])))
  expect_rel(coef(fit), coef(lavaan_fit), 1e-5)
  expect_rel(sqrt(diag(vcov(fit))), sqrt(diag(vcov(lavaan_fit))), 0.02)
  expect_true(isSymmetric(unname(vcov(fit))))
  expect_true(all(is.finite(vcov(fit))))
  expect_true(all(is.finite(predict(fit))))

  alias <- svycfa(cfa_model, design, estimator = "MLR", missing = "fiml")
  expect_rel(coef(alias), coef(fit), 1e-8)
})

test_that("MLR refuses ordered data, groups, pairwise missing, and replicate controls", {
  design <- simulate_design(n_psu = 6, n_per = 5, n_strata = 2, seed = 10)
  expect_error(
    svycfa(cfa_model, design, estimator = "MLR", ordered = "y1"),
    class = "svysem_error_data"
  )
  expect_error(
    svycfa(cfa_model, design, estimator = "MLR", group = "strat"),
    class = "svysem_error_data"
  )
  expect_error(
    svycfa(cfa_model, design, estimator = "MLR", missing = "pairwise"),
    class = "svysem_error_data"
  )
  expect_error(
    svycfa(cfa_model, design, estimator = "MLR", conditional.x = TRUE),
    class = "svysem_error_data"
  )
  expect_error(
    svycfa(cfa_model, design, estimator = "MLR", replicates = 4),
    class = "svysem_error_design"
  )
  replicate <- survey::as.svrepdesign(design, type = "auto")
  replicated <- svycfa(cfa_model, replicate, estimator = "MLR")
  expect_identical(replicated@variance_method, "replication")
  expect_true(all(is.finite(sqrt(diag(vcov(replicated))))))
  expect_error(
    svycfa(cfa_model, replicate, estimator = "MLR", replicate.type = "JKn"),
    class = "svysem_error_design"
  )
})

test_that("multiple imputation accepts an MLR fit", {
  skip_if_not_installed("mitools")
  frames <- lapply(1:2, function(i) {
    simulate_design(n_psu = 4, n_per = 6, n_strata = 2, nest = TRUE, seed = 40L + i)$variables
  })
  design <- survey::svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, nest = TRUE,
    data = mitools::imputationList(frames)
  )
  fit <- svycfa("f =~ y1 + y2 + y3", design, estimator = "MLR")
  expect_s4_class(fit, "svysem_mi")
  expect_identical(fit@estimator, "MLR")
  expect_length(fit@fits, 2)
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(predict(fit))))
})

replicate_se <- function(design, model, missing = NULL, ...) {
  replicate <- survey::as.svrepdesign(design, type = "JKn", compress = FALSE)
  statistic <- survey::withReplicates(replicate, function(w, data) {
    data$.w <- as.numeric(w)
    args <- list(
      model = model,
      data = data,
      estimator = "ML",
      se = "none",
      sampling.weights = ".w",
      ...
    )
    if (!is.null(missing)) {
      args$missing <- missing
    }
    lavaan::coef(do.call(lavaan::sem, args))
  })
  variance <- stats::vcov(statistic)
  sqrt(diag(variance))
}

test_that("linearization matches a jackknife of the same design", {
  skip_on_cran()
  # 160 primary sampling units. At 100, a full-information residual
  # variance is still about 2.1 percent below its jackknife standard error.
  # At 160 the same comparison is under 1 percent.
  design <- simulate_design(
    n_psu = 80, n_per = 4, n_strata = 2, nest = TRUE, seed = 5
  )
  regression <- "y ~ x1 + x2"
  mlm <- svysem(regression, design, fixed.x = FALSE, meanstructure = TRUE)
  mlm_oracle <- replicate_se(
    design, regression, fixed.x = FALSE, meanstructure = TRUE
  )
  expect_rel(sqrt(diag(vcov(mlm)))[names(mlm_oracle)], mlm_oracle, 0.02)

  cfa_fit <- svycfa(cfa_model, design, estimator = "MLR")
  cfa_oracle <- replicate_se(design, cfa_model)
  expect_rel(sqrt(diag(vcov(cfa_fit)))[names(cfa_oracle)], cfa_oracle, 0.02)

  design$variables$y2[seq(1L, nrow(design$variables), by = 6L)] <- NA
  design$variables$y4[seq(2L, nrow(design$variables), by = 8L)] <- NA
  fiml <- svycfa(cfa_model, design, estimator = "MLR", missing = "ml")
  fiml_oracle <- replicate_se(design, cfa_model, missing = "ml")
  expect_equal(as.numeric(lavInspect(fiml, "nobs")), nrow(design$variables))
  expect_rel(sqrt(diag(vcov(fiml)))[names(fiml_oracle)], fiml_oracle, 0.02)
})

test_that("the MLR covariance is the hand-built sandwich and numerical scores match lavScores", {
  data <- iid_frame(80L, 4L)
  design <- iid_design(data)
  fit <- svycfa(cfa_model, design, estimator = "MLR")
  pieces <- svysem:::svysem_sandwich_pieces(fit, design, numerical = FALSE)
  bread <- solve(pieces$H)
  built <- bread %*% pieces$meat %*% bread
  expect_equal((built + t(built)) / 2, vcov(fit), tolerance = 1e-8)

  numerical <- svysem:::svysem_numerical_scores(fit)
  analytic <- lavaan::lavScores(fit, remove.empty.cases = FALSE)
  expect_equal(dim(numerical), dim(analytic))
  if (!is.null(colnames(analytic))) {
    colnames(numerical) <- colnames(analytic)
  }
  weights <- svysem:::svysem_lavaan_weight_vector(fit, nrow(data))
  # lavaan 0.7 already returns weight * score. Earlier versions do not.
  expected <- if (isTRUE(svysem:::svysem_lavscores_weighted())) {
    analytic
  } else {
    analytic * weights
  }
  expect_equal(numerical, expected, tolerance = 1e-6)
})

test_that("unstable FIML scores are an error and do not fall back to complete cases", {
  # A converged score sums to about zero, so the guard is checked with a
  # constructed total. The full-information path calls this before svytotal
  # and has no complete-case branch.
  unbalanced <- matrix(c(1, 1, 1, 1), ncol = 1)
  expect_error(
    svysem:::svysem_check_score_balance(unbalanced),
    "does not fall back to complete cases",
    class = "svysem_error_data"
  )
  balanced <- matrix(c(1, -1), ncol = 1)
  expect_null(svysem:::svysem_check_score_balance(balanced))
})

test_that("listwise MLR keeps incomplete rows and matches complete-case estimates", {
  design <- simulate_design(
    n_psu = 6, n_per = 8, n_strata = 2, nest = TRUE, seed = 4
  )
  design$variables$y1[c(1, 9)] <- NA
  fit <- svycfa(cfa_model, design, estimator = "MLR")
  expect_s4_class(fit, "svysem")
  expect_identical(fit@estimator, "MLR")
  expect_identical(fit@variance_method, "linearization")
  expect_true(isSymmetric(vcov(fit)))
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))

  oracle_data <- design$variables
  oracle_data$w_oracle <- as.numeric(weights(design, "sampling"))
  oracle <- cfa(
    cfa_model,
    data = oracle_data,
    estimator = "MLR",
    sampling.weights = "w_oracle",
    missing = "listwise"
  )
  expect_rel(coef(fit), coef(oracle), 1e-6)
  expect_equal(lavaan::lavInspect(fit, "nobs"), lavaan::lavInspect(oracle, "nobs"))

  replicated <- svycfa(
    cfa_model, survey::as.svrepdesign(design), estimator = "MLR"
  )
  expect_identical(replicated@variance_method, "replication")
  expect_rel(coef(replicated), coef(fit), 1e-6)
  expect_true(all(is.finite(sqrt(diag(vcov(replicated))))))
})

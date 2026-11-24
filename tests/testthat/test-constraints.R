library(lavaan)

svycfa <- function(...) suppressWarnings(svysem::svycfa(...))

srs_sample <- function(n, seed) {
  data <- simulate_design(n_psu = n, n_per = 1L, n_strata = 0L, seed = seed)$variables
  list(
    data = data,
    design = survey::svydesign(ids = ~0, weights = ~w, data = data)
  )
}

# coef() repeats a label under an equality constraint. expect_rel() would
# keep only the first column of that name, so the comparison stays positional.
expect_positional <- function(got, expected, tolerance) {
  expect_rel(unname(as.numeric(got)), unname(as.numeric(expected)), tolerance)
}

expect_constraint_mlr <- function(model, ..., seed = 21L) {
  sample <- srs_sample(300L, seed)
  fit <- svycfa(model, sample$design, estimator = "MLR", ...)
  oracle <- lavaan::cfa(
    model, data = sample$data, estimator = "MLR",
    sampling.weights = "w", ...
  )
  expect_s4_class(fit, "svysem")
  expect_positional(coef(fit), coef(oracle), 1e-6)
  expect_positional(sqrt(diag(vcov(fit))), sqrt(diag(vcov(oracle))), 0.01)
  expect_true(isSymmetric(unname(vcov(fit))))
  measures <- lavaan::fitMeasures(fit, c("chisq.scaled", "df.scaled"))
  oracle_measures <- lavaan::fitMeasures(oracle, c("chisq.scaled", "df.scaled"))
  expect_true(all(is.finite(measures)))
  expect_positional(measures, oracle_measures, 0.05)
  estimates <- parameterEstimates(fit, se = TRUE)
  free <- estimates$free > 0L
  expect_equal(
    unname(estimates$se[free]),
    unname(sqrt(diag(vcov(fit)))[estimates$free[free]]),
    tolerance = 1e-8
  )
}

test_that("MLR standard errors match lavaan when a label is shared", {
  expect_constraint_mlr("f =~ y1 + a*y2 + a*y3 + y4")
})

test_that("MLR standard errors match lavaan for an explicit equality", {
  expect_constraint_mlr("f =~ y1 + a*y2 + b*y3 + y4\na == b")
})

test_that("MLR standard errors match lavaan for a nonlinear equality", {
  expect_constraint_mlr("f =~ y1 + a*y2 + b*y3 + y4\na == b^2")
})

test_that("MLR standard errors match lavaan for effect coding and std.lv", {
  expect_constraint_mlr("f =~ y1 + y2 + y3 + y4", effect.coding = TRUE, seed = 22L)
  expect_constraint_mlr("f =~ y1 + y2 + y3 + y4", std.lv = TRUE, seed = 23L)
})

test_that("a defined parameter keeps the delta method beside an equality", {
  sample <- srs_sample(240L, 24L)
  model <- "f =~ y1 + a*y2 + a*y3 + y4\ny4 ~~ v*y4\ntwice := 2*v"
  fit <- svycfa(model, sample$design, estimator = "MLR")
  estimates <- parameterEstimates(fit, se = TRUE)
  expect_equal(
    estimates$se[estimates$lhs == "twice"],
    2 * estimates$se[estimates$label == "v"],
    tolerance = 1e-8
  )
  expect_true(is.finite(lavTestWald(fit, "twice == 0")$stat))
})

# Test 11. svyimputationList method: parallel fits, stop, and omit.

mi_design <- function(n_imp = 4, fail_at = 0L, zero_weight = FALSE, shift_strata = FALSE) {
  frames <- lapply(seq_len(n_imp), function(i) {
    n <- 20
    strata <- if (shift_strata && i == 2L) rep(1:4, each = 5) else rep(1:2, each = 10)
    w <- rep(1, n)
    if (zero_weight) {
      w[[1]] <- 0
    }
    data.frame(
      y = (seq_len(n) / 10) + sin(seq_len(n)) + i * 0.05,
      x = seq_len(n) / 5,
      w = w,
      id = seq_len(n),
      strata = strata,
      fail = i == fail_at
    )
  })
  survey::svydesign(
    ids = ~id, strata = ~strata, weights = ~w,
    data = mitools::imputationList(frames)
  )
}

numeric_stub <- function(model, design, estimator = "MLM", ...) {
  vars <- design$variables
  if (isTRUE(vars$fail[[1]])) {
    stop("stub explosion", call. = FALSE)
  }
  w <- as.numeric(stats::weights(design))
  x <- vars$x
  y <- vars$y
  b <- sum(w * x * y) / sum(w * x * x)
  s2 <- sum(w * (y - b * x)^2) / sum(w)
  cf <- c(b1 = unname(b), b2 = unname(s2))
  V <- diag(c(1 / sum(w), 1 / sum(w)^2), nrow = 2)
  dimnames(V) <- list(names(cf), names(cf))
  structure(list(coefficients = cf, variance = V), class = "svysem_fake")
}

test_that("parallel and sequential fits are identical", {
  skip_if_not_installed("mitools")
  skip_if_not_installed("parallel")
  des <- mi_design()
  restore <- register_design_method("svysem", numeric_stub)
  on.exit(restore(), add = TRUE)
  sequential <- svysem("y ~ x", des, estimator = "MLM")
  parallel <- svysem(
    "y ~ x", des, estimator = "MLM",
    parallel = function(X, FUN) parallel::mclapply(X, FUN, mc.cores = 2)
  )
  expect_s4_class(sequential, "svysem_mi")
  expect_equal(coef(parallel), coef(sequential), tolerance = 0)
  expect_equal(vcov(parallel), vcov(sequential), tolerance = 0)
  expect_equal(nobs(sequential), 20)
  expect_equal(sequential@df, survey::degf(des$designs[[1]]))
  expect_identical(
    names(sequential@status),
    c("no.conv", "no.se", "no.npd", "message")
  )
  expect_false(any(sequential@status$no.conv))
  expect_false(any(sequential@status$no.se))
  expect_true(all(sequential@status$message == ""))
  expect_true(is.call(sequential@svysem_call))
  pe <- pooled_estimates(sequential)
  manual <- pooled_estimates(sequential@fits, df.com = sequential@df)
  expect_equal(pe, manual, tolerance = 1e-12)
  over <- pooled_estimates(sequential, df.com = NULL)
  expect_equal(over$se, pe$se, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(over$df, pe$df)))
  expect_equal(unname(sqrt(diag(vcov(sequential)))), unname(pe$se), tolerance = 1e-10)
  expect_output(show(sequential), "svysem_mi")
  expect_output(print(summary(sequential)), "no.conv")
  expect_equal(lavaan::coef(sequential), coef(sequential))
  expect_equal(lavaan::nobs(sequential), nobs(sequential))
})

test_that("a zero weight is not an observation", {
  skip_if_not_installed("mitools")
  des <- mi_design(zero_weight = TRUE)
  restore <- register_design_method("svysem", numeric_stub)
  on.exit(restore(), add = TRUE)
  fit <- svysem("y ~ x", des)
  expect_equal(nobs(fit), 19)
})

test_that("on_fail stop signals the imputation index", {
  skip_if_not_installed("mitools")
  des <- mi_design(fail_at = 2L)
  restore <- register_design_method("svysem", numeric_stub)
  on.exit(restore(), add = TRUE)
  err <- expect_error(
    svysem("y ~ x", des, on_fail = "stop"),
    class = "svysem_error_imputation"
  )
  expect_s3_class(err, "svysem_error")
  expect_match(conditionMessage(err), "2")
  expect_match(conditionMessage(err), "stub explosion")
})

test_that("on_fail omit drops the failure and pools the rest", {
  skip_if_not_installed("mitools")
  des <- mi_design(n_imp = 4, fail_at = 2L)
  restore <- register_design_method("svysem", numeric_stub)
  on.exit(restore(), add = TRUE)
  fit <- svysem("y ~ x", des, on_fail = "omit")
  expect_equal(length(fit@fits), 3)
  expect_equal(nrow(fit@status), 4)
  expect_true(fit@status$no.conv[[2]])
  expect_match(fit@status$message[[2]], "stub explosion")
  expect_equal(fit@status$message[[1]], "")
  manual <- lapply(des$designs[c(1, 3, 4)], function(d) numeric_stub("y ~ x", d))
  expect_equal(
    unname(coef(fit)),
    unname(pooled_estimates(manual, df.com = fit@df)$est),
    tolerance = 1e-12
  )
  expect_equal(unique(pooled_estimates(fit)$m), 3)
})

test_that("omit with fewer than two fits cannot pool", {
  skip_if_not_installed("mitools")
  des <- mi_design(n_imp = 2, fail_at = 1L)
  restore <- register_design_method("svysem", numeric_stub)
  on.exit(restore(), add = TRUE)
  err <- expect_error(
    svysem("y ~ x", des, on_fail = "omit"),
    class = "svysem_error_pooling"
  )
  expect_s3_class(err, "svysem_error")
})

test_that("completed designs must share degrees of freedom", {
  skip_if_not_installed("mitools")
  des <- mi_design(n_imp = 2, shift_strata = TRUE)
  restore <- register_design_method("svysem", numeric_stub)
  on.exit(restore(), add = TRUE)
  expect_error(svysem("y ~ x", des), class = "svysem_error_design")
})

test_that("estimators other than MLM are rejected", {
  skip_if_not_installed("mitools")
  des <- mi_design(n_imp = 2)
  err <- expect_error(
    svysem("y ~ x", des, estimator = "WLSMV"),
    class = "svysem_error_estimator"
  )
  expect_s3_class(err, "svysem_error")
  expect_error(
    svycfa("y ~ x", des, estimator = "MLR"),
    class = "svysem_error_estimator"
  )
})

test_that("svycfa and svygrowth call their own generics", {
  skip_if_not_installed("mitools")
  des <- mi_design(n_imp = 2)
  tag_stub <- function(tag) {
    function(model, design, estimator = "MLM", ...) {
      cf <- c(tag = tag)
      V <- matrix(0.01, 1, 1, dimnames = list("tag", "tag"))
      structure(list(coefficients = cf, variance = V), class = "svysem_fake")
    }
  }
  restore <- list(
    register_design_method("svysem", tag_stub(1)),
    register_design_method("svycfa", tag_stub(2)),
    register_design_method("svygrowth", tag_stub(3))
  )
  on.exit(lapply(restore, function(fn) fn()), add = TRUE)
  expect_equal(unname(coef(svysem("y ~ x", des))), 1)
  expect_equal(unname(coef(svycfa("y ~ x", des))), 2)
  expect_equal(unname(coef(svygrowth("y ~ x", des))), 3)
})

test_that("a lavaan stub records convergence and is pooled with lavaan coef", {
  skip_if_not_installed("mitools")
  des <- mi_design(n_imp = 3)
  lavaan_stub <- function(model, design, estimator = "MLM", ...) {
    lavaan::sem("y ~ x", data = design$variables)
  }
  restore <- register_design_method("svysem", lavaan_stub)
  on.exit(restore(), add = TRUE)
  fit <- svysem("y ~ x", des)
  expect_false(any(fit@status$no.conv))
  expect_false(any(fit@status$no.se))
  expect_false(any(fit@status$no.npd))
  manual <- lapply(des$designs, function(d) lavaan_stub("y ~ x", d))
  expect_equal(
    unname(coef(fit)),
    unname(pooled_estimates(manual, df.com = fit@df)$est),
    tolerance = 1e-10
  )
  nose <- function(model, design, estimator = "MLM", ...) {
    lavaan::sem("y ~ x", data = design$variables, se = "none")
  }
  restore_nose <- register_design_method("svysem", nose)
  on.exit(restore_nose(), add = TRUE)
  quiet <- svysem("y ~ x", des)
  expect_true(all(quiet@status$no.se))
  expect_false(any(quiet@status$no.conv))
  expect_equal(quiet@status$message, rep("", 3))
})

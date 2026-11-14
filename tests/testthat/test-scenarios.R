library(survey)
library(lavaan)

# Saturated MLM fits make lavaan warn that a parameter-covariance eigenvalue
# is about -1e-17. The standard errors still match svyglm; the warning is
# numerical noise, so the wrappers below keep it out of the test log.
svysem <- function(...) suppressWarnings(svysem::svysem(...))
svycfa <- function(...) suppressWarnings(svysem::svycfa(...))
svygrowth <- function(...) suppressWarnings(svysem::svygrowth(...))

expect_rel_equal <- function(got, expected, tolerance = 1e-6) {
  got <- as.numeric(got)
  expected <- as.numeric(expected)
  expect_identical(length(got), length(expected))
  scale <- pmax(abs(got), abs(expected), .Machine$double.eps)
  expect_lte(max(abs(got - expected) / scale), tolerance)
}

expect_classed <- function(expr, class) {
  err <- expect_error(expr, class = class)
  expect_s3_class(err, "svysem_error")
  expect_identical(class(err)[[1L]], class)
  err
}

expect_matches_glm <- function(design) {
  fit <- svysem(
    "y ~ x1 + x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  got <- regression_parts(fit)
  oracle <- glm_parts(design)
  expect_rel_equal(got$est, oracle$est, 1e-6)
  expect_rel_equal(got$se, oracle$se, 1e-6)
  expect_true(all(is.finite(got$est)))
  expect_true(all(is.finite(got$se)))
  invisible(fit)
}

small_data <- function(seed) {
  simulate_design(
    n_psu = 2, n_per = 5, n_strata = 4, nest = TRUE, seed = seed
  )$variables
}

cluster_design <- function(seed, nest = TRUE) {
  simulate_design(
    n_psu = 2, n_per = 5, n_strata = 4, nest = nest, seed = seed
  )
}

mi_frames <- function(m, seed, bad = 0L, zero = FALSE) {
  lapply(seq_len(m), function(i) {
    set.seed(seed + i)
    n <- 40L
    eta <- stats::rnorm(n)
    x <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    w <- stats::runif(n, 0.7, 1.6)
    if (zero) {
      w[c(1L, 2L)] <- 0
    }
    if (i == bad) {
      w[[1L]] <- -1
    }
    data.frame(
      y = 0.5 * x - 0.2 * x2 + 0.05 * i + stats::rnorm(n),
      x = x,
      x2 = x2,
      y1 = eta + stats::rnorm(n, 0, 0.3),
      y2 = eta + stats::rnorm(n, 0, 0.3),
      y3 = eta + stats::rnorm(n, 0, 0.3),
      y4 = eta + stats::rnorm(n, 0, 0.3),
      w = w,
      id = rep(1:8, each = 5),
      strata = rep(1:4, each = 10)
    )
  })
}

mi_design <- function(m, seed, bad = 0L, zero = FALSE) {
  survey::svydesign(
    ids = ~id, strata = ~strata, weights = ~w, nest = TRUE,
    data = mitools::imputationList(mi_frames(m, seed, bad = bad, zero = zero))
  )
}

test_that("an unclustered design with no strata matches svyglm, including the intercept", {
  data <- small_data(1)
  for (ids in list(~0, ~1)) {
    design <- svydesign(ids = ids, weights = ~w, data = data)
    expect_matches_glm(design)
  }
})

test_that("a stratified unclustered design matches svyglm", {
  data <- small_data(2)
  design <- svydesign(ids = ~1, strata = ~strat, weights = ~w, data = data)
  expect_matches_glm(design)
})

test_that("nest = TRUE with primary sampling unit labels reused across strata matches svyglm", {
  design <- cluster_design(3, nest = TRUE)
  expect_true(length(unique(design$variables$psu)) < nrow(unique(design$variables[, c("strat", "psu")])))
  expect_matches_glm(design)
})

test_that("unbalanced cluster sizes match svyglm", {
  set.seed(4)
  spec <- data.frame(
    strat = rep(1:4, each = 2),
    psu = rep(1:2, 4),
    n = rep(c(2L, 8L), 4)
  )
  pieces <- lapply(seq_len(nrow(spec)), function(i) {
    k <- spec$n[[i]]
    data.frame(
      strat = spec$strat[[i]],
      psu = spec$psu[[i]],
      x1 = stats::rnorm(k),
      x2 = stats::rnorm(k),
      w = stats::runif(k, 0.5, 2)
    )
  })
  data <- do.call(rbind, pieces)
  data$y <- stats::rnorm(nrow(data)) + 0.3 * data$x1
  design <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = data, nest = TRUE
  )
  expect_identical(sort(unique(as.integer(table(data$strat, data$psu)))), c(2L, 8L))
  expect_matches_glm(design)
})

test_that("strata stored as a factor and as a character give the same estimates", {
  data <- small_data(5)
  as_character <- data
  as_character$strat <- as.character(as_character$strat)
  as_factor <- data
  as_factor$strat <- factor(
    as.character(as_factor$strat),
    levels = as.character(sort(unique(data$strat)))
  )
  char_design <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = as_character, nest = TRUE
  )
  factor_design <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = as_factor, nest = TRUE
  )
  char_fit <- svysem("y ~ x1", char_design, fixed.x = FALSE, meanstructure = TRUE)
  factor_fit <- svysem("y ~ x1", factor_design, fixed.x = FALSE, meanstructure = TRUE)
  expect_rel_equal(coef(char_fit), coef(factor_fit), 1e-8)
  expect_rel_equal(sqrt(diag(vcov(char_fit))), sqrt(diag(vcov(factor_fit))), 1e-8)
})

test_that("a domain via subset() has finite standard errors that differ from a rebuilt design", {
  # Three PSUs per stratum: dropping one label leaves two, so the rebuilt
  # design is not a lonely-PSU design.
  base <- simulate_design(n_psu = 3, n_per = 4, n_strata = 4, nest = TRUE, seed = 6)
  design <- subset(base, psu != 1)
  kept <- base$variables[base$variables$psu != 1, , drop = FALSE]
  rebuilt <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = kept, nest = TRUE
  )
  fit_domain <- svysem("y ~ x1 + x2", design, fixed.x = FALSE, meanstructure = TRUE)
  fit_rebuilt <- svysem("y ~ x1 + x2", rebuilt, fixed.x = FALSE, meanstructure = TRUE)
  se_domain <- regression_parts(fit_domain)$se
  se_rebuilt <- regression_parts(fit_rebuilt)$se
  expect_true(all(is.finite(se_domain)))
  relative <- abs(se_domain - se_rebuilt) / pmax(abs(se_rebuilt), 1e-12)
  expect_gt(max(relative), 1e-3)
})

test_that("dropping one stratum leaves a lonely PSU and does not change survey.lonely.psu", {
  set.seed(7)
  blocks <- list(
    c(strat = 1, psu = 1, n = 5),
    c(strat = 2, psu = 1, n = 4),
    c(strat = 2, psu = 2, n = 4),
    c(strat = 3, psu = 1, n = 4),
    c(strat = 3, psu = 2, n = 4),
    c(strat = 4, psu = 1, n = 4),
    c(strat = 4, psu = 2, n = 4)
  )
  pieces <- lapply(blocks, function(block) {
    k <- block[["n"]]
    data.frame(
      strat = block[["strat"]],
      psu = block[["psu"]],
      y = stats::rnorm(k),
      x1 = stats::rnorm(k),
      w = stats::runif(k, 0.5, 2)
    )
  })
  data <- do.call(rbind, pieces)
  design <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = data, nest = TRUE
  )
  sub <- subset(design, strat != 4)
  previous <- options(survey.lonely.psu = getOption("survey.lonely.psu"))
  on.exit(options(previous), add = TRUE)
  n_ok <- 0L
  for (setting in list(NULL, "fail", "remove", "adjust", "average", "certainty")) {
    options(survey.lonely.psu = setting)
    means <- tryCatch(svymean(~y, sub), error = function(e) e)
    if (inherits(means, "error")) {
      expect_identical(getOption("survey.lonely.psu"), setting)
      next
    }
    fit <- svysem("y ~ 1", sub, fixed.x = FALSE, meanstructure = TRUE)
    expect_rel_equal(coef(fit)[["y~1"]], coef(means), 1e-6)
    expect_rel_equal(sqrt(vcov(fit)["y~1", "y~1"]), SE(means), 1e-6)
    expect_identical(getOption("survey.lonely.psu"), setting)
    n_ok <- n_ok + 1L
  }
  expect_gt(n_ok, 0L)
})

test_that("weights of very different scales are invariant to multiplying every weight", {
  set.seed(8)
  n <- 40L
  data <- data.frame(
    y = stats::rnorm(n),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    strat = rep(1:4, each = 10),
    psu = rep(rep(1:2, each = 5), 4),
    w = c(stats::runif(20, 5e-4, 2e-3), stats::runif(20, 5e2, 2e3))
  )
  data$y <- data$y + 0.5 * data$x1
  fit_one <- function(scale) {
    scaled <- data
    scaled$w <- scaled$w * scale
    design <- svydesign(
      ids = ~psu, strata = ~strat, weights = ~w, data = scaled, nest = TRUE
    )
    svysem("y ~ x1 + x2", design, fixed.x = FALSE, meanstructure = TRUE)
  }
  base <- fit_one(1)
  expect_true(all(is.finite(coef(base))))
  expect_true(all(is.finite(vcov(base))))
  for (scale in c(1e-4, 1e4)) {
    scaled <- fit_one(scale)
    expect_rel_equal(coef(scaled), coef(base), 1e-8)
    expect_rel_equal(sqrt(diag(vcov(scaled))), sqrt(diag(vcov(base))), 1e-8)
  }
})

test_that("zero weights do not change point estimates, and imputation nobs ignores them", {
  data <- small_data(9)
  data$w[c(1, 2, 6, 11)] <- 0
  with_zeros <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = data, nest = TRUE
  )
  kept <- data[data$w > 0, , drop = FALSE]
  dropped <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = kept, nest = TRUE
  )
  fit_zeros <- svysem("y ~ x1", with_zeros, fixed.x = FALSE, meanstructure = TRUE)
  fit_dropped <- svysem("y ~ x1", dropped, fixed.x = FALSE, meanstructure = TRUE)
  expect_rel_equal(coef(fit_zeros), coef(fit_dropped), 1e-8)
  n_pos <- sum(data$w > 0)
  # The single-design object has no nobs() method of its own. lavaan's nobs
  # is the one that runs. Observed with lavaan 0.6-21: it equals the number
  # of rows (40 here), so the four zero-weight rows are counted. nobs() on
  # an imputation object counts only positive weights.
  observed <- as.numeric(lavaan::nobs(fit_zeros))
  expect_true(is.finite(observed))
  expect_gte(observed, n_pos)
  expect_equal(observed, as.numeric(nrow(data)))

  skip_if_not_installed("mitools")
  frames <- mi_frames(2, seed = 90, zero = TRUE)
  imp <- mi_design(2, seed = 90, zero = TRUE)
  fit_imp <- svysem("y ~ x", imp, fixed.x = FALSE)
  expect_equal(nobs(fit_imp), sum(frames[[1]]$w > 0))
})

test_that("an unstratified cluster sample saturated mean model matches svymean", {
  data <- small_data(10)
  design <- svydesign(ids = ~psu, weights = ~w, data = data)
  fit <- svysem(
    "y ~ 1\nx1 ~ 1\nx2 ~ 1", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  means <- svymean(~ y + x1 + x2, design)
  expect_rel_equal(coef(fit)[c("y~1", "x1~1", "x2~1")], coef(means), 1e-6)
  expect_rel_equal(
    sqrt(diag(vcov(fit)))[c("y~1", "x1~1", "x2~1")],
    SE(means),
    1e-6
  )
})

test_that("a variance matches the mean of squared residuals, not svyvar", {
  design <- cluster_design(11)
  mu <- as.numeric(coef(svymean(~y, design)))
  variables <- data.frame(design$variables, check.names = FALSE, stringsAsFactors = FALSE)
  variables$z_sq <- (variables$y - mu)^2
  moment_design <- design
  moment_design$variables <- variables
  moment <- svymean(~z_sq, moment_design)
  fit <- svysem("y ~~ y", design, fixed.x = FALSE, meanstructure = TRUE)
  expect_rel_equal(pe_est(fit, "y", "~~", "y"), coef(moment), 1e-6)
  expect_rel_equal(pe_se(fit, "y", "~~", "y"), SE(moment), 1e-6)
  unbiased <- as.numeric(svyvar(~y, design))
  expect_gt(abs(pe_est(fit, "y", "~~", "y") - unbiased) / abs(unbiased), 1e-4)
})

test_that("fixed.x TRUE and FALSE both return a finite covariance for the slope", {
  design <- cluster_design(12)
  fit_true <- svysem("y ~ x1 + x2", design, fixed.x = TRUE, meanstructure = TRUE)
  fit_false <- svysem("y ~ x1 + x2", design, fixed.x = FALSE, meanstructure = TRUE)
  for (fit in list(fit_true, fit_false)) {
    V <- vcov(fit)
    expect_identical(dim(V), c(length(coef(fit)), length(coef(fit))))
    expect_true("y~x1" %in% names(coef(fit)))
    expect_true(all(is.finite(V)))
    expect_gt(sqrt(V["y~x1", "y~x1"]), 0)
  }
})

test_that("meanstructure FALSE adds no intercepts and TRUE matches svymean", {
  design <- cluster_design(13)
  fit_false <- svysem("y ~ x1", design, fixed.x = FALSE, meanstructure = FALSE)
  pe <- parameterEstimates(fit_false, se = TRUE)
  expect_equal(sum(pe$op == "~1"), 0L)
  expect_false(any(grepl("~1", names(coef(fit_false)), fixed = TRUE)))

  fit_true <- svysem(
    "y ~ 1\nx1 ~ 1", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  means <- svymean(~ y + x1, design)
  expect_rel_equal(coef(fit_true)[c("y~1", "x1~1")], coef(means)[c("y", "x1")], 1e-6)
  expect_rel_equal(
    sqrt(diag(vcov(fit_true)))[c("y~1", "x1~1")],
    SE(means)[c("y", "x1")],
    1e-6
  )
})

test_that("svycfa and svygrowth return design-based svysem fits", {
  design <- cluster_design(14)
  fits <- list(
    svycfa("f =~ y1 + y2 + y3", design),
    svygrowth(
      "i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4\n s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4",
      design
    )
  )
  for (fit in fits) {
    expect_s4_class(fit, "svysem")
    expect_true(inherits(fit, "lavaan"))
    expect_equal(fit@df, as.numeric(degf(design)))
    expect_identical(fit@estimator, "MLM")
    expect_identical(fit@variance_method, "linearization")
    expect_true(is.call(fit@svysem_call))
    V <- vcov(fit)
    expect_equal(V, t(V))
    expect_true(all(is.finite(V)))
    expect_true(is.finite(as.numeric(fitMeasures(fit, "chisq.scaled"))))
  }
})

test_that("parameterEstimates(se = TRUE) returns the design-based standard errors", {
  design <- cluster_design(15)
  fit <- svysem("y ~ x1 + x2", design, fixed.x = FALSE, meanstructure = TRUE)
  pe0 <- parameterEstimates(fit)
  ss0 <- standardizedSolution(fit)
  # lavaan 0.6-21 omits se unless se=TRUE when class[1] is not "lavaan".
  # svysem has to stay the first class so show() and summary() dispatch.
  expect_false("se" %in% names(pe0))
  expect_false("se" %in% names(ss0))

  pe <- parameterEstimates(fit, se = TRUE)
  names_pe <- paste0(pe$lhs, pe$op, pe$rhs)
  free <- names(coef(fit))
  got <- pe$se[match(free, names_pe)]
  want <- sqrt(diag(vcov(fit)))[free]
  expect_false(anyNA(got))
  expect_rel_equal(got, want, 1e-8)
  expect_true(all(is.finite(standardizedSolution(fit, se = TRUE)$se)))
})

test_that("lavTestWald is finite and a defined mean difference matches the delta method", {
  design <- cluster_design(16)
  labeled <- svysem(
    "y ~ b*x1 + x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  wald <- lavTestWald(labeled, "b == 0")
  expect_true(is.finite(wald$stat))

  fit <- svysem(
    "y ~ my*1\nx1 ~ mx*1\ndiff := my - mx", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  pe <- parameterEstimates(fit, se = TRUE)
  diff_se <- pe$se[pe$lhs == "diff" & pe$op == ":="]
  V <- vcov(fit)
  delta <- sqrt(V["my", "my"] + V["mx", "mx"] - 2 * V["my", "mx"])
  expect_rel_equal(diff_se, delta, 1e-6)
})

test_that("summary and show mention the estimator or the design degrees of freedom", {
  design <- cluster_design(17)
  fit <- svysem("y ~ x1", design, fixed.x = FALSE, meanstructure = TRUE)
  expect_output(show(fit), "MLM|Design df")
  expect_output(print(summary(fit)), "MLM|degrees of freedom")
})

test_that("dropping survey.design2 makes survey fail, so the design method is called directly", {
  design <- cluster_design(18)
  old <- design
  class(old) <- "survey.design"
  expect_identical(class(old), "survey.design")
  survey_err <- tryCatch(suppressWarnings(svymean(~y, old)), error = function(e) e)
  expect_s3_class(survey_err, "error")
  expect_identical(svysem.survey.design, svysem.survey.design2)
  fit <- suppressWarnings(svysem.survey.design(
    "y ~ x1", design,
    fixed.x = FALSE, meanstructure = TRUE
  ))
  direct <- svysem(
    "y ~ x1", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  expect_s4_class(fit, "svysem")
  expect_rel_equal(coef(fit), coef(direct), 1e-8)
})

test_that("replicate designs, non-designs, and reserved estimators are classed errors", {
  design <- cluster_design(19)
  replicate <- as.svrepdesign(design)
  err <- expect_classed(svysem("y ~ x1", replicate), "svysem_error_design")
  expect_match(conditionMessage(err), "this design is not supported yet", fixed = TRUE)

  expect_classed(
    svysem("y ~ x1", data.frame(y = 1, x1 = 1)),
    "svysem_error_design"
  )
  expect_classed(svysem("y ~ x1", design, estimator = "WLSMV"), "svysem_error_estimator")
  expect_classed(svysem("y ~ x1", design, estimator = "MLR"), "svysem_error_estimator")
  missing_err <- expect_classed(
    svysem("y ~ x1", design, missing = "ml"),
    "svysem_error_missing"
  )
  expect_match(conditionMessage(missing_err), "MLR", fixed = TRUE)
})

test_that("listwise missing rows stay in the design and the standard errors differ from dropping them", {
  # One incomplete row inside a multi-row PSU leaves the same cluster totals
  # as dropping that row. Removing every row of one PSU keeps the PSU, at
  # weight 0, in the design variance; rebuilding drops the PSU.
  base <- simulate_design(n_psu = 3, n_per = 4, n_strata = 4, nest = TRUE, seed = 15)
  variables <- data.frame(base$variables, check.names = FALSE, stringsAsFactors = FALSE)
  variables$y[variables$strat == 1 & variables$psu == 1] <- NA
  incomplete <- base
  incomplete$variables <- variables
  fit <- svysem("y ~ x1", incomplete, fixed.x = FALSE, meanstructure = TRUE)
  expect_true(all(is.finite(coef(fit))))
  kept <- variables[stats::complete.cases(variables[, c("y", "x1")]), , drop = FALSE]
  rebuilt <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = kept, nest = TRUE
  )
  fit_rebuilt <- svysem("y ~ x1", rebuilt, fixed.x = FALSE, meanstructure = TRUE)
  se_domain <- sqrt(diag(vcov(fit)))
  se_rebuilt <- sqrt(diag(vcov(fit_rebuilt)))
  relative <- abs(se_domain - se_rebuilt) / pmax(abs(se_rebuilt), 1e-12)
  expect_gt(max(relative), 1e-3)
})

test_that("a missing variable, a negative weight, and all-zero weights are classed errors", {
  design <- cluster_design(21)
  expect_classed(svysem("y ~ not_in_design", design), "svysem_error_variables")

  negative <- design
  negative$prob <- design$prob * 1
  negative$prob[1] <- -1
  expect_classed(svysem("y ~ x1", negative), "svysem_error_weights")

  zeros <- design
  zeros$prob <- rep(Inf, length(design$prob))
  expect_classed(svysem("y ~ x1", zeros), "svysem_error_weights")
})

test_that("unused columns are ignored when the model uses two of five numeric columns", {
  set.seed(22)
  n <- 40L
  full <- data.frame(
    y = stats::rnorm(n),
    x1 = stats::rnorm(n),
    z1 = stats::rnorm(n),
    z2 = stats::rnorm(n),
    z3 = stats::rnorm(n),
    w = stats::runif(n, 0.5, 2),
    strat = rep(1:4, each = 10),
    psu = rep(rep(1:2, each = 5), 4)
  )
  full$y <- full$y + 0.4 * full$x1
  reduced <- full[, c("y", "x1", "w", "strat", "psu")]
  fit_full <- svysem(
    "y ~ x1",
    svydesign(ids = ~psu, strata = ~strat, weights = ~w, data = full, nest = TRUE),
    fixed.x = FALSE, meanstructure = TRUE
  )
  fit_reduced <- svysem(
    "y ~ x1",
    svydesign(ids = ~psu, strata = ~strat, weights = ~w, data = reduced, nest = TRUE),
    fixed.x = FALSE, meanstructure = TRUE
  )
  expect_rel_equal(coef(fit_full), coef(fit_reduced), 1e-8)
  expect_rel_equal(sqrt(diag(vcov(fit_full))), sqrt(diag(vcov(fit_reduced))), 1e-8)
})

test_that("sampling.weights.normalization none and total give the same coefficients", {
  design <- cluster_design(23)
  fit_none <- svysem(
    "y ~ x1", design,
    fixed.x = FALSE, meanstructure = TRUE,
    sampling.weights.normalization = "none"
  )
  fit_total <- svysem(
    "y ~ x1", design,
    fixed.x = FALSE, meanstructure = TRUE,
    sampling.weights.normalization = "total"
  )
  expect_true(all(is.finite(coef(fit_none))))
  expect_true(all(is.finite(coef(fit_total))))
  expect_rel_equal(coef(fit_none), coef(fit_total), 1e-8)
})

test_that("group, cluster, and ordered arguments are classed data errors", {
  design <- cluster_design(24)
  expect_classed(svysem("y ~ x1", design, group = "strat"), "svysem_error_data")
  expect_classed(svysem("y ~ x1", design, cluster = "psu"), "svysem_error_data")
  expect_classed(svysem("y ~ x1", design, ordered = "y"), "svysem_error_data")
})

test_that("imputed designs pool sequential and parallel fits and record failures", {
  skip_if_not_installed("mitools")
  skip_if_not_installed("parallel")
  par_fn <- function(X, FUN) parallel::mclapply(X, FUN, mc.cores = 2)
  for (m in c(2L, 5L)) {
    design <- mi_design(m, seed = 250 + m)
    sequential <- svysem("y ~ x", design, fixed.x = FALSE)
    parallel <- svysem("y ~ x", design, fixed.x = FALSE, parallel = par_fn)
    expect_s4_class(sequential, "svysem_mi")
    expect_equal(length(sequential@fits), m)
    expect_equal(coef(parallel), coef(sequential), tolerance = 0)
    expect_equal(vcov(parallel), vcov(sequential), tolerance = 0)
  }

  bad <- mi_design(3, seed = 260, bad = 2L)
  omitted <- svysem("y ~ x", bad, fixed.x = FALSE, on_fail = "omit")
  expect_equal(length(omitted@fits), 2L)
  expect_equal(nrow(omitted@status), 3L)
  expect_true(omitted@status$no.conv[[2]])
  expect_true(nzchar(omitted@status$message[[2]]))
  expect_identical(omitted@status$message[[1]], "")
  expect_identical(omitted@status$message[[3]], "")
  manual <- lapply(bad$designs[-2L], function(design) {
    svysem("y ~ x", design, fixed.x = FALSE)
  })
  expect_rel_equal(
    coef(omitted),
    pooled_estimates(manual, df.com = omitted@df)$est,
    1e-8
  )

  stopped <- expect_classed(
    svysem("y ~ x", bad, fixed.x = FALSE, on_fail = "stop"),
    "svysem_error_imputation"
  )
  expect_match(conditionMessage(stopped), "Imputation 2", fixed = TRUE)
})

test_that("pooled_estimates on a real imputation fit follows Rubin's rules", {
  skip_if_not_installed("mitools")
  design <- mi_design(5, seed = 270)
  fit <- svysem("y ~ b1*x + b2*x2", design, fixed.x = FALSE, meanstructure = TRUE)
  pe <- pooled_estimates(fit)
  expect_identical(
    names(pe),
    c("est", "se", "t", "df", "p", "ci.lower", "ci.upper", "riv", "lambda", "fmi", "m")
  )
  expect_rel_equal(pe$se, sqrt(diag(vcov(fit)))[rownames(pe)], 1e-8)
  expect_rel_equal(pe$t, pe$est / pe$se, 1e-8)
  tcrit <- stats::qt(1 - (1 - 0.95) / 2, df = pe$df)
  expect_rel_equal(pe$ci.lower, pe$est - tcrit * pe$se, 1e-8)
  expect_rel_equal(pe$ci.upper, pe$est + tcrit * pe$se, 1e-8)

  over <- pooled_estimates(fit, df.com = 40)
  expect_rel_equal(over$est, pe$est, 1e-12)
  expect_false(isTRUE(all.equal(over$df, pe$df)))

  standardized <- pooled_estimates(fit, standardized = TRUE)
  expect_gt(nrow(standardized), 0L)
  expect_true(all(is.finite(standardized$est)))
  expect_true(all(is.finite(standardized$se)))
})

test_that("pooled_wald on two labeled regressions returns a finite F", {
  skip_if_not_installed("mitools")
  design <- mi_design(5, seed = 280)
  fit <- svysem("y ~ b1*x + b2*x2", design, fixed.x = FALSE, meanstructure = TRUE)
  one <- pooled_wald(fit, "b1 == b2")
  expect_true(is.finite(one$F))
  expect_equal(one$df1, 1)
  expect_true(one$pvalue > 0 && one$pvalue <= 1)

  two <- pooled_wald(fit, c("b1 == b2", "b1 == 0"))
  expect_equal(two$df1, 2)
  expect_true(is.finite(two$F))
  expect_true(is.finite(two$pvalue))
  expect_true(two$pvalue > 0 && two$pvalue <= 1)
})

test_that("fitMeasures on an imputation fit returns finite scaled and robust indices", {
  skip_if_not_installed("mitools")
  # Four indicators so the factor model has positive degrees of freedom.
  # Three indicators are just identified and the D2 rule rejects a zero df.
  design <- mi_design(2, seed = 290)
  fit <- svycfa("f =~ y1 + y2 + y3 + y4", design)
  measures <- fitMeasures(fit)
  expect_type(measures, "double")
  expect_false(is.null(names(measures)))
  needed <- c(
    "chisq.scaled",
    "cfi.robust.mean", "cfi.robust.min", "cfi.robust.max",
    "rmsea.robust.mean", "rmsea.robust.min", "rmsea.robust.max"
  )
  expect_true(all(needed %in% names(measures)))
  expect_true(all(is.finite(measures)))
})

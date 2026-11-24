library(lavaan)

svysem <- function(...) suppressWarnings(svysem::svysem(...))

# Finding 11. A block-diagonal NACOV drops the covariance between groups
# when primary sampling units contain both of them. Multigroup MLM therefore
# uses replicate weights. The slow test records that the block-diagonal
# contrast is the one that misses the jackknife.

spanning_design <- function(n_psu, n_per, seed, shock = 1.5) {
  set.seed(seed)
  n_strata <- 2L
  psu <- rep(rep(seq_len(n_psu), each = n_per), times = n_strata)
  strat <- rep(seq_len(n_strata), each = n_psu * n_per)
  key <- interaction(strat, psu, drop = TRUE)
  cluster <- stats::rnorm(nlevels(key))[key]
  g <- rep(rep(c(1L, 2L), each = n_per / 2L), times = n_strata * n_psu)
  x <- stats::rnorm(length(g)) + 0.3 * cluster
  y <- 0.3 * (g == 2L) + 0.6 * x + 0.15 * (g == 2L) * x +
    shock * cluster + stats::rnorm(length(g), 0, 0.4)
  frame <- data.frame(
    y = y,
    x = x,
    g = factor(g),
    w = stats::runif(length(g), 0.5, 2),
    psu = psu,
    strat = strat
  )
  survey::svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = frame, nest = TRUE
  )
}

contrast_model <- "y ~ c(b1, b2)*x\ndiff := b1 - b2"

contrast_se <- function(fit) {
  covariance <- vcov(fit)
  sqrt(
    covariance["b1", "b1"] + covariance["b2", "b2"] -
      2 * covariance["b1", "b2"]
  )
}

test_that("multigroup MLM chi-square is the maximum-likelihood statistic", {
  set.seed(4)
  n_psu <- 6L
  n_per <- 8L
  n_strata <- 2L
  psu <- rep(rep(seq_len(n_psu), each = n_per), times = n_strata)
  strat <- rep(seq_len(n_strata), each = n_psu * n_per)
  key <- interaction(strat, psu, drop = TRUE)
  cluster <- stats::rnorm(nlevels(key))[key]
  g <- rep(rep(c(1L, 2L), each = n_per / 2L), times = n_strata * n_psu)
  factor_score <- stats::rnorm(length(g)) + 0.4 * cluster
  frame <- data.frame(
    y1 = 0.8 * factor_score + stats::rnorm(length(g), 0, 0.5),
    y2 = 0.7 * factor_score + stats::rnorm(length(g), 0, 0.5),
    y3 = 0.6 * factor_score + stats::rnorm(length(g), 0, 0.5),
    y4 = 0.5 * factor_score + 0.2 * (g == 2L) + stats::rnorm(length(g), 0, 0.5),
    g = factor(g),
    w = stats::runif(length(g), 0.5, 2),
    psu = psu,
    strat = strat
  )
  design <- survey::svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = frame, nest = TRUE
  )
  model <- "f =~ y1 + y2 + y3 + y4"
  fit <- svycfa(model, design, group = "g")
  oracle <- lavaan::cfa(
    model, data = frame, group = "g", sampling.weights = "w",
    estimator = "ML", test = "standard", se = "none"
  )
  expect_equal(
    as.numeric(fitMeasures(fit, "chisq")),
    as.numeric(lavaan::fitMeasures(oracle, "chisq")),
    tolerance = 1e-8
  )
  expect_false("chisq.scaled" %in% names(fitMeasures(fit)))
  expect_false(any(vapply(fit@test, function(element) {
    identical(element$test, "yuan.bentler.mplus")
  }, logical(1))))
})

test_that("multigroup MLM returns replicate-weight standard errors", {
  design <- spanning_design(5L, 4L, seed = 3L)
  fit <- svysem(
    contrast_model, design,
    group = "g", fixed.x = FALSE, meanstructure = TRUE
  )
  expect_s4_class(fit, "svysem")
  expect_identical(fit@estimator, "MLM")
  expect_identical(fit@variance_method, "replication")
  expect_identical(lavInspect(fit, "ngroups"), 2L)
  expect_true(isTRUE(lavInspect(fit, "converged")))
  covariance <- vcov(fit)
  expect_true(isSymmetric(covariance))
  expect_true(all(is.finite(covariance)))
  expect_true(all(c("b1", "b2") %in% colnames(covariance)))

  estimates <- parameterEstimates(fit, se = TRUE)
  expect_equal(
    estimates$se[estimates$lhs == "diff"],
    unname(contrast_se(fit)),
    tolerance = 1e-8
  )
  wald <- lavTestWald(fit, "diff == 0")
  expect_true(is.finite(wald$stat))
  expect_true(is.finite(unname(fitMeasures(fit, "chisq"))))
  expect_error(
    svysem(
      contrast_model, design,
      group = "g", replicate.type = "bootstrap"
    ),
    class = "svysem_error_design"
  )
})

test_that("a between-group contrast matches the jackknife and not the block diagonal", {
  skip_on_cran()
  design <- spanning_design(30L, 8L, seed = 13L, shock = 0.4)
  fit <- svysem(
    contrast_model, design,
    group = "g", fixed.x = FALSE, meanstructure = TRUE
  )
  replicate <- survey::as.svrepdesign(design, type = "JKn", compress = FALSE)
  oracle <- survey::withReplicates(replicate, function(w, data) {
    data$w <- as.numeric(w)
    replicated <- lavaan::sem(
      contrast_model, data,
      group = "g", sampling.weights = "w", estimator = "ML", se = "none",
      fixed.x = FALSE, meanstructure = TRUE
    )
    # se = "none" leaves the := row out of coef(). b1 and b2 are the slopes.
    estimate <- lavaan::coef(replicated)
    unname(estimate[["b1"]] - estimate[["b2"]])
  })
  oracle_se <- sqrt(as.numeric(vcov(oracle)))
  expect_rel(contrast_se(fit), oracle_se, 1e-5)

  group_one <- svysem(
    "y ~ x", subset(design, g == "1"),
    fixed.x = FALSE, meanstructure = TRUE
  )
  group_two <- svysem(
    "y ~ x", subset(design, g == "2"),
    fixed.x = FALSE, meanstructure = TRUE
  )
  expect_equal(unname(coef(fit)["b1"]), unname(coef(group_one)["y~x"]), tolerance = 1e-6)
  expect_equal(unname(coef(fit)["b2"]), unname(coef(group_two)["y~x"]), tolerance = 1e-6)
  block <- sqrt(vcov(group_one)["y~x", "y~x"] + vcov(group_two)["y~x", "y~x"])
  block_gap <- abs(block - oracle_se) / oracle_se
  replicate_gap <- abs(contrast_se(fit) - oracle_se) / oracle_se
  expect_gt(block_gap, 0.05)
  expect_lt(replicate_gap, block_gap)
})

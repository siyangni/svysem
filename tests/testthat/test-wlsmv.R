library(lavaan)

svysem <- function(...) suppressWarnings(svysem::svysem(...))
svycfa <- function(...) suppressWarnings(svysem::svycfa(...))
svygrowth <- function(...) suppressWarnings(svysem::svygrowth(...))

cfa_model <- "f =~ y1 + y2 + y3 + y4"

weighted_lavaan <- function(data, model, ...) {
  args <- list(
    model = model, data = data, estimator = "WLSMV",
    sampling.weights = "w", ...
  )
  # lavaan 0.7 requires ordered= even when every indicator is continuous.
  if (is.null(args$ordered)) {
    args$ordered <- FALSE
  }
  suppressWarnings(do.call(lavaan::cfa, args))
}

clustered <- function(seed) {
  simulate_design(n_psu = 6, n_per = 8, n_strata = 2, nest = TRUE, seed = seed)
}

test_that("WLSMV point estimates match lavaan and clustering changes the variance", {
  design <- clustered(3)
  fit <- svycfa(cfa_model, design, estimator = "WLSMV")
  expect_s4_class(fit, "svysem")
  expect_identical(fit@estimator, "WLSMV")
  expect_identical(fit@variance_method, "replication")
  expect_true(isTRUE(lavInspect(fit, "converged")))

  lavaan_fit <- weighted_lavaan(design$variables, cfa_model)
  expect_rel(coef(fit), coef(lavaan_fit), 1e-6)

  independent <- design$variables
  independent$id <- seq_len(nrow(independent))
  independent_design <- survey::svydesign(
    ids = ~id, weights = ~w, data = independent
  )
  independent_fit <- svycfa(cfa_model, independent_design, estimator = "WLSMV")
  expect_rel(coef(independent_fit), coef(lavaan_fit), 1e-6)

  clustered_se <- sqrt(diag(vcov(fit)))
  independent_se <- sqrt(diag(vcov(independent_fit)))
  expect_rel_gap <- max(abs(clustered_se - independent_se) / pmax(abs(independent_se), 1e-8))
  expect_gt(expect_rel_gap, 0.02)
  clustered_test <- fitMeasures(fit, "chisq.scaled")
  independent_test <- fitMeasures(independent_fit, "chisq.scaled")
  expect_gt(abs(unname(clustered_test) - unname(independent_test)), 1e-3)
  expect_true(all(is.finite(clustered_se)))
})

test_that("listwise deletion keeps design-based standard errors", {
  design <- clustered(4)
  design$variables$y1[1:4] <- NA
  design$variables$y3[5:7] <- NA
  listwise <- svycfa(
    cfa_model, design, estimator = "WLSMV", missing = "listwise", se = "none"
  )
  expect_rel(
    coef(listwise),
    coef(weighted_lavaan(design$variables, cfa_model, missing = "listwise")),
    1e-6
  )
  expect_true(all(is.finite(sqrt(diag(vcov(listwise))))))
  expect_lt(as.numeric(lavInspect(listwise, "nobs")), nrow(design$variables))
})

test_that("ordered listwise and pairwise WLSMV match lavaan", {
  # Continuous pairwise DWLS is rejected by lavaan 0.6-21. Pairwise is the
  # ordered-data estimator.
  design <- clustered(4)
  cut3 <- function(x) ordered(findInterval(x, c(-0.4, 0.4)))
  for (column in c("y1", "y2", "y3", "y4")) {
    design$variables[[column]] <- cut3(design$variables[[column]])
  }
  design$variables$y1[1:4] <- NA
  design$variables$y3[5:7] <- NA
  listwise <- svycfa(cfa_model, design, estimator = "WLSMV", missing = "listwise")
  pairwise <- svycfa(cfa_model, design, estimator = "WLSMV", missing = "pairwise")
  expect_true(isTRUE(lavInspect(listwise, "converged")))
  expect_true(isTRUE(lavInspect(pairwise, "converged")))
  expect_rel(
    coef(listwise),
    coef(weighted_lavaan(design$variables, cfa_model, missing = "listwise")),
    1e-5
  )
  expect_rel(
    coef(pairwise),
    coef(weighted_lavaan(design$variables, cfa_model, missing = "pairwise")),
    1e-5
  )
  expect_gt(max(abs(unname(coef(listwise)) - unname(coef(pairwise)))), 1e-3)
})

test_that("a grouped WLSMV fit matches lavaan's coefficients and drops the cross-group block", {
  design <- clustered(5)
  design$variables$g <- ifelse(design$variables$strat == 1, "a", "b")
  fit <- svycfa(cfa_model, design, estimator = "WLSMV", group = "g")
  lavaan_fit <- weighted_lavaan(design$variables, cfa_model, group = "g")
  expect_true(isTRUE(lavInspect(fit, "converged")))
  expect_rel(coef(fit), coef(lavaan_fit), 1e-5)
  expect_identical(lavInspect(fit, "ngroups"), 2L)
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))

  restricted <- svycfa(
    "f =~ y1 + y2 + y3 + y4\ny2 ~~ a*y2\ny3 ~~ a*y3",
    design, estimator = "WLSMV"
  )
  comparison <- lavTestLRT(
    svycfa(cfa_model, design, estimator = "WLSMV"),
    restricted
  )
  expect_s3_class(comparison, "anova")
  expect_true(any(is.finite(unlist(comparison[2L, vapply(comparison, is.numeric, logical(1))]))))
})

test_that("configural, metric, and scalar WLSMV models are nested", {
  design <- clustered(6)
  design$variables$g <- ifelse(design$variables$strat == 1, "a", "b")
  intercepts <- "y1 ~ 1\ny2 ~ 1\ny3 ~ 1\ny4 ~ 1"
  configural <- paste("f =~ y1 + y2 + y3 + y4", intercepts, sep = "\n")
  metric <- paste("f =~ y1 + l2*y2 + l3*y3 + l4*y4", intercepts, sep = "\n")
  scalar <- paste(
    "f =~ y1 + l2*y2 + l3*y3 + l4*y4",
    "y1 ~ i1*1\ny2 ~ i2*1\ny3 ~ i3*1\ny4 ~ i4*1",
    sep = "\n"
  )
  fits <- lapply(list(configural, metric, scalar), function(model) {
    svycfa(
      model, design, estimator = "WLSMV", group = "g", meanstructure = TRUE
    )
  })
  expect_true(all(vapply(fits, function(fit) isTRUE(lavInspect(fit, "converged")), logical(1))))
  expect_true(all(vapply(fits, function(fit) lavInspect(fit, "ngroups") == 2L, logical(1))))
  comparison <- lavTestLRT(fits[[1L]], fits[[2L]], fits[[3L]])
  expect_s3_class(comparison, "anova")
  numeric_columns <- comparison[vapply(comparison, is.numeric, logical(1))]
  expect_true(any(is.finite(unlist(numeric_columns[-1L, , drop = FALSE]))))
})

test_that("replicate controls are checked and an existing replicate design is reused", {
  design <- clustered(6)
  expect_error(
    svycfa(cfa_model, design, estimator = "WLSMV", replicate.type = "JK1"),
    class = "svysem_error_design"
  )
  expect_error(
    svycfa(cfa_model, design, estimator = "WLSMV", replicates = 8),
    class = "svysem_error_design"
  )
  expect_error(
    svycfa(cfa_model, design, estimator = "WLSMV", replicate.type = "bootstrap"),
    class = "svysem_error_design"
  )
  expect_error(
    svycfa(cfa_model, design, estimator = "WLSMV", missing = "ml"),
    class = "svysem_error_missing"
  )
  set.seed(1)
  bootstrap <- svycfa(
    cfa_model, design, estimator = "WLSMV",
    replicate.type = "bootstrap", replicates = 6
  )
  expect_identical(bootstrap@variance_method, "replication")
  expect_true(all(is.finite(coef(bootstrap))))

  replicate <- survey::as.svrepdesign(design, type = "JKn", compress = FALSE)
  from_replicate <- svycfa(cfa_model, replicate, estimator = "WLSMV")
  expect_identical(from_replicate@estimator, "WLSMV")
  expect_identical(from_replicate@variance_method, "replication")
  expect_error(
    svycfa(cfa_model, replicate, estimator = "WLSMV", replicate.type = "JKn"),
    class = "svysem_error_design"
  )

  seen <- 0L
  fit <- svycfa(
    cfa_model, design, estimator = "WLSMV",
    parallel = function(X, FUN) {
      seen <<- length(X)
      lapply(X, FUN)
    }
  )
  expect_gte(seen, 2L)
  expect_s4_class(fit, "svysem")
  expect_error(
    svycfa(
      cfa_model, design, estimator = "WLSMV",
      parallel = function(X, FUN) lapply(X, FUN)[0]
    ),
    class = "svysem_error_design"
  )
})

test_that("ordinal WLSMV with one case per PSU matches lavaan's standard errors", {
  skip_on_cran()
  # Equal cutpoints leave the factor variance and one loading about 4 percent
  # away from lavaan's analytic gamma. Uneven cuts, including one binary
  # indicator, keep every standard error inside 1 percent on this seed.
  set.seed(4)
  n <- 300L
  eta <- stats::rnorm(n)
  noise <- function() stats::rnorm(n, 0, 0.4)
  data <- data.frame(
    y1 = ordered(findInterval(eta + noise(), c(-0.5, 0.5))),
    y2 = ordered(findInterval(0.9 * eta + noise(), c(0))),
    y3 = ordered(findInterval(0.8 * eta + noise(), c(-0.4, 0.6))),
    w = stats::runif(n, 0.5, 2),
    id = seq_len(n)
  )
  expect_identical(
    vapply(data[c("y1", "y2", "y3")], function(column) length(unique(column)), integer(1)),
    c(y1 = 3L, y2 = 2L, y3 = 3L)
  )
  design <- survey::svydesign(ids = ~id, weights = ~w, data = data)
  model <- "f =~ y1 + y2 + y3"
  fit <- svycfa(model, design, estimator = "WLSMV")
  lavaan_fit <- lavaan::cfa(
    model, data = data, estimator = "WLSMV", sampling.weights = "w", ordered = c("y1", "y2", "y3")
  )
  expect_true(isTRUE(lavInspect(fit, "converged")))
  expect_true(isTRUE(lavInspect(lavaan_fit, "converged")))
  expect_rel(coef(fit), coef(lavaan_fit), 1e-5)
  expect_rel(sqrt(diag(vcov(fit))), sqrt(diag(vcov(lavaan_fit))), 0.02)
})

test_that("mixed ordered and continuous WLSMV matches lavaan's coefficients", {
  design <- clustered(8)
  design$variables$y1 <- ordered(findInterval(design$variables$y1, c(-0.2, 0.4)))
  design$variables$y2 <- ordered(findInterval(design$variables$y2, 0))
  model <- "f =~ y1 + y2 + y3 + y4"
  fit <- svycfa(model, design, estimator = "WLSMV", ordered = c("y1", "y2"))
  lavaan_fit <- weighted_lavaan(
    design$variables, model, ordered = c("y1", "y2")
  )
  expect_true(isTRUE(lavInspect(fit, "converged")))
  expect_true(isTRUE(lavInspect(lavaan_fit, "converged")))
  expect_rel(coef(fit), coef(lavaan_fit), 1e-5)
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
})

test_that("measEq.syntax invariance models fit through svycfa", {
  skip_if_not_installed("semTools")
  design <- clustered(7)
  design$variables$g <- rep(c("a", "b"), length.out = nrow(design$variables))
  syntax_for <- function(constraints) {
    semTools::measEq.syntax(
      configural.model = cfa_model,
      data = design$variables,
      group = "g",
      group.equal = constraints,
      return.fit = FALSE
    )
  }
  fits <- lapply(
    list("", "loadings", c("loadings", "intercepts")),
    function(constraints) {
      svycfa(
        as.character(syntax_for(constraints)),
        design,
        estimator = "WLSMV",
        group = "g"
      )
    }
  )
  expect_true(all(vapply(fits, function(fit) {
    isTRUE(lavInspect(fit, "converged")) && lavInspect(fit, "ngroups") == 2L
  }, logical(1))))
  comparison <- lavTestLRT(fits[[1L]], fits[[2L]], fits[[3L]])
  expect_s3_class(comparison, "anova")
  numeric_columns <- comparison[vapply(comparison, is.numeric, logical(1))]
  expect_true(any(is.finite(unlist(numeric_columns[-1L, , drop = FALSE]))))
})

library(survey)
library(lavaan)

expect_regression <- function(design) {
  fit <- svysem(
    "y ~ x1 + x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  got <- regression_parts(fit)
  oracle <- glm_parts(design)
  expect_equal(got$est, oracle$est, tolerance = 1e-6)
  expect_equal(got$se, oracle$se, tolerance = 1e-6)
  fit
}

test_that("an unstratified regression matches svyglm, including the intercept", {
  design <- simulate_design(n_psu = 20, n_per = 6, n_strata = 0, seed = 1)
  expect_regression(design)
})

test_that("a stratified regression matches svyglm", {
  design <- simulate_design(n_psu = 8, n_per = 5, n_strata = 4, seed = 2)
  expect_regression(design)
})

test_that("nest = TRUE with reused primary sampling unit labels matches svyglm", {
  design <- simulate_design(
    n_psu = 8, n_per = 5, n_strata = 4, nest = TRUE, seed = 3
  )
  expect_regression(design)
})

test_that("a saturated mean model matches svymean", {
  design <- simulate_design(n_psu = 8, n_per = 5, n_strata = 4, seed = 4)
  fit <- svysem(
    "y ~ 1\nx1 ~ 1\nx2 ~ 1", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  means <- svymean(~ y + x1 + x2, design)
  est <- c(
    pe_est(fit, "y", "~1", ""),
    pe_est(fit, "x1", "~1", ""),
    pe_est(fit, "x2", "~1", "")
  )
  se <- c(
    pe_se(fit, "y", "~1", ""),
    pe_se(fit, "x1", "~1", ""),
    pe_se(fit, "x2", "~1", "")
  )
  expect_equal(est, unname(coef(means)[c("y", "x1", "x2")]), tolerance = 1e-6)
  expect_equal(se, unname(SE(means)[c("y", "x1", "x2")]), tolerance = 1e-6)
})

test_that("a variance matches the mean of squared residuals, not svyvar", {
  design <- simulate_design(n_psu = 8, n_per = 5, n_strata = 4, seed = 5)
  mu <- as.numeric(coef(svymean(~y, design)))
  variables <- data.frame(design$variables, check.names = FALSE, stringsAsFactors = FALSE)
  variables$z_sq <- (variables$y - mu)^2
  z_design <- design
  z_design$variables <- variables
  moment <- svymean(~z_sq, z_design)
  fit <- svysem("y ~~ y", design, fixed.x = FALSE, meanstructure = TRUE)
  expect_equal(pe_est(fit, "y", "~~", "y"), as.numeric(coef(moment)), tolerance = 1e-6)
  expect_equal(pe_se(fit, "y", "~~", "y"), as.numeric(SE(moment)), tolerance = 1e-6)
  unbiased <- as.numeric(svyvar(~y, design))
  expect_gt(abs(pe_est(fit, "y", "~~", "y") - unbiased) / unbiased, 1e-4)
})

test_that("a domain via subset() does not match a design rebuilt on the kept rows", {
  base <- simulate_design(
    n_psu = 8, n_per = 5, n_strata = 4, nest = TRUE, seed = 6
  )
  design <- subset(base, psu != 1)
  kept <- base$variables[base$variables$psu != 1, , drop = FALSE]
  rebuilt <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = kept, nest = TRUE
  )
  fit_domain <- svysem(
    "y ~ x1 + x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  fit_rebuilt <- svysem(
    "y ~ x1 + x2", rebuilt,
    fixed.x = FALSE, meanstructure = TRUE
  )
  se_domain <- regression_parts(fit_domain)$se
  se_rebuilt <- regression_parts(fit_rebuilt)$se
  relative <- abs(se_domain - se_rebuilt) / abs(se_rebuilt)
  expect_gt(max(relative), 1e-3)
})

test_that("lonely PSU options are read from survey and are not set by the fit", {
  design <- simulate_design(n_psu = 4, n_per = 5, lonely = TRUE, seed = 7)
  previous <- options(survey.lonely.psu = getOption("survey.lonely.psu"))
  on.exit(options(previous), add = TRUE)
  for (setting in c("adjust", "average", "remove", "certainty")) {
    options(survey.lonely.psu = setting)
    means <- tryCatch(svymean(~y, design), error = function(e) e)
    if (inherits(means, "error")) {
      next
    }
    fit <- svysem("y ~ 1", design, fixed.x = FALSE, meanstructure = TRUE)
    expect_equal(pe_est(fit, "y", "~1", ""), as.numeric(coef(means)), tolerance = 1e-6)
    expect_equal(pe_se(fit, "y", "~1", ""), as.numeric(SE(means)), tolerance = 1e-6)
    expect_identical(getOption("survey.lonely.psu"), setting)
  }
})

test_that("fixed.x = TRUE returns finite standard errors and lavaan's gamma dimension", {
  design <- simulate_design(n_psu = 8, n_per = 5, n_strata = 3, seed = 8)
  fit <- svysem("y ~ x1 + x2", design)
  gamma <- lavInspect(fit, "gamma")
  expect_equal(nrow(gamma), 6L)
  expect_equal(ncol(gamma), 6L)
  expect_true(all(is.finite(parameterEstimates(fit, se = TRUE)$se)))
  expect_true(all(is.finite(vcov(fit))))
})

test_that("listwise missing rows stay in the design as a zero-weight domain", {
  design <- simulate_design(n_psu = 8, n_per = 5, n_strata = 3, seed = 9)
  original <- as.numeric(weights(design, "sampling"))
  variables <- data.frame(design$variables, check.names = FALSE, stringsAsFactors = FALSE)
  variables$y[c(2, 9)] <- NA
  variables$x1[14] <- NA
  incomplete <- design
  incomplete$variables <- variables
  fit <- svysem(
    "y ~ x1 + x2", incomplete,
    fixed.x = FALSE, meanstructure = TRUE
  )
  expect_equal(as.numeric(weights(design, "sampling")), original)
  expect_equal(as.numeric(weights(incomplete, "sampling")), original)

  na_rows <- !stats::complete.cases(variables[, c("y", "x1", "x2")])
  zeroed <- incomplete
  zeroed$prob[na_rows] <- Inf
  filled <- data.frame(variables, check.names = FALSE, stringsAsFactors = FALSE)
  filled$y[na_rows] <- 0
  filled$x1[na_rows] <- 0
  filled$x2[na_rows] <- 0
  zeroed$variables <- filled
  got <- regression_parts(fit)
  oracle <- suppressWarnings(glm_parts(zeroed))
  expect_equal(got$est, oracle$est, tolerance = 1e-6)
  expect_equal(got$se, oracle$se, tolerance = 1e-6)
})

test_that("lavaan generics see the design-based fit", {
  design <- simulate_design(n_psu = 8, n_per = 5, n_strata = 2, seed = 10)
  fit <- svysem(
    "y ~ x1 + x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  labeled <- svysem(
    "y ~ b*x1 + x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  expect_s4_class(fit, "svysem")
  expect_true(is(fit, "lavaan"))
  expect_identical(fit@estimator, "MLM")
  expect_identical(fit@variance_method, "linearization")
  expect_equal(fit@df, as.numeric(degf(design)))
  expect_true(is.call(fit@svysem_call))
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(vcov(fit))))
  expect_s3_class(parameterEstimates(fit), "data.frame")
  expect_s3_class(standardizedSolution(fit), "data.frame")
  measures <- fitMeasures(fit, c("chisq", "df"))
  expect_true(all(is.finite(measures)))
  wald <- lavTestWald(labeled, "b == 0")
  expect_true(is.finite(wald$stat))
  expect_output(show(fit), "Estimator: MLM")
  expect_output(show(fit), "Design df:")
  expect_output(print(summary(fit)), "Design degrees of freedom:")
  expect_output(print(summary(fit)), "Estimator: MLM")

  cfa_fit <- svycfa("f =~ y1 + y2 + y3 + y4", design)
  growth_fit <- svygrowth(
    "i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4\n s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4",
    design
  )
  expect_s4_class(cfa_fit, "svysem")
  expect_s4_class(growth_fit, "svysem")
  expect_true(all(is.finite(coef(cfa_fit))))
  expect_true(all(is.finite(coef(growth_fit))))

  delegated <- svysem.survey.design(
    "y ~ x1 + x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  expect_equal(unname(coef(delegated)), unname(coef(fit)), tolerance = 1e-8)
})

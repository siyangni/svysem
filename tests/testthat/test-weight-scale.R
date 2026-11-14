library(survey)
library(lavaan)

test_that("multiplying every weight by a constant leaves means, coefficients, SEs and the scaled chi-square unchanged", {
  design <- simulate_design(
    n_psu = 8, n_per = 5, n_strata = 4, nest = TRUE, seed = 21
  )
  scaled_data <- data.frame(design$variables, check.names = FALSE, stringsAsFactors = FALSE)
  scaled_data$w <- scaled_data$w * 10
  scaled <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = scaled_data, nest = TRUE
  )

  fit <- svysem(
    "y ~ x1 + x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  fit_scaled <- svysem(
    "y ~ x1 + x2", scaled,
    fixed.x = FALSE, meanstructure = TRUE
  )
  expect_equal(regression_parts(fit)$est, regression_parts(fit_scaled)$est, tolerance = 1e-10)
  expect_equal(regression_parts(fit)$se, regression_parts(fit_scaled)$se, tolerance = 1e-10)

  means <- svysem(
    "y ~ 1\nx1 ~ 1\nx2 ~ 1", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
  means_scaled <- svysem(
    "y ~ 1\nx1 ~ 1\nx2 ~ 1", scaled,
    fixed.x = FALSE, meanstructure = TRUE
  )
  mean_est <- function(object) {
    c(
      pe_est(object, "y", "~1", ""),
      pe_est(object, "x1", "~1", ""),
      pe_est(object, "x2", "~1", "")
    )
  }
  mean_se <- function(object) {
    c(
      pe_se(object, "y", "~1", ""),
      pe_se(object, "x1", "~1", ""),
      pe_se(object, "x2", "~1", "")
    )
  }
  expect_equal(mean_est(means), mean_est(means_scaled), tolerance = 1e-10)
  expect_equal(mean_se(means), mean_se(means_scaled), tolerance = 1e-10)

  cfa_model <- "f =~ y1 + y2 + y3 + y4"
  cfa_fit <- svycfa(cfa_model, design)
  cfa_scaled <- svycfa(cfa_model, scaled)
  expect_equal(unname(coef(cfa_fit)), unname(coef(cfa_scaled)), tolerance = 1e-10)
  chisq <- fitMeasures(cfa_fit, "chisq.scaled")
  chisq_scaled <- fitMeasures(cfa_scaled, "chisq.scaled")
  expect_gt(as.numeric(chisq), 0)
  expect_equal(as.numeric(chisq), as.numeric(chisq_scaled), tolerance = 1e-10)
  se <- parameterEstimates(cfa_fit, se = TRUE)$se
  se_scaled <- parameterEstimates(cfa_scaled, se = TRUE)$se
  expect_equal(se, se_scaled, tolerance = 1e-10)
})

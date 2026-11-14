library(survey)
library(lavaan)

# lavaan.survey standard errors use svyvar products, which carry n / (n - 1),
# and a block-diagonal Gamma. ~~ follows that factor; ~1 does not.
se_ratio_ok <- function(oracle_se, fit_se, factor) {
  free <- is.finite(oracle_se) & is.finite(fit_se) & oracle_se > 0 & fit_se > 0
  if (!any(free)) {
    return(TRUE)
  }
  ratio <- oracle_se[free] / fit_se[free]
  all(abs(ratio / factor - 1) < 1e-6)
}

test_that("MLM matches the lavaan.survey fixtures up to the svyvar factor", {
  # lavaan warns about the defined parameter dmean. The fit is unchanged.
  withCallingHandlers(
    warning = function(w) {
      if (grepl("dmean", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    },
    expr = {
  paths <- list.files(
    testthat::test_path("fixtures"),
    pattern = "\\.rds$",
    full.names = TRUE
  )
  expect_length(paths, 3L)
  for (path in paths) {
    oracle <- readRDS(path)
    design <- svydesign(
      ids = stats::reformulate(oracle$ids),
      strata = stats::reformulate(oracle$strata),
      weights = stats::reformulate(oracle$weights),
      data = oracle$data,
      nest = oracle$nest
    )
    args <- oracle$lavaan_args
    # A just-identified covariance can estimate a negative variance. The
    # comparison does not require that variance to be positive.
    fit <- suppressWarnings(if (oracle$fitter == "growth") {
      do.call(svygrowth, c(list(model = oracle$model, design = design), args))
    } else {
      do.call(svysem, c(list(model = oracle$model, design = design), args))
    })
    expect_lt(max(abs(unname(coef(fit)) - unname(oracle$coef))), 1e-8)

    chisq <- as.numeric(fitMeasures(fit, "chisq.scaled"))
    expect_lt(abs(chisq - oracle$chisq_scaled) / abs(oracle$chisq_scaled), 1e-8)

    pe <- parameterEstimates(fit, se = TRUE)
    merged <- merge(
      oracle$parameters,
      pe[, c("lhs", "op", "rhs", "est", "se")],
      by = c("lhs", "op", "rhs"),
      suffixes = c(".oracle", ".fit"),
      sort = FALSE
    )
    expect_equal(nrow(merged), nrow(oracle$parameters))
    nobs <- sum(weights(design, "sampling") != 0)
    factor <- nobs / (nobs - 1)

    covariance <- merged$op == "~~"
    means <- merged$op == "~1"
    expect_true(se_ratio_ok(merged$se.oracle[covariance], merged$se.fit[covariance], factor))
    expect_true(se_ratio_ok(merged$se.oracle[means], merged$se.fit[means], 1))

    defined <- merged$op == ":="
    expect_true(any(defined))
    for (kind in c("mean", "variance")) {
      rows <- defined & merged$kind == kind
      expect_true(any(rows))
      kind_factor <- if (kind == "mean") 1 else factor
      expect_true(se_ratio_ok(merged$se.oracle[rows], merged$se.fit[rows], kind_factor))
    }

    std <- standardizedSolution(fit, se = TRUE)
    std_merged <- merge(
      oracle$standardized,
      std[, c("lhs", "op", "rhs", "se")],
      by = c("lhs", "op", "rhs"),
      suffixes = c(".oracle", ".fit"),
      sort = FALSE
    )
    std_covariance <- std_merged$op == "~~"
    expect_true(se_ratio_ok(
      std_merged$se.oracle[std_covariance],
      std_merged$se.fit[std_covariance],
      factor
    ))
  }
    }
  )
})

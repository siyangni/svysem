# Test 10. Pooling against mice, mitml and lavaan.mi.

pool_fixture <- function(seed = 20261113, df.com = 80) {
  set.seed(seed)
  p <- 3
  m <- 8
  means <- c(0.4, -0.2, 0.7)
  Qhat <- vapply(seq_len(m), function(i) {
    means + stats::rnorm(p, sd = c(0.05, 0.08, 0.03))
  }, numeric(p))
  rownames(Qhat) <- c("y~x", "y~z", "y~~y")
  Uhat <- array(NA_real_, c(p, p, m), dimnames = list(rownames(Qhat), rownames(Qhat), NULL))
  for (i in seq_len(m)) {
    a <- matrix(stats::rnorm(p * p, sd = 0.3), p, p)
    Uhat[, , i] <- crossprod(a) + diag(0.02, p)
  }
  list(Qhat = Qhat, Uhat = Uhat, df.com = df.com, fits = fits_from_qu(Qhat, Uhat))
}

test_that("pooled estimates match mice::pool.scalar and mitml::testEstimates", {
  skip_if_not_installed("mice")
  skip_if_not_installed("mitml")
  fx <- pool_fixture()
  got <- pooled_estimates(fx$fits, df.com = fx$df.com)
  mice_est <- mice_se <- mice_df <- mice_r <- mice_fmi <- numeric(nrow(fx$Qhat))
  for (j in seq_len(nrow(fx$Qhat))) {
    ps <- mice::pool.scalar(
      fx$Qhat[j, ], fx$Uhat[j, j, ],
      n = fx$df.com + 1, k = 1
    )
    mice_est[[j]] <- ps$qbar
    mice_se[[j]] <- sqrt(ps$t)
    mice_df[[j]] <- ps$df
    mice_r[[j]] <- ps$r
    mice_fmi[[j]] <- ps$fmi
  }
  expect_equal(got$est, mice_est, tolerance = 1e-10)
  expect_equal(got$se, mice_se, tolerance = 1e-10)
  expect_equal(got$df, mice_df, tolerance = 1e-10)
  expect_equal(got$riv, mice_r, tolerance = 1e-10)
  expect_equal(got$fmi, mice_fmi, tolerance = 1e-10)

  te <- mitml::testEstimates(qhat = fx$Qhat, uhat = fx$Uhat, df.com = fx$df.com)
  oracle <- te$estimates
  expect_equal(unname(got$est), unname(oracle[, "Estimate"]), tolerance = 1e-10)
  expect_equal(unname(got$se), unname(oracle[, "Std.Error"]), tolerance = 1e-10)
  expect_equal(unname(got$df), unname(oracle[, "df"]), tolerance = 1e-10)
  expect_equal(unname(got$riv), unname(oracle[, "RIV"]), tolerance = 1e-10)
  expect_equal(unname(got$fmi), unname(oracle[, "FMI"]), tolerance = 1e-10)
  expect_equal(got$t, got$est / got$se, tolerance = 1e-10)
  expect_equal(got$lambda, got$riv / (1 + got$riv), tolerance = 1e-10)
  tcrit <- stats::qt(1 - (1 - 0.95) / 2, df = got$df)
  expect_equal(got$ci.lower, got$est - tcrit * got$se, tolerance = 1e-10)
  expect_equal(got$ci.upper, got$est + tcrit * got$se, tolerance = 1e-10)
  expect_identical(
    names(got),
    c("est", "se", "t", "df", "p", "ci.lower", "ci.upper", "riv", "lambda", "fmi", "m")
  )
  expect_identical(rownames(got), rownames(fx$Qhat))

  abs_mice <- max(
    abs(got$est - mice_est), abs(got$se - mice_se),
    abs(got$df - mice_df), abs(got$riv - mice_r), abs(got$fmi - mice_fmi)
  )
  abs_mitml <- max(
    abs(unname(got$est) - unname(oracle[, "Estimate"])),
    abs(unname(got$se) - unname(oracle[, "Std.Error"])),
    abs(unname(got$df) - unname(oracle[, "df"])),
    abs(unname(got$riv) - unname(oracle[, "RIV"])),
    abs(unname(got$fmi) - unname(oracle[, "FMI"]))
  )
  expect_lte(abs_mice, 1e-10)
  expect_lte(abs_mitml, 1e-10)
})

test_that("riv of zero uses the Barnard-Rubin limit, matching mice", {
  skip_if_not_installed("mice")
  Qhat <- matrix(1.5, nrow = 1, ncol = 5, dimnames = list("a", NULL))
  Uhat <- array(0.2, dim = c(1, 1, 5))
  got <- pooled_estimates(fits_from_qu(Qhat, Uhat), df.com = 80)
  ps <- mice::pool.scalar(Qhat[1, ], Uhat[1, 1, ], n = 81, k = 1)
  expect_equal(got$riv, 0, tolerance = 1e-12)
  expect_equal(got$lambda, 0, tolerance = 1e-12)
  expect_equal(got$se, sqrt(ps$t), tolerance = 1e-10)
  expect_equal(got$df, ps$df, tolerance = 1e-10)
  expect_equal(got$fmi, ps$fmi, tolerance = 1e-10)
  got_inf <- pooled_estimates(fits_from_qu(Qhat, Uhat), df.com = NULL)
  expect_equal(got_inf$df, Inf)
  expect_equal(got_inf$fmi, 0)
})

test_that("pooled_wald matches mitml::testConstraints on lhs-op-rhs names", {
  skip_if_not_installed("mitml")
  fx <- pool_fixture()
  one <- "`y~x` - `y~z`"
  two <- c("`y~x` - `y~z`", "`y~x` + `y~z`")
  for (constraints in list(one, two)) {
    got <- pooled_wald(fx$fits, constraints = constraints, df.com = fx$df.com)
    oracle <- mitml::testConstraints(
      qhat = fx$Qhat, uhat = fx$Uhat,
      constraints = constraints, method = "D1", df.com = fx$df.com
    )$test
    expect_equal(unname(got$F), unname(oracle[, "F.value"]), tolerance = 1e-8)
    expect_equal(unname(got$df1), unname(oracle[, "df1"]), tolerance = 1e-8)
    expect_equal(unname(got$df2), unname(oracle[, "df2"]), tolerance = 1e-8)
    expect_equal(unname(got$pvalue), unname(oracle[, "P(>F)"]), tolerance = 1e-8)
    abs_wald <- max(abs(c(
      unname(got$F) - unname(oracle[, "F.value"]),
      unname(got$df1) - unname(oracle[, "df1"]),
      unname(got$df2) - unname(oracle[, "df2"]),
      unname(got$pvalue) - unname(oracle[, "P(>F)"])
    )))
    expect_lte(abs_wald, 1e-8)
  }
  expect_error(
    pooled_wald(fx$fits, "y~x > 0", df.com = fx$df.com),
    class = "svysem_error_constraint"
  )
  expect_error(
    pooled_wald(fx$fits, "not_a_parameter - 1", df.com = fx$df.com),
    class = "svysem_error_constraint"
  )
})

test_that("label and := constraints match mitml on rewritten parameter names", {
  skip_if_not_installed("mitml")
  set.seed(11)
  model <- "
    y1 ~ b2h * x
    y2 ~ b2c * x
    diff := b2h - b2c
    ratio := b2h / b2c
  "
  fits <- lapply(seq_len(6), function(i) {
    n <- 160
    x <- stats::rnorm(n)
    y1 <- (1 + 0.05 * i) * x + stats::rnorm(n)
    y2 <- (0.6 + 0.04 * i) * x + stats::rnorm(n)
    lavaan::sem(model, data = data.frame(y1 = y1, y2 = y2, x = x))
  })
  label <- pooled_wald(fits, "b2h == b2c", df.com = 80)
  defined <- pooled_wald(fits, "diff == 0", df.com = 80)
  rewritten <- "`y1~x` - `y2~x`"
  oracle <- mitml::testConstraints(
    fits, constraints = rewritten, method = "D1", df.com = 80
  )$test
  expect_equal(unname(label$F), unname(oracle[, "F.value"]), tolerance = 1e-8)
  expect_equal(unname(label$df1), unname(oracle[, "df1"]), tolerance = 1e-8)
  expect_equal(unname(label$df2), unname(oracle[, "df2"]), tolerance = 1e-8)
  expect_equal(unname(label$pvalue), unname(oracle[, "P(>F)"]), tolerance = 1e-8)
  expect_equal(defined$F, label$F, tolerance = 1e-8)
  expect_equal(defined$df2, label$df2, tolerance = 1e-8)
  backtick <- pooled_wald(fits, rewritten, df.com = 80)
  expect_equal(backtick$F, label$F, tolerance = 1e-8)

  ratio <- pooled_wald(fits, "ratio == 1", df.com = 80)
  ratio_oracle <- mitml::testConstraints(
    fits, constraints = "`y1~x` / `y2~x` - 1", method = "D1", df.com = 80
  )$test
  expect_equal(unname(ratio$F), unname(ratio_oracle[, "F.value"]), tolerance = 1e-8)
  expect_equal(unname(ratio$df1), unname(ratio_oracle[, "df1"]), tolerance = 1e-8)
  expect_equal(unname(ratio$df2), unname(ratio_oracle[, "df2"]), tolerance = 1e-8)
  expect_equal(unname(ratio$pvalue), unname(ratio_oracle[, "P(>F)"]), tolerance = 1e-8)
})

# Independent copy of lavaan.mi::calculate.D2, written from the published formula.
d2_reference <- function(w, DF, asymptotic = FALSE) {
  n_imps <- sum(!is.na(w))
  w_bar <- mean(w, na.rm = TRUE)
  ariv <- (1 + 1 / n_imps) * stats::var(sqrt(w), na.rm = TRUE)
  f_stat <- (w_bar / DF - ((n_imps + 1) * ariv / (n_imps - 1))) / (1 + ariv)
  if (is.finite(f_stat) && f_stat < 0) {
    f_stat <- 0
  }
  fmi <- ariv / (1 + ariv)
  if (asymptotic) {
    return(c(
      chisq = f_stat * DF, df = DF,
      pvalue = stats::pchisq(f_stat * DF, DF, lower.tail = FALSE),
      ariv = ariv, fmi = fmi
    ))
  }
  df2 <- DF^(-3 / n_imps) * (n_imps - 1) * (1 + 1 / ariv)^2
  c(
    F = f_stat, df1 = DF, df2 = df2,
    pvalue = stats::pf(f_stat, DF, df2, lower.tail = FALSE),
    ariv = ariv, fmi = fmi
  )
}

test_that("D2 matches the hand formula and lavaan.mi::calculate.D2", {
  # w = c(3.2, 4.1, 2.7, 5, 3.8), DF = 2.
  # mean(w) = 3.76, ariv = 1.2 * var(sqrt(w)), F floored at 0.
  w <- c(3.2, 4.1, 2.7, 5.0, 3.8)
  DF <- 2
  hand_f <- c(
    F = 1.68396541222381, df1 = 2, df2 = 784.530779480324,
    pvalue = 0.186306678080509, ariv = 0.0615693207669844,
    fmi = 0.0579983987503525
  )
  hand_chisq <- c(
    chisq = 3.36793082444762, df = 2, pvalue = 0.185636389779922,
    ariv = 0.0615693207669844, fmi = 0.0579983987503525
  )
  rubin_d2 <- getFromNamespace("rubin_d2", "svysem")
  got <- rubin_d2(w, DF, asymptotic = FALSE)
  got_chi <- rubin_d2(w, DF, asymptotic = TRUE)
  ref <- d2_reference(w, DF, asymptotic = FALSE)
  ref_chi <- d2_reference(w, DF, asymptotic = TRUE)
  expect_equal(unname(got$F), unname(hand_f["F"]), tolerance = 1e-12)
  expect_equal(unname(got$df2), unname(hand_f["df2"]), tolerance = 1e-12)
  expect_equal(unname(got$pvalue), unname(hand_f["pvalue"]), tolerance = 1e-12)
  expect_equal(unname(got$ariv), unname(hand_f["ariv"]), tolerance = 1e-12)
  expect_equal(unname(got$fmi), unname(hand_f["fmi"]), tolerance = 1e-12)
  expect_equal(unname(got_chi$chisq), unname(hand_chisq["chisq"]), tolerance = 1e-12)
  expect_equal(unname(got_chi$pvalue), unname(hand_chisq["pvalue"]), tolerance = 1e-12)
  expect_equal(unname(got$F), unname(ref["F"]), tolerance = 1e-12)
  expect_equal(unname(got_chi$chisq), unname(ref_chi["chisq"]), tolerance = 1e-12)

  skip_if_not_installed("lavaan.mi")
  oracle <- unclass(lavaan.mi::calculate.D2(w, DF, asymptotic = FALSE))
  oracle_chi <- unclass(lavaan.mi::calculate.D2(w, DF, asymptotic = TRUE))
  expect_equal(unname(got$F), unname(oracle["F"]), tolerance = 1e-10)
  expect_equal(unname(got$df2), unname(oracle["df2"]), tolerance = 1e-10)
  expect_equal(unname(got$pvalue), unname(oracle["pvalue"]), tolerance = 1e-10)
  expect_equal(unname(got$ariv), unname(oracle["ariv"]), tolerance = 1e-10)
  expect_equal(unname(got$fmi), unname(oracle["fmi"]), tolerance = 1e-10)
  expect_equal(unname(got_chi$chisq), unname(oracle_chi["chisq"]), tolerance = 1e-10)
  expect_equal(unname(got_chi$pvalue), unname(oracle_chi["pvalue"]), tolerance = 1e-10)

  # A negative raw F is floored at 0, as in calculate.D2.
  w_neg <- c(0.1, 10, 0.2, 12)
  floored <- rubin_d2(w_neg, DF = 5, asymptotic = TRUE)
  oracle_neg <- unclass(lavaan.mi::calculate.D2(w_neg, 5, asymptotic = TRUE))
  expect_equal(unname(floored$chisq), 0)
  expect_equal(unname(floored$chisq), unname(oracle_neg["chisq"]), tolerance = 1e-10)
  expect_equal(unname(floored$ariv), unname(oracle_neg["ariv"]), tolerance = 1e-10)
})

test_that("standardized correlations are pooled on the Fisher z scale", {
  set.seed(12)
  model <- "
    f1 =~ y1 + y2 + y3
    f2 =~ y4 + y5 + y6
  "
  rhos <- c(0.2, 0.5, 0.8)
  fits <- lapply(rhos, function(rho) {
    n <- 250
    f1 <- stats::rnorm(n)
    f2 <- rho * f1 + sqrt(1 - rho^2) * stats::rnorm(n)
    d <- data.frame(
      y1 = f1 + stats::rnorm(n, sd = 0.5),
      y2 = 0.8 * f1 + stats::rnorm(n, sd = 0.5),
      y3 = 0.6 * f1 + stats::rnorm(n, sd = 0.5),
      y4 = f2 + stats::rnorm(n, sd = 0.5),
      y5 = 0.8 * f2 + stats::rnorm(n, sd = 0.5),
      y6 = 0.7 * f2 + stats::rnorm(n, sd = 0.5)
    )
    lavaan::cfa(model, data = d)
  })
  got <- pooled_estimates(fits, standardized = TRUE)
  pull <- function(lhs, op, rhs, column) {
    vapply(fits, function(fit) {
      ss <- lavaan::standardizedSolution(fit, type = "std.all")
      ss[[column]][ss$lhs == lhs & ss$op == op & ss$rhs == rhs]
    }, numeric(1))
  }
  r <- pull("f1", "~~", "f2", "est.std")
  se_r <- pull("f1", "~~", "f2", "se")
  z <- atanh(r)
  se_z <- se_r / (1 - r^2)
  m <- length(z)
  ubar <- mean(se_z^2)
  b <- stats::var(z)
  se_z_pool <- sqrt(ubar + (1 + 1 / m) * b)
  r_bar <- tanh(mean(z))
  se_r_pool <- se_z_pool * (1 - r_bar^2)
  expect_equal(got["f1~~f2", "est"], r_bar, tolerance = 1e-8)
  expect_equal(got["f1~~f2", "se"], se_r_pool, tolerance = 1e-8)
  expect_gt(abs(got["f1~~f2", "est"] - mean(r)), 1e-4)

  load <- pull("f1", "=~", "y2", "est.std")
  expect_equal(got["f1=~y2", "est"], mean(load), tolerance = 1e-10)
  resid <- pull("y1", "~~", "y1", "est.std")
  expect_equal(got["y1~~y1", "est"], mean(resid), tolerance = 1e-10)
})

test_that("fitMeasures pools the scaled chi-square with D2 and averages indices", {
  skip_if_not_installed("lavaan.mi")
  set.seed(13)
  fits <- lapply(seq_len(4), function(i) {
    n <- 180
    f <- stats::rnorm(n)
    d <- data.frame(
      y1 = f + stats::rnorm(n),
      y2 = (0.8 + 0.02 * i) * f + stats::rnorm(n),
      y3 = 0.7 * f + stats::rnorm(n),
      y4 = 0.6 * f + stats::rnorm(n),
      w = stats::runif(n, 0.5, 1.5)
    )
    lavaan::cfa(
      "f =~ y1 + y2 + y3 + y4", data = d,
      estimator = "MLM", sampling.weights = "w"
    )
  })
  status <- data.frame(
    no.conv = rep(FALSE, 4),
    no.se = rep(FALSE, 4),
    no.npd = rep(FALSE, 4),
    message = rep("", 4)
  )
  obj <- methods::new(
    "svysem_mi",
    fits = fits, df = 60, nobs_used = 180, status = status,
    estimator = "MLM", svysem_call = quote(svysem(model, design))
  )
  fm <- lavaan::fitMeasures(obj)
  w <- vapply(fits, function(fit) {
    lavaan::fitMeasures(fit)[["chisq.scaled"]]
  }, numeric(1))
  DF <- unname(lavaan::fitMeasures(fits[[1]])[["df.scaled"]])
  oracle <- unclass(lavaan.mi::calculate.D2(w, DF, asymptotic = TRUE))
  expect_equal(unname(fm[["chisq.scaled"]]), unname(oracle[["chisq"]]), tolerance = 1e-8)
  expect_equal(unname(fm[["df.scaled"]]), DF, tolerance = 1e-8)
  expect_equal(unname(fm[["pvalue.scaled"]]), unname(oracle[["pvalue"]]), tolerance = 1e-8)
  expect_equal(unname(fm[["ariv"]]), unname(oracle[["ariv"]]), tolerance = 1e-8)
  expect_equal(unname(fm[["fmi"]]), unname(oracle[["fmi"]]), tolerance = 1e-8)
  cfi <- vapply(fits, function(fit) lavaan::fitMeasures(fit)[["cfi.robust"]], numeric(1))
  srmr <- vapply(fits, function(fit) lavaan::fitMeasures(fit)[["srmr"]], numeric(1))
  expect_equal(unname(fm[["cfi.robust.mean"]]), mean(cfi), tolerance = 1e-10)
  expect_equal(unname(fm[["cfi.robust.min"]]), min(cfi), tolerance = 1e-10)
  expect_equal(unname(fm[["cfi.robust.max"]]), max(cfi), tolerance = 1e-10)
  expect_equal(unname(fm[["srmr.mean"]]), mean(srmr), tolerance = 1e-10)
  expect_equal(unname(fm[["srmr.min"]]), min(srmr), tolerance = 1e-10)
  expect_equal(unname(fm[["srmr.max"]]), max(srmr), tolerance = 1e-10)
  as_df <- lavaan::fitMeasures(obj, output = "data.frame")
  expect_s3_class(as_df, "data.frame")
  expect_equal(nrow(as_df), 1)
  expect_equal(unname(as_df$chisq.scaled), unname(fm[["chisq.scaled"]]), tolerance = 1e-12)
})

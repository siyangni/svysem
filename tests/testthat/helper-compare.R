pe_row <- function(fit, lhs, op, rhs) {
  # lavaan skips standard errors unless se = TRUE when class[1] is not "lavaan".
  pe <- lavaan::parameterEstimates(fit, se = TRUE)
  keep <- pe$lhs == lhs & pe$op == op & pe$rhs == rhs
  if (sum(keep) != 1L) {
    stop("Expected one parameter row for ", lhs, op, rhs, call. = FALSE)
  }
  pe[keep, , drop = FALSE]
}

pe_est <- function(fit, lhs, op, rhs) {
  pe_row(fit, lhs, op, rhs)$est
}

pe_se <- function(fit, lhs, op, rhs) {
  pe_row(fit, lhs, op, rhs)$se
}

regression_parts <- function(fit) {
  list(
    est = c(
      pe_est(fit, "y", "~", "x1"),
      pe_est(fit, "y", "~", "x2"),
      pe_est(fit, "y", "~1", "")
    ),
    se = c(
      pe_se(fit, "y", "~", "x1"),
      pe_se(fit, "y", "~", "x2"),
      pe_se(fit, "y", "~1", "")
    )
  )
}

expect_rel <- function(got, expected, tolerance = 1e-6) {
  if (is.null(dim(got)) && is.null(dim(expected)) &&
        !is.null(names(got)) && !is.null(names(expected))) {
    got <- got[names(expected)]
  }
  got <- as.numeric(got)
  expected <- as.numeric(expected)
  expect_identical(length(got), length(expected))
  scale <- pmax(abs(got), abs(expected), 1e-8)
  expect_lte(max(abs(got - expected) / scale), tolerance)
}

glm_parts <- function(design) {
  fitted <- survey::svyglm(y ~ x1 + x2, design = design)
  list(
    est = unname(stats::coef(fitted)[c("x1", "x2", "(Intercept)")]),
    se = unname(sqrt(diag(stats::vcov(fitted)))[c("x1", "x2", "(Intercept)")])
  )
}

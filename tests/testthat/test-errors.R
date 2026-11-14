expect_svysem_error <- function(expr, class) {
  err <- tryCatch(expr, error = function(e) e)
  expect_s3_class(err, class)
  expect_s3_class(err, "svysem_error")
  expect_identical(class(err)[[1L]], class)
  err
}

test_that("an unimplemented estimator is a classed error", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  expect_svysem_error(
    svysem("y ~ x1", design, estimator = "WLSMV"),
    "svysem_error_estimator"
  )
  expect_snapshot(error = TRUE, svysem("y ~ x1", design, estimator = "WLSMV"))
  expect_svysem_error(
    svycfa("f =~ y1 + y2 + y3", design, estimator = "MLR"),
    "svysem_error_estimator"
  )
})

test_that("missing = ml names MLR and is refused for MLM", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  expect_svysem_error(
    svysem("y ~ x1", design, missing = "ml"),
    "svysem_error_missing"
  )
  expect_snapshot(error = TRUE, svysem("y ~ x1", design, missing = "ml"))
  expect_svysem_error(
    svygrowth("y1 ~ 1", design, missing = "fiml"),
    "svysem_error_missing"
  )
})

test_that("an object that is not a survey design is refused", {
  expect_svysem_error(
    svysem("y ~ x1", data.frame(y = 1, x1 = 1)),
    "svysem_error_design"
  )
  expect_snapshot(error = TRUE, svysem("y ~ x1", data.frame(y = 1, x1 = 1)))
  expect_svysem_error(
    svycfa("y1 ~ y2", data.frame(y1 = 1, y2 = 1)),
    "svysem_error_design"
  )
  expect_svysem_error(
    svygrowth("y1 ~ y2", data.frame(y1 = 1, y2 = 1)),
    "svysem_error_design"
  )
})

test_that("a model variable missing from the design is refused", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  expect_svysem_error(
    svysem("y ~ not_in_design", design),
    "svysem_error_variables"
  )
  expect_snapshot(error = TRUE, svysem("y ~ not_in_design", design))
})

test_that("negative weights and a design with no positive weight are refused", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  negative <- design
  negative$prob <- design$prob * 1
  negative$prob[1] <- -1
  expect_svysem_error(
    svysem("y ~ x1", negative),
    "svysem_error_weights"
  )
  expect_snapshot(error = TRUE, svysem("y ~ x1", negative))

  zeros <- design
  zeros$prob <- rep(Inf, length(design$prob))
  expect_svysem_error(svysem("y ~ x1", zeros), "svysem_error_weights")
})

test_that("ordered indicators are refused", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  expect_svysem_error(
    svysem("y1 ~ y2", design, ordered = "y1"),
    "svysem_error_data"
  )
  expect_snapshot(error = TRUE, svysem("y1 ~ y2", design, ordered = "y1"))

  ordered_data <- data.frame(design$variables, check.names = FALSE, stringsAsFactors = FALSE)
  ordered_data$y1 <- ordered(ordered_data$y1)
  ordered_design <- design
  ordered_design$variables <- ordered_data
  expect_svysem_error(
    svysem("y1 ~ y2", ordered_design),
    "svysem_error_data"
  )
})

test_that("a group= argument is refused", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  expect_svysem_error(
    svysem("y ~ x1", design, group = "strat"),
    "svysem_error_data"
  )
  expect_snapshot(error = TRUE, svysem("y ~ x1", design, group = "strat"))
})

test_that("a cluster= argument is refused", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  expect_svysem_error(
    svysem("y ~ x1", design, cluster = "psu"),
    "svysem_error_data"
  )
  expect_snapshot(error = TRUE, svysem("y ~ x1", design, cluster = "psu"))
})

test_that("multiple groups are refused", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  expect_svysem_error(
    svysem("y ~ x1", design, ngroups = 2),
    "svysem_error_data"
  )
  expect_snapshot(error = TRUE, svysem("y ~ x1", design, ngroups = 2))
})

test_that("missing other than listwise is refused", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  expect_svysem_error(
    svysem("y ~ x1", design, missing = "pairwise"),
    "svysem_error_data"
  )
  expect_snapshot(error = TRUE, svysem("y ~ x1", design, missing = "pairwise"))
})

test_that("conditional.x = TRUE is refused", {
  design <- simulate_design(n_psu = 4, n_per = 4, n_strata = 2, seed = 1)
  expect_svysem_error(
    svysem(
      "y ~ x1", design,
      conditional.x = TRUE, fixed.x = TRUE, meanstructure = TRUE
    ),
    "svysem_error_data"
  )
  expect_snapshot(
    error = TRUE,
    svysem(
      "y ~ x1", design,
      conditional.x = TRUE, fixed.x = TRUE, meanstructure = TRUE
    )
  )
})

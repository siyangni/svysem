svysem_sim_outcomes <- function(data) {
  n <- nrow(data)
  if ("strat" %in% names(data)) {
    key <- interaction(data$strat, data$psu, drop = TRUE)
  } else {
    key <- factor(data$psu)
  }
  shock <- stats::rnorm(nlevels(key))
  data$x1 <- stats::rnorm(n)
  data$x2 <- stats::rnorm(n)
  data$y <- 1 + 0.5 * data$x1 - 0.3 * data$x2 + shock[key] * 0.4 +
    stats::rnorm(n)
  eta <- stats::rnorm(n) + shock[key] * 0.3
  slope <- stats::rnorm(n, 0, 0.25)
  data$y1 <- eta + stats::rnorm(n, 0, 0.4)
  data$y2 <- eta + slope + stats::rnorm(n, 0, 0.4)
  data$y3 <- eta + 2 * slope + stats::rnorm(n, 0, 0.4)
  data$y4 <- eta + 3 * slope + stats::rnorm(n, 0, 0.4)
  data$w <- stats::runif(n, 0.5, 2)
  data
}

svysem_sim_frame <- function(n_psu, n_per, n_strata, nest, lonely) {
  if (isTRUE(lonely)) {
    n_psu <- max(as.integer(n_psu), 2L)
    n_per <- as.integer(n_per)
    psu <- c(
      rep(1L, n_per),
      rep(seq_len(n_psu) + 1L, each = n_per),
      rep(seq_len(n_psu) + n_psu + 1L, each = n_per)
    )
    strat <- c(
      rep(1L, n_per),
      rep(2L, n_psu * n_per),
      rep(3L, n_psu * n_per)
    )
    return(data.frame(psu = psu, strat = strat))
  }
  if (n_strata <= 1L) {
    return(data.frame(psu = rep(seq_len(n_psu), each = n_per)))
  }
  if (isTRUE(nest)) {
    psu <- rep(rep(seq_len(n_psu), each = n_per), times = n_strata)
  } else {
    psu <- rep(seq_len(n_psu * n_strata), each = n_per)
  }
  strat <- rep(seq_len(n_strata), each = n_psu * n_per)
  data.frame(psu = psu, strat = strat)
}

#' Simulate a small survey design.
#'
#' Draws a single-stage cluster sample with optional strata. The columns are
#' a regression outcome `y`, covariates `x1` and `x2`, four repeated measures
#' `y1` to `y4` from a linear growth curve, a positive weight `w`, and the
#' sampling unit identifiers. Intended for examples and tests, not as a
#' population generator for applied analyses.
#'
#' @param n_psu Primary sampling units in each stratum. For `lonely = TRUE`
#'   this is the size of each non-lonely stratum, and at least 2.
#' @param n_per Observations in each primary sampling unit.
#' @param n_strata Number of strata. A value below 2 gives an unstratified
#'   design, except when `lonely = TRUE`, which always uses three strata.
#' @param nest If `TRUE`, primary sampling unit labels restart in every
#'   stratum and [survey::svydesign()] is called with `nest = TRUE`. Requires
#'   at least two strata.
#' @param lonely If `TRUE`, the first stratum contains one primary sampling
#'   unit, so [survey::svymean()] follows `survey.lonely.psu`.
#' @param domain If `TRUE`, the design is passed through [subset()] so that
#'   every row with the smallest primary sampling unit label is out of the
#'   domain. With `nest = TRUE` that drops one unit from each stratum and
#'   keeps the original sample size of units in the variance.
#' @param seed Optional seed. The random generator is not changed when
#'   `seed` is `NULL`.
#' @return A `survey.design2` object.
#' @references Oberski D (2014). "lavaan.survey: An R Package for Complex
#'   Survey Analysis of Structural Equation Models." Journal of Statistical
#'   Software, 57(1), 1--27. \doi{10.18637/jss.v057.i01}
#' @seealso [svysem()], [survey::svydesign()]
#' @examples
#' design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 1)
#' survey::svymean(~y, design)
#' @keywords internal
#' @export
simulate_design <- function(n_psu = 10L, n_per = 5L, n_strata = 1L,
                            nest = FALSE, lonely = FALSE, domain = FALSE,
                            seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_psu <- as.integer(n_psu)
  n_per <- as.integer(n_per)
  n_strata <- as.integer(n_strata)
  if (length(n_psu) != 1L || length(n_per) != 1L || length(n_strata) != 1L ||
        is.na(n_psu) || is.na(n_per) || is.na(n_strata) ||
        n_psu < 1L || n_per < 1L || n_strata < 0L) {
    svysem_stop(
      "n_psu and n_per must be positive integers, and n_strata must be non-negative.",
      "svysem_error_data"
    )
  }
  if (isTRUE(nest) && n_strata < 2L && !isTRUE(lonely)) {
    svysem_stop(
      "nest = TRUE requires at least two strata.",
      "svysem_error_data"
    )
  }
  data <- svysem_sim_outcomes(svysem_sim_frame(
    n_psu = n_psu,
    n_per = n_per,
    n_strata = n_strata,
    nest = nest,
    lonely = lonely
  ))
  stratified <- "strat" %in% names(data)
  if (stratified) {
    design <- survey::svydesign(
      ids = ~psu,
      strata = ~strat,
      weights = ~w,
      data = data,
      nest = isTRUE(nest) && !isTRUE(lonely)
    )
  } else {
    design <- survey::svydesign(ids = ~psu, weights = ~w, data = data)
  }
  if (isTRUE(domain)) {
    keep <- data$psu != min(data$psu)
    design <- subset(design, keep)
  }
  design
}

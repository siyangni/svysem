# Oracles from lavaan.survey 1.1.3.1. The test reads the rds files and does
# not call lavaan.survey. Run from anywhere:
#   Rscript data-raw/make_fixtures.R
# The package root is the parent of this script. The working directory is
# not changed.

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
  }
  normalizePath("data-raw/make_fixtures.R", mustWork = TRUE)
}

root <- normalizePath(file.path(dirname(script_path()), ".."), mustWork = TRUE)
out_dir <- file.path(root, "tests", "testthat", "fixtures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(lavaan)
library(survey)
library(lavaan.survey)

versions <- list(
  R = paste(R.version$major, R.version$minor, sep = "."),
  lavaan = as.character(packageVersion("lavaan")),
  survey = as.character(packageVersion("survey")),
  lavaan.survey = as.character(packageVersion("lavaan.survey"))
)

cluster_frame <- function(seed, n_strata = 5L, n_psu = 8L, n_per = 5L) {
  set.seed(seed)
  psu <- rep(rep(seq_len(n_psu), each = n_per), times = n_strata)
  strat <- rep(seq_len(n_strata), each = n_psu * n_per)
  key <- interaction(strat, psu, drop = TRUE)
  list(
    psu = psu,
    strat = strat,
    key = key,
    n = length(psu),
    n_key = nlevels(key)
  )
}

add_weights <- function(data) {
  data$w <- stats::runif(nrow(data), 0.4, 2.5)
  data
}

piecewise_data <- function(seed) {
  base <- cluster_frame(seed)
  key <- base$key
  i_c <- stats::rnorm(base$n_key, 0, 0.6)
  s1_c <- stats::rnorm(base$n_key, 0, 0.2)
  s2_c <- stats::rnorm(base$n_key, 0, 0.15)
  i <- i_c[key] + stats::rnorm(base$n, 0, 0.3)
  s1 <- 0.25 + s1_c[key] + stats::rnorm(base$n, 0, 0.08)
  s2 <- -0.1 + s2_c[key] + stats::rnorm(base$n, 0, 0.08)
  noise <- 0.8
  data <- data.frame(
    y1 = i + stats::rnorm(base$n, 0, noise),
    y2 = i + s1 + stats::rnorm(base$n, 0, noise),
    y3 = i + 2 * s1 + s2 + stats::rnorm(base$n, 0, noise),
    y4 = i + 2 * s1 + 2 * s2 + stats::rnorm(base$n, 0, noise),
    psu = base$psu,
    strat = base$strat
  )
  add_weights(data)
}

bivariate_data <- function(seed) {
  base <- cluster_frame(seed)
  key <- base$key
  iy_c <- stats::rnorm(base$n_key)
  ix_c <- 0.4 * iy_c + stats::rnorm(base$n_key, 0, 0.5)
  sy_c <- stats::rnorm(base$n_key, 0, 0.3)
  sx_c <- 0.3 * sy_c + stats::rnorm(base$n_key, 0, 0.25)
  iy <- iy_c[key] + stats::rnorm(base$n, 0, 0.35)
  ix <- ix_c[key] + stats::rnorm(base$n, 0, 0.35)
  sy <- 0.25 + sy_c[key] + stats::rnorm(base$n, 0, 0.12)
  sx <- -0.15 + sx_c[key] + stats::rnorm(base$n, 0, 0.12)
  data <- data.frame(
    y1 = iy + stats::rnorm(base$n, 0, 0.3),
    y2 = iy + sy + stats::rnorm(base$n, 0, 0.3),
    y3 = iy + 2 * sy + stats::rnorm(base$n, 0, 0.3),
    x1 = ix + stats::rnorm(base$n, 0, 0.3),
    x2 = ix + sx + stats::rnorm(base$n, 0, 0.3),
    x3 = ix + 2 * sx + stats::rnorm(base$n, 0, 0.3),
    psu = base$psu,
    strat = base$strat
  )
  add_weights(data)
}

riclpm_data <- function(seed) {
  base <- cluster_frame(seed)
  key <- base$key
  rix_c <- stats::rnorm(base$n_key)
  riy_c <- 0.45 * rix_c + stats::rnorm(base$n_key, 0, 0.6)
  rix <- rix_c[key] + stats::rnorm(base$n, 0, 0.3)
  riy <- riy_c[key] + stats::rnorm(base$n, 0, 0.3)
  wx <- matrix(0, base$n, 3)
  wy <- matrix(0, base$n, 3)
  wx[, 1] <- stats::rnorm(base$n)
  wy[, 1] <- 0.35 * wx[, 1] + stats::rnorm(base$n)
  wx[, 2] <- 0.3 * wx[, 1] + 0.2 * wy[, 1] + stats::rnorm(base$n, 0, 0.6)
  wy[, 2] <- 0.25 * wy[, 1] + 0.15 * wx[, 1] + stats::rnorm(base$n, 0, 0.6)
  wx[, 3] <- 0.3 * wx[, 2] + 0.2 * wy[, 2] + stats::rnorm(base$n, 0, 0.6)
  wy[, 3] <- 0.25 * wy[, 2] + 0.15 * wx[, 2] + stats::rnorm(base$n, 0, 0.6)
  data <- data.frame(
    x1 = rix + wx[, 1],
    x2 = rix + wx[, 2],
    x3 = rix + wx[, 3],
    y1 = riy + wy[, 1],
    y2 = riy + wy[, 2],
    y3 = riy + wy[, 3],
    psu = base$psu,
    strat = base$strat
  )
  add_weights(data)
}

piecewise_model <- "
  i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
  s1 =~ 0*y1 + 1*y2 + 2*y3 + 2*y4
  s2 =~ 0*y1 + 0*y2 + 1*y3 + 2*y4
  i ~ mi*1
  s1 ~ ms1*1
  y1 ~~ v1*y1
  y4 ~~ v4*y4
  dmean := mi - ms1
  dvar := v1 - v4
"

bivariate_model <- "
  iy =~ 1*y1 + 1*y2 + 1*y3
  sy =~ 0*y1 + 1*y2 + 2*y3
  ix =~ 1*x1 + 1*x2 + 1*x3
  sx =~ 0*x1 + 1*x2 + 2*x3
  iy ~ miy*1
  ix ~ mix*1
  y1 ~~ v1*y1
  x1 ~~ u1*x1
  y1 ~~ x1
  y2 ~~ x2
  y3 ~~ x3
  y1 ~~ x2
  y2 ~~ x3
  dmean := miy - mix
  dvar := v1 - u1
"

riclpm_model <- "
  ri_x =~ 1*x1 + 1*x2 + 1*x3
  ri_y =~ 1*y1 + 1*y2 + 1*y3
  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3
  x1 ~~ 0*x1
  x2 ~~ 0*x2
  x3 ~~ 0*x3
  y1 ~~ 0*y1
  y2 ~~ 0*y2
  y3 ~~ 0*y3
  ri_x ~~ ri_y
  wx1 ~~ wy1
  wx2 ~~ wy2
  wx3 ~~ wy3
  wx2 ~ a*wx1 + b*wy1
  wy2 ~ c*wy1 + d*wx1
  wx3 ~ a*wx2 + b*wy2
  wy3 ~ c*wy2 + d*wx2
  x1 ~ mx*1
  x2 ~ mx*1
  x3 ~ mx*1
  y1 ~ my*1
  y2 ~ my*1
  y3 ~ my*1
  x1 ~~ y2
  wx1 ~~ vx*wx1
  wy1 ~~ vy*wy1
  dmean := my - mx
  dvar := vx - vy
"

kind_of <- function(lhs, op) {
  ifelse(op != ":=", "", ifelse(lhs == "dmean", "mean", "variance"))
}

one_oracle <- function(name, model, data, fitter, lavaan_args) {
  design <- svydesign(
    ids = ~psu, strata = ~strat, weights = ~w, data = data, nest = TRUE
  )
  if (fitter == "growth") {
    naive <- do.call(
      growth,
      c(list(model = model, data = data, sampling.weights = "w"), lavaan_args)
    )
  } else {
    naive <- do.call(
      sem,
      c(list(model = model, data = data, sampling.weights = "w"), lavaan_args)
    )
  }
  if (!lavInspect(naive, "converged")) {
    stop("Naive fit did not converge: ", name, call. = FALSE)
  }
  surveyed <- lavaan.survey(naive, design, estimator = "MLM")
  if (!lavInspect(surveyed, "converged")) {
    stop("lavaan.survey fit did not converge: ", name, call. = FALSE)
  }
  chisq_scaled <- as.numeric(fitMeasures(surveyed, "chisq.scaled"))
  if (!is.finite(chisq_scaled) || chisq_scaled <= 0) {
    stop("Scaled chi-square is not finite and positive: ", name, call. = FALSE)
  }
  pe_check <- parameterEstimates(surveyed)
  defined <- pe_check$op == ":="
  if (any(!is.finite(pe_check$se[defined])) || any(pe_check$se[defined] <= 0)) {
    stop("Defined-parameter standard errors are not positive: ", name, call. = FALSE)
  }
  pe <- parameterEstimates(surveyed)
  std <- standardizedSolution(surveyed)
  parameters <- data.frame(
    lhs = pe$lhs,
    op = pe$op,
    rhs = pe$rhs,
    est = pe$est,
    se = pe$se,
    kind = kind_of(pe$lhs, pe$op),
    stringsAsFactors = FALSE
  )
  standardized <- data.frame(
    lhs = std$lhs,
    op = std$op,
    rhs = std$rhs,
    est.std = std$est.std,
    se = std$se,
    stringsAsFactors = FALSE
  )
  list(
    name = name,
    versions = versions,
    note = paste(
      "Point estimates are weighted ML.",
      "lavaan.survey standard errors use svyvar (n / (n - 1)) and a",
      "block-diagonal Gamma, so ~~ standard errors differ from svysem MLM",
      "by n / (n - 1) and ~1 standard errors do not."
    ),
    model = model,
    fitter = fitter,
    lavaan_args = lavaan_args,
    data = data,
    ids = "psu",
    strata = "strat",
    weights = "w",
    nest = TRUE,
    n_nonzero = sum(weights(design, "sampling") != 0),
    coef = coef(surveyed),
    vcov = vcov(surveyed),
    chisq_scaled = chisq_scaled,
    parameters = parameters,
    standardized = standardized
  )
}

growth_args <- list(fixed.x = FALSE)
sem_args <- list(fixed.x = FALSE, meanstructure = TRUE)

oracles <- list(
  piecewise_growth = one_oracle(
    "piecewise_growth", piecewise_model, piecewise_data(20261113), "growth", growth_args
  ),
  bivariate_growth = one_oracle(
    "bivariate_growth", bivariate_model, bivariate_data(20261114), "growth", growth_args
  ),
  riclpm = one_oracle(
    "riclpm", riclpm_model, riclpm_data(20261115), "sem", sem_args
  )
)

for (name in names(oracles)) {
  saveRDS(oracles[[name]], file.path(out_dir, paste0(name, ".rds")))
  message(
    name, ": chisq.scaled = ", signif(oracles[[name]]$chisq_scaled, 6),
    ", n = ", oracles[[name]]$n_nonzero,
    ", converged coef length ", length(oracles[[name]]$coef)
  )
}

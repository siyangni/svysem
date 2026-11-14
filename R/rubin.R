# Rubin's rules, Barnard-Rubin degrees of freedom, and the D1 and D2
# combining rules. The algebra matches mice 3.19 pool.scalar / barnard.rubin
# and mitml 0.4-5 .pool.estimates / .D1. D2 matches lavaan.mi::calculate.D2.
# Nothing here is exported. svysem_stop() lives in checks.R.

# NULL and Inf both mean "no finite complete-data df".
svysem_check_df_com <- function(df.com) {
  if (is.null(df.com) || is.infinite(df.com)) {
    return(NULL)
  }
  if (length(df.com) != 1L || !is.finite(df.com) || df.com <= 0) {
    svysem_stop(
      "df.com must be NULL, Inf, or one positive finite number.",
      "svysem_error_pooling"
    )
  }
  as.numeric(df.com)
}

# Qhat is p by m (columns are imputations). Uhat is p by p by m.
# When riv is 0 the Barnard-Rubin limit is used: df is infinite if df.com
# is null, otherwise vobs = df.com * (df.com + 1) / (df.com + 3), and
# fmi is 2 / (df + 3). That is mice::barnard.rubin, not a hard zero.
rubin_pool <- function(Qhat, Uhat, df.com = NULL) {
  Qhat <- as.matrix(Qhat)
  m <- ncol(Qhat)
  p <- nrow(Qhat)
  if (m < 2L) {
    svysem_stop(
      "At least two imputations are required for pooling.",
      "svysem_error_pooling"
    )
  }
  if (length(dim(Uhat)) != 3L || any(dim(Uhat)[1:3] != c(p, p, m))) {
    svysem_stop(
      "Uhat must be a p by p by m array matching Qhat.",
      "svysem_error_pooling"
    )
  }
  if (!all(is.finite(Qhat)) || !all(is.finite(Uhat))) {
    svysem_stop(
      "Estimates and covariance matrices must be finite.",
      "svysem_error_pooling"
    )
  }
  df.com <- svysem_check_df_com(df.com)
  nms <- rownames(Qhat)
  Qbar <- rowMeans(Qhat)
  names(Qbar) <- nms
  Ubar <- apply(Uhat, c(1, 2), mean)
  if (!is.matrix(Ubar)) {
    Ubar <- matrix(Ubar, p, p)
  }
  dev <- sweep(Qhat, 1, Qbar, "-")
  B <- tcrossprod(dev) / (m - 1)
  total <- Ubar + (1 + 1 / m) * B
  if (!is.null(nms)) {
    dimnames(Ubar) <- list(nms, nms)
    dimnames(B) <- list(nms, nms)
    dimnames(total) <- list(nms, nms)
  }
  se <- sqrt(diag(total))
  bdiag <- diag(B)
  udiag <- diag(Ubar)
  riv <- (1 + 1 / m) * bdiag / udiag
  # 0/0 is NaN. No between-imputation variance is riv = 0 either way.
  zero_b <- is.finite(bdiag) & bdiag == 0
  riv[zero_b] <- 0
  if (any(!is.finite(riv))) {
    svysem_stop(
      "A within-imputation variance is zero while the estimates still vary.",
      "svysem_error_pooling"
    )
  }
  lambda <- riv / (1 + riv)
  lambda[riv == 0] <- 0
  vm <- rep(Inf, p)
  positive <- riv > 0
  vm[positive] <- (m - 1) * (1 + 1 / riv[positive])^2
  df <- vm
  if (!is.null(df.com)) {
    vobs <- (1 - lambda) * ((df.com + 1) / (df.com + 3)) * df.com
    both <- is.finite(vm) & is.finite(vobs) & vobs > 0
    df[both] <- 1 / (1 / vm[both] + 1 / vobs[both])
    only_vobs <- !is.finite(vm) & is.finite(vobs)
    df[only_vobs] <- vobs[only_vobs]
  }
  names(se) <- nms
  names(df) <- nms
  names(riv) <- nms
  names(lambda) <- nms
  tstat <- Qbar / se
  fmi <- (riv + 2 / (df + 3)) / (riv + 1)
  pval <- 2 * (1 - stats::pt(abs(tstat), df = df))
  list(
    est = Qbar, se = se, t = tstat, df = df, p = pval,
    riv = riv, lambda = lambda, fmi = fmi, m = m,
    total = total, Ubar = Ubar, B = B
  )
}

# Qhat is k by m of already transformed constraint values.
# Uhat is k by k by m. r = 0 means the constraints do not vary.
rubin_d1 <- function(Qhat, Uhat, df.com = NULL) {
  Qhat <- as.matrix(Qhat)
  k <- nrow(Qhat)
  m <- ncol(Qhat)
  if (m < 2L) {
    svysem_stop(
      "At least two imputations are required for pooling.",
      "svysem_error_pooling"
    )
  }
  if (k < 1L) {
    svysem_stop("At least one constraint is required.", "svysem_error_constraint")
  }
  if (!all(is.finite(Qhat)) || !all(is.finite(Uhat))) {
    svysem_stop(
      "Constraint estimates and covariances must be finite.",
      "svysem_error_pooling"
    )
  }
  df.com <- svysem_check_df_com(df.com)
  Qbar <- rowMeans(Qhat)
  Ubar <- apply(Uhat, c(1, 2), mean)
  if (!is.matrix(Ubar)) {
    Ubar <- matrix(Ubar, k, k)
  }
  B <- stats::cov(t(Qhat))
  if (!is.matrix(B)) {
    B <- matrix(B, k, k)
  }
  Uinv <- tryCatch(solve(Ubar), error = function(e) NULL)
  if (is.null(Uinv)) {
    svysem_stop(
      "The within-imputation covariance of the constraints is singular.",
      "svysem_error_pooling"
    )
  }
  r <- (1 + 1 / m) * sum(diag(B %*% Uinv)) / k
  if (!is.finite(r) || r == 0) {
    f_stat <- as.numeric(t(Qbar) %*% Uinv %*% Qbar / k)
    return(list(
      F = f_stat, df1 = as.numeric(k), df2 = Inf,
      pvalue = stats::pf(f_stat, k, Inf, lower.tail = FALSE),
      riv = 0
    ))
  }
  t_tilde <- (1 + r) * Ubar
  f_stat <- as.numeric(t(Qbar) %*% solve(t_tilde) %*% Qbar / k)
  t_df <- k * (m - 1)
  if (is.null(df.com)) {
    if (t_df > 4) {
      v <- 4 + (t_df - 4) * (1 + (1 - 2 / t_df) / r)^2
    } else {
      v <- t_df * (1 + 1 / k) * (1 + 1 / r)^2 / 2
    }
  } else {
    # mitml 0.4-5 small-sample denominator, used whenever df.com is supplied.
    a <- r * t_df / (t_df - 2)
    vstar <- ((df.com + 1) / (df.com + 3)) * df.com
    c0 <- 1 / (t_df - 4)
    c1 <- vstar - 2 * (1 + a)
    c2 <- vstar - 4 * (1 + a)
    z <- 1 / c2 +
      c0 * (a^2 * c1 / ((1 + a)^2 * c2)) +
      c0 * (8 * a^2 * c1 / ((1 + a) * c2^2) + 4 * a^2 / ((1 + a) * c2)) +
      c0 * (4 * a^2 / (c2 * c1) + 16 * a^2 * c1 / c2^3) +
      c0 * (8 * a^2 / c2^2)
    v <- 4 + 1 / z
  }
  list(
    F = f_stat, df1 = as.numeric(k), df2 = as.numeric(v),
    pvalue = stats::pf(f_stat, k, v, lower.tail = FALSE),
    riv = r
  )
}

# w is a vector of chi-square statistics. The fmi here is ariv / (1 + ariv),
# the D2 fraction of missing information, not the Barnard-Rubin fmi.
rubin_d2 <- function(w, DF, asymptotic = FALSE) {
  w <- as.numeric(w)
  if (length(DF) != 1L || !is.finite(DF) || DF <= 0) {
    svysem_stop(
      "D2 requires one positive degrees-of-freedom value.",
      "svysem_error_pooling"
    )
  }
  n_imps <- sum(!is.na(w))
  if (n_imps < 2L) {
    svysem_stop(
      "At least two imputations are required for pooling.",
      "svysem_error_pooling"
    )
  }
  w_bar <- mean(w, na.rm = TRUE)
  ariv <- (1 + 1 / n_imps) * stats::var(sqrt(w), na.rm = TRUE)
  f_stat <- (w_bar / DF - ((n_imps + 1) * ariv / (n_imps - 1))) / (1 + ariv)
  if (is.finite(f_stat) && f_stat < 0) {
    f_stat <- 0
  }
  fmi <- ariv / (1 + ariv)
  if (isTRUE(asymptotic)) {
    chisq <- f_stat * DF
    return(list(
      chisq = chisq, df = DF,
      pvalue = stats::pchisq(chisq, df = DF, lower.tail = FALSE),
      ariv = ariv, fmi = fmi
    ))
  }
  df2 <- DF^(-3 / n_imps) * (n_imps - 1) * (1 + 1 / ariv)^2
  list(
    F = f_stat, df1 = DF, df2 = df2,
    pvalue = stats::pf(f_stat, df1 = DF, df2 = df2, lower.tail = FALSE),
    ariv = ariv, fmi = fmi
  )
}

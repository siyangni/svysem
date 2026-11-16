# Multiply imputed design-based SEM

Fits of one model to each completed survey design, pooled afterwards
with Rubin's rules. The design degrees of freedom are the complete-data
degrees of freedom (`df.com`) unless the user overrides them in
[`pooled_estimates()`](https://siyangni.github.io/svysem/reference/pooled_estimates.md)
or
[`pooled_wald()`](https://siyangni.github.io/svysem/reference/pooled_wald.md).

## Usage

``` r
# S4 method for class 'svysem_mi'
show(object)

# S3 method for class 'svysem_mi'
coef(object, ...)

# S4 method for class 'svysem_mi'
coef(object, ...)

# S3 method for class 'svysem_mi'
vcov(object, ...)

# S4 method for class 'svysem_mi'
vcov(object, ...)

# S3 method for class 'svysem_mi'
nobs(object, ...)

# S3 method for class 'svysem_mi'
summary(object, ...)

# S3 method for class 'svysem_mi.summary'
print(x, ...)

# S4 method for class 'svysem_mi'
fitMeasures(
  object,
  fit.measures = "all",
  baseline.model = NULL,
  h1.model = NULL,
  fm.args = list(standard.test = "default", scaled.test = "default", rmsea.ci.level =
    0.9, rmsea.close.h0 = 0.05, rmsea.notclose.h0 = 0.08, robust = TRUE, cat.check.pd =
    TRUE),
  output = "vector",
  ...
)
```

## Arguments

- object:

  A `svysem_mi` object.

- ...:

  Further arguments, ignored.

- x:

  A summary of a `svysem_mi` object.

- fit.measures:

  `"all"` or a character vector of measure names to keep.

- baseline.model:

  Accepted for compatibility with lavaan's generic. Not used: each
  imputation already carries its own baseline.

- h1.model:

  Accepted for compatibility with lavaan's generic. Not used.

- fm.args:

  Passed to
  [`lavaan::fitMeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html)
  on each imputation.

- output:

  `"vector"` (default), or `"matrix"` / `"data.frame"` for a one-row
  data frame.

## Value

An object of class `svysem_mi`.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns the pooled point
estimates, [`vcov()`](https://rdrr.io/r/stats/vcov.html) their total
covariance, and [`summary()`](https://rdrr.io/r/base/summary.html) an
object of class `svysem_mi.summary`.

## Details

The object is returned by
[`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md),
[`svycfa()`](https://siyangni.github.io/svysem/reference/svysem.md) or
[`svygrowth()`](https://siyangni.github.io/svysem/reference/svysem.md)
when `design` is a `svyimputationList`. Each completed design is fitted
separately. `on_fail = "stop"` stops at the first imputation that
errors. `on_fail = "omit"` records that imputation in `status` and pools
the rest, provided at least two fits remain.
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) are the Rubin pooled
estimate and its total covariance. Both are S4 methods, so they still
dispatch after lavaan has turned those generics into S4 methods.

## Slots

- `fits`:

  Successful fits, in imputation order, with failures removed.

- `df`:

  Design degrees of freedom,
  [`survey::degf()`](https://rdrr.io/pkg/survey/man/svychisq.html) of
  the completed designs. This is the default `df.com` when pooling.

- `nobs_used`:

  Number of rows with a positive sampling weight in the first completed
  design.

- `status`:

  One row per imputation, including failures. Columns are `no.conv`,
  `no.se`, `no.npd` and `message` (an empty string when the imputation
  did not error).

- `estimator`:

  Estimator name. `"MLM"`, `"WLSMV"`, or `"MLR"`.

- `svysem_call`:

  The call that created the object.

## References

Rubin D (1987). Multiple Imputation for Nonresponse in Surveys. Wiley,
New York.
[doi:10.1002/9780470316696](https://doi.org/10.1002/9780470316696)

Barnard J, Rubin DB (1999). "Small-Sample Degrees of Freedom with
Multiple Imputation." Biometrika, 86(4), 948–959.
[doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948)

## See also

[`pooled_estimates()`](https://siyangni.github.io/svysem/reference/pooled_estimates.md),
[`pooled_wald()`](https://siyangni.github.io/svysem/reference/pooled_wald.md),
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)

## Examples

``` r
if (requireNamespace("mitools", quietly = TRUE)) {
  base <- simulate_design(
    n_psu = 6, n_per = 4, n_strata = 2, seed = 1
  )$variables
  alt <- base
  alt$y <- alt$y + 0.2
  des <- survey::svydesign(
    ids = ~psu, strata = ~strat, weights = ~w,
    data = mitools::imputationList(list(base, alt))
  )
  fit <- svysem(
    "y ~ x1 + x2", des,
    fixed.x = FALSE, meanstructure = TRUE
  )
  coef(fit)
  fit
}
#> svysem_mi (multiple imputation)
#> Estimator: MLM
#> Imputations: 2 used / 2 supplied
#> Design degrees of freedom: 10
#> Observations: 48
#> Not converged: 0; no SE: 0; not positive definite: 0
```

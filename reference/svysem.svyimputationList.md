# Fit each completed design

[`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md),
[`svycfa()`](https://siyangni.github.io/svysem/reference/svysem.md) and
[`svygrowth()`](https://siyangni.github.io/svysem/reference/svysem.md)
dispatch on a `svyimputationList`. Every completed design is passed to
the same generic, which selects
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html) or
[`lavaan::growth()`](https://rdrr.io/pkg/lavaan/man/growth.html).
`estimator` may be `"MLM"`, `"WLSMV"`, or `"MLR"`. `parallel` is any
function of `(X, FUN)` with the same contract as
[`lapply()`](https://rdrr.io/r/base/lapply.html); `NULL` uses
[`lapply()`](https://rdrr.io/r/base/lapply.html). It loops imputations.
Replicate weights inside one WLSMV fit stay sequential, so the two loops
are not nested. Tests and examples use at most two cores.

## Usage

``` r
# S3 method for class 'svyimputationList'
svysem(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...,
  on_fail = c("stop", "omit")
)

# S3 method for class 'svyimputationList'
svycfa(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...,
  on_fail = c("stop", "omit")
)

# S3 method for class 'svyimputationList'
svygrowth(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...,
  on_fail = c("stop", "omit")
)
```

## Arguments

- model:

  Model syntax, as for
  [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html),
  [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html) or
  [`lavaan::growth()`](https://rdrr.io/pkg/lavaan/man/growth.html).

- design:

  A `svyimputationList`, from
  `svydesign(data = mitools::imputationList(...))`.

- estimator:

  Variance estimator. `"MLM"`, `"WLSMV"`, or `"MLR"`.

- replicates, replicate.type:

  Passed to each completed design. See
  [`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md).
  Replicate weights inside an imputation are not parallelized.

- parallel:

  `NULL` to use [`lapply()`](https://rdrr.io/r/base/lapply.html), or a
  function called as `parallel(X, FUN)`. One example is
  `function(X, FUN) parallel::mclapply(X, FUN, mc.cores = 2)`.

- ...:

  Passed to the fitter for each completed design.

- on_fail:

  `"stop"` (default) or `"omit"`. `"omit"` drops a failed imputation,
  records it in the status table, and continues.

## Value

A
[svysem_mi](https://siyangni.github.io/svysem/reference/svysem_mi-class.md)
object. Its `df` slot is the design degrees of freedom and is the
default `df.com` for pooling.

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
}
#>        y~x1        y~x2        y~~y      x1~~x1      x1~~x2      x2~~x2 
#>  0.44063045 -0.48330670  0.88005445  0.72560561 -0.08520233  0.95830528 
#>         y~1        x1~1        x2~1 
#>  1.00160774  0.08784715  0.10371083 
```

# Design-based SEM fit

An object of class `svysem` is a
[lavaan::lavaan](https://rdrr.io/pkg/lavaan/man/lavaan-class.html) fit
whose standard errors use the survey design. Where the estimator has a
scaled test, that test uses the design too.
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`lavaan::parameterEstimates()`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html),
[`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html),
[`lavaan::lavTestWald()`](https://rdrr.io/pkg/lavaan/man/lavTestWald.html)
and
[`lavaan::fitMeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html)
therefore report design-based results with no further step.

## Usage

``` r
# S3 method for class 'svysem'
coef(object, ...)

# S3 method for class 'svysem'
vcov(object, ...)

# S4 method for class 'svysem'
coef(object, ...)

# S4 method for class 'svysem'
vcov(object, ...)

# S4 method for class 'svysem'
show(object)

# S4 method for class 'svysem'
summary(object, ...)

# S3 method for class 'svysem.summary'
print(x, ...)
```

## Arguments

- object:

  A `svysem` fit.

- ...:

  Passed to the lavaan summary.

- x:

  A summary of a `svysem` fit.

## Value

The class extends `lavaan`. See
[`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md) for
the fitting functions.

## Slots

- `df`:

  Design degrees of freedom,
  [`survey::degf()`](https://rdrr.io/pkg/survey/man/svychisq.html) of
  the design passed to
  [`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md),
  [`svycfa()`](https://siyangni.github.io/svysem/reference/svysem.md) or
  [`svygrowth()`](https://siyangni.github.io/svysem/reference/svysem.md).

- `estimator`:

  Estimator used for the variance. `"MLM"`, `"WLSMV"`, or `"MLR"`.

- `variance_method`:

  How the variance was estimated. `"linearization"` for single-group MLM
  or for MLR on a `survey.design`, and `"replication"` for WLSMV, for
  MLM with `group=`, or for any estimator on a `svyrep.design`.

- `svysem_call`:

  The call to
  [`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md),
  [`svycfa()`](https://siyangni.github.io/svysem/reference/svysem.md) or
  [`svygrowth()`](https://siyangni.github.io/svysem/reference/svysem.md).

## References

Oberski D (2014). "lavaan.survey: An R Package for Complex Survey
Analysis of Structural Equation Models." Journal of Statistical
Software, 57(1), 1–27.
[doi:10.18637/jss.v057.i01](https://doi.org/10.18637/jss.v057.i01)

## See also

[`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md),
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html),
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html)

## Examples

``` r
design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 1)
fit <- svysem(
  "y ~ x1 + x2", design,
  fixed.x = FALSE, meanstructure = TRUE
)
fit@df
#> [1] 10
fit@estimator
#> [1] "MLM"
fit@variance_method
#> [1] "linearization"
```

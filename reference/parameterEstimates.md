# Parameter estimates

Design-based standard errors and confidence intervals for a
[svysem](https://siyangni.github.io/svysem/reference/svysem-class.md)
fit.

## Usage

``` r
parameterEstimates(object, ...)

# S3 method for class 'svysem'
parameterEstimates(object, se = TRUE, ...)

# Default S3 method
parameterEstimates(object, ...)
```

## Arguments

- object:

  A fitted model. A `svysem` fit uses the method below. Any other object
  is passed to
  [`lavaan::parameterEstimates()`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html).

- ...:

  Passed to
  [`lavaan::parameterEstimates()`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html).
  `zstat`, `pvalue`, and `ci` keep lavaan's defaults.

- se:

  Include standard errors. The default is `TRUE`.

## Value

The data frame from
[`lavaan::parameterEstimates()`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html).

## Details

[`lavaan::parameterEstimates()`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html)
omits standard errors, z statistics, p-values, and confidence intervals
when the object's first class is not `"lavaan"`, unless `se = TRUE` is
supplied. A svysem fit keeps the class `svysem` so that
[`summary()`](https://rdrr.io/r/base/summary.html) and `show()` can
report the design. This method supplies `se = TRUE`, so the default
table includes `se`, `z`, `pvalue`, `ci.lower`, and `ci.upper`. Those
standard errors match `sqrt(diag(vcov(object)))`. `se = FALSE` returns
the point estimates only.

[`lavaan::parameterEstimates()`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html)
still applies lavaan's class check. After
[`library(lavaan)`](https://lavaan.ugent.be) is called later than
[`library(svysem)`](https://github.com/siyangni/svysem), call
`svysem::parameterEstimates()`.

## See also

[standardizedSolution()](https://siyangni.github.io/svysem/reference/standardizedSolution.md),
[`lavaan::parameterEstimates()`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html)

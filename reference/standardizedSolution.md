# Standardized solution

Design-based standard errors for the standardized estimates of a
[svysem](https://siyangni.github.io/svysem/reference/svysem-class.md)
fit.

## Usage

``` r
standardizedSolution(object, ...)

# S3 method for class 'svysem'
standardizedSolution(object, se = TRUE, ...)

# Default S3 method
standardizedSolution(object, ...)
```

## Arguments

- object:

  A fitted model. A `svysem` fit uses the method below. Any other object
  is passed to
  [`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html).

- ...:

  Passed to
  [`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html).

- se:

  Include standard errors. The default is `TRUE`.

## Value

The data frame from
[`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html).

## Details

[`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)
uses the same class check as
[`lavaan::parameterEstimates()`](https://rdrr.io/pkg/lavaan/man/parameterEstimates.html).
This method supplies `se = TRUE` for a svysem fit. `se = FALSE` returns
the standardized point estimates only. After
[`library(lavaan)`](https://lavaan.ugent.be) is called later than
[`library(svysem)`](https://github.com/siyangni/svysem), call
`svysem::standardizedSolution()`.

## See also

[parameterEstimates()](https://siyangni.github.io/svysem/reference/parameterEstimates.md),
[`lavaan::standardizedSolution()`](https://rdrr.io/pkg/lavaan/man/standardizedSolution.html)

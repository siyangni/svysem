# Factor scores and plausible values

Posterior means of the latent variables in a
[svysem](https://siyangni.github.io/svysem/reference/svysem-class.md)
fit, using the model-implied moments and the survey-weighted estimates.
`draws = 0` returns one row per case. A positive `draws` returns that
many plausible-value matrices, drawn from the posterior covariance of
the same regression. Cases are scored from the indicators they actually
have, so a full-information fit can leave some indicators missing.

## Usage

``` r
# S3 method for class 'svysem'
predict(object, newdata = NULL, draws = 0, seed = NULL, ...)

# S4 method for class 'svysem'
predict(object, newdata = NULL, draws = 0, seed = NULL, ...)

# S3 method for class 'svysem_mi'
predict(
  object,
  newdata = NULL,
  draws = 0,
  seed = NULL,
  per.imputation = FALSE,
  ...
)

# S4 method for class 'svysem_mi'
predict(
  object,
  newdata = NULL,
  draws = 0,
  seed = NULL,
  per.imputation = FALSE,
  ...
)
```

## Arguments

- object:

  A `svysem` or `svysem_mi` fit with continuous outcomes and one group.

- newdata:

  A data frame with the indicators used in the measurement model. `NULL`
  scores the cases stored in the fit.

- draws:

  Number of plausible-value data sets. `0` returns the posterior means.

- seed:

  Optional seed, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the draws.
  The previous random-generator state is not restored. `NULL` leaves it
  alone.

- ...:

  Unused. Present so the method accepts the lavaan arguments it does not
  use.

- per.imputation:

  If `TRUE`, return one result per imputation instead of the average
  posterior mean. Required when `draws` is positive.

## Value

For `draws = 0`, a numeric matrix with one row per case and one column
per latent variable. For `draws > 0`, a list of that many matrices. A
multiple-imputation fit returns the average of the posterior means, or a
list of per-imputation results when `per.imputation = TRUE`.

## Details

[`lavaan::lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html)
on a weighted fit in lavaan 0.6-21 does not return these values. The
comparison in the examples is the unweighted formula, which is what
[`predict()`](https://rdrr.io/r/stats/predict.html) implements.

For a
[svysem_mi](https://siyangni.github.io/svysem/reference/svysem_mi-class.md)
object the posterior means are averaged across imputations. Plausible
values are not averaged: set `per.imputation = TRUE` to get one set of
matrices per imputation. `newdata` is scored under each imputation's
estimates. Without `newdata`, each imputation is scored on the data
lavaan stored for that fit, and those rows are averaged in imputation
order, so the completed designs need the same row order.

## See also

[`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md),
[`lavaan::lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html)

## Examples

``` r
design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 1)
fit <- svycfa("f =~ y1 + y2 + y3 + y4", design)
predict(fit)[1, , drop = FALSE]
#>           f
#> 1 0.1295933
```

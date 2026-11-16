# Simulate a small survey design.

Draws a single-stage cluster sample with optional strata. The columns
are a regression outcome `y`, covariates `x1` and `x2`, four repeated
measures `y1` to `y4` from a linear growth curve, a positive weight `w`,
and the sampling unit identifiers. Intended for examples and tests, not
as a population generator for applied analyses.

## Usage

``` r
simulate_design(
  n_psu = 10L,
  n_per = 5L,
  n_strata = 1L,
  nest = FALSE,
  lonely = FALSE,
  domain = FALSE,
  seed = NULL
)
```

## Arguments

- n_psu:

  Primary sampling units in each stratum. For `lonely = TRUE` this is
  the size of each non-lonely stratum, and at least 2.

- n_per:

  Observations in each primary sampling unit.

- n_strata:

  Number of strata. A value below 2 gives an unstratified design, except
  when `lonely = TRUE`, which always uses three strata.

- nest:

  If `TRUE`, primary sampling unit labels restart in every stratum and
  [`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)
  is called with `nest = TRUE`. Requires at least two strata.

- lonely:

  If `TRUE`, the first stratum contains one primary sampling unit, so
  [`survey::svymean()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
  follows `survey.lonely.psu`.

- domain:

  If `TRUE`, the design is passed through
  [`subset()`](https://rdrr.io/r/base/subset.html) so that every row
  with the smallest primary sampling unit label is out of the domain.
  With `nest = TRUE` that drops one unit from each stratum and keeps the
  original sample size of units in the variance.

- seed:

  Optional seed. The random generator is not changed when `seed` is
  `NULL`.

## Value

A `survey.design2` object.

## References

Oberski D (2014). "lavaan.survey: An R Package for Complex Survey
Analysis of Structural Equation Models." Journal of Statistical
Software, 57(1), 1–27.
[doi:10.18637/jss.v057.i01](https://doi.org/10.18637/jss.v057.i01)

## See also

[`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md),
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)

## Examples

``` r
design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 1)
survey::svymean(~y, design)
#>      mean     SE
#> y 0.89019 0.2112
```

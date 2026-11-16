# Pooled parameter estimates

Combine point estimates and covariance matrices from a list of
completed-data fits with Rubin's rules. Degrees of freedom use the
Barnard-Rubin adjustment when `df.com` is finite. For a
[svysem_mi](https://siyangni.github.io/svysem/reference/svysem_mi-class.md)
object, `df.com` defaults to the design degrees of freedom.

## Usage

``` r
pooled_estimates(object, standardized = FALSE, level = 0.95, df.com = NULL)
```

## Arguments

- object:

  A list of objects with
  [`stats::coef()`](https://rdrr.io/r/stats/coef.html) and
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html), or a
  [svysem_mi](https://siyangni.github.io/svysem/reference/svysem_mi-class.md)
  object.

- standardized:

  If `TRUE`, pool `std.all` standardized estimates. Correlations use the
  Fisher z scale.

- level:

  Confidence level for `ci.lower` and `ci.upper`.

- df.com:

  Complete-data degrees of freedom. For a `svysem_mi` object the default
  is the design degrees of freedom. `NULL` or `Inf` uses the
  large-sample Rubin degrees of freedom.

## Value

A data frame with one row per parameter and columns `est`, `se`, `t`,
`df`, `p`, `ci.lower`, `ci.upper`, `riv`, `lambda`, `fmi` and `m`.

## Details

`t` is the t statistic, the pooled estimate divided by its standard
error, not the total variance. With `standardized = TRUE`, each fit must
be a lavaan object. Standard errors on that table are requested with
`se = TRUE`, because lavaan omits them when the fit's class is not
exactly `lavaan`. `se = TRUE` matches `sqrt(diag(vcov()))`. Loadings,
regressions and variances are pooled on the `std.all` scale.
Correlations (`~~` with different left- and right-hand sides) are pooled
on the Fisher z scale, `z = atanh(r)` with `se_z = se_r / (1 - r^2)`.
Reported `est`, `se` and confidence limits are transformed back
(`est = tanh(z)`, `se = se_z * (1 - est^2)`, limits
`tanh(z +/- t * se_z)`). For those rows, `t`, `df`, `p`, `riv`, `lambda`
and `fmi` stay on the z scale.

## References

Rubin D (1987). Multiple Imputation for Nonresponse in Surveys. Wiley,
New York.
[doi:10.1002/9780470316696](https://doi.org/10.1002/9780470316696)

Barnard J, Rubin DB (1999). "Small-Sample Degrees of Freedom with
Multiple Imputation." Biometrika, 86(4), 948–959.
[doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948)

## See also

[`pooled_wald()`](https://siyangni.github.io/svysem/reference/pooled_wald.md),
[`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md)

## Examples

``` r
fit_once <- function(seed) {
  design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = seed)
  svysem(
    "y ~ x1 + x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
}
fits <- list(fit_once(1), fit_once(2))
pooled_estimates(fits, df.com = fits[[1]]@df)
#>                est        se           t        df           p     ci.lower
#> y~x1    0.43776222 0.1731844  2.52772350 8.4545273 0.033908647   0.04210485
#> y~x2   -0.52854889 0.1396682 -3.78431789 3.6823825 0.022497424  -0.92988395
#> y~~y    0.89548647 0.1811540  4.94323226 8.2449794 0.001034507   0.47989519
#> x1~~x1  1.07242967 0.6386310  1.67926341 0.5529025 0.461142789 -61.94565379
#> x1~~x2  0.02269523 0.2505618  0.09057736 1.7366016 0.937169858  -1.22829334
#> x2~~x2  1.20146675 0.4620148  2.60049411 0.7195076 0.305149805 -14.12323045
#> y~1     1.04165645 0.2862732  3.63868009 1.0700012 0.157764570  -2.07894494
#> x1~1    0.07388410 0.1736993  0.42535633 8.2717101 0.681427588  -0.32438851
#> x2~1    0.01186196 0.2655119  0.04467580 3.1922618 0.967015091  -0.80504230
#>          ci.upper          riv       lambda       fmi m
#> y~x1    0.8334196 0.0008235493 0.0008228716 0.1752826 2
#> y~x2   -0.1272138 0.4593951153 0.3147846052 0.5198658 2
#> y~~y    1.3110777 0.0222551340 0.0217706258 0.1957557 2
#> x1~~x1 64.0905131 7.6797249378 0.8847889758 0.9496436 2
#> x1~~x2  1.2736838 1.2538117920 0.5563072287 0.7436537 2
#> x2~~x2 16.5261639 4.9170253434 0.8309961607 0.9218705 2
#> y~1     4.1622578 2.5459697437 0.7179896975 0.8565697 2
#> x1~1    0.4721567 0.0197691479 0.0193859050 0.1933815 2
#> x2~1    0.8287662 0.5600757137 0.3590054693 0.5660363 2
```

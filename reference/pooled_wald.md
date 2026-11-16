# Wald test of pooled constraints

Combine a set of constraints with the D1 rule. A constraint is an R
expression in the free-parameter names that is tested against zero.
`a == b` is the scalar `a - b`. Labels and `:=` names from a lavaan fit
are rewritten onto the free parameters and differentiated with the delta
method inside each imputation. Inequality constraints are not supported.

## Usage

``` r
pooled_wald(object, constraints, df.com = NULL)
```

## Arguments

- object:

  A list of objects with
  [`stats::coef()`](https://rdrr.io/r/stats/coef.html) and
  [`stats::vcov()`](https://rdrr.io/r/stats/vcov.html), or a
  [svysem_mi](https://siyangni.github.io/svysem/reference/svysem_mi-class.md)
  object.

- constraints:

  Character vector of constraints. Lavaan names that are not syntactic
  are written in backticks, for example `` "`y~x` - `y~z`" ``.

- df.com:

  Complete-data degrees of freedom. For a `svysem_mi` object the default
  is the design degrees of freedom. `NULL` uses the large-sample D1
  denominator.

## Value

A one-row data frame with columns `F`, `df1`, `df2`, `pvalue` and `riv`.

## Details

For a
[svysem_mi](https://siyangni.github.io/svysem/reference/svysem_mi-class.md)
object, `df.com` defaults to the design degrees of freedom. The
small-sample denominator follows mitml when that value is positive. When
it is undefined (`k * (m - 1) <= 4`) or negative, one constraint uses
the scalar Barnard-Rubin degrees of freedom and several constraints use
a positive harmonic mean of the large-sample denominator and the
observed-data degrees of freedom.

## References

Rubin D (1987). Multiple Imputation for Nonresponse in Surveys. Wiley,
New York.
[doi:10.1002/9780470316696](https://doi.org/10.1002/9780470316696)

Barnard J, Rubin DB (1999). "Small-Sample Degrees of Freedom with
Multiple Imputation." Biometrika, 86(4), 948–959.
[doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948)

## See also

[`pooled_estimates()`](https://siyangni.github.io/svysem/reference/pooled_estimates.md),
[`svysem()`](https://siyangni.github.io/svysem/reference/svysem.md)

## Examples

``` r
fit_once <- function(seed) {
  design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = seed)
  svysem(
    "y ~ b1*x1 + b2*x2", design,
    fixed.x = FALSE, meanstructure = TRUE
  )
}
fits <- list(fit_once(1), fit_once(2))
pooled_wald(fits, constraints = "b1 - b2", df.com = fits[[1]]@df)
#>          F df1      df2      pvalue      riv
#> 1 21.10291   1 6.694175 0.002808073 0.138613
```

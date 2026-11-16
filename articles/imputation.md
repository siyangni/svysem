# Multiply imputed survey data

Each completed data set is fitted on its own design. The fits are pooled
with Rubin’s rules. The design degrees of freedom are the complete-data
degrees of freedom.

``` r

library(svysem)
base <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 2)$variables
second <- base
second$y1 <- second$y1 + 0.15
second$y3 <- second$y3 - 0.10
design <- survey::svydesign(
  ids = ~psu, strata = ~strat, weights = ~w,
  data = mitools::imputationList(list(base, second))
)
```

The model labels a residual covariance and defines its Fisher *z* with
[`atanh()`](https://rdrr.io/r/base/Hyperbolic.html). That defined
parameter is pooled on the *z* scale.

``` r

model <- "
  f =~ y1 + y2 + y3 + y4
  y1 ~~ v1*y1
  y2 ~~ v2*y2
  y1 ~~ c*y2
  z := atanh(c / sqrt(v1 * v2))
"
fit <- svycfa(model, design)
coef(fit)
#>     f=~y2     f=~y3     f=~y4        v1        v2         c    y3~~y3    y4~~y4 
#> 1.1549870 1.1927962 1.4056068 0.4098036 0.1791252 0.1547116 0.1738689 0.2163105 
#>      f~~f 
#> 0.6157304
```

[`pooled_estimates()`](https://siyangni.github.io/svysem/reference/pooled_estimates.md)
reports the Rubin estimate, its standard error, and the Barnard-Rubin
degrees of freedom. `standardized = TRUE` pools correlations on the
Fisher *z* scale and transforms them back.
[`pooled_wald()`](https://siyangni.github.io/svysem/reference/pooled_wald.md)
takes the same constraints as
[`lavaan::lavTestWald()`](https://rdrr.io/pkg/lavaan/man/lavTestWald.html).
`fitMeasures()` pools the scaled chi-square and summarizes the robust
indices by their mean, minimum, and maximum. Those indices are not
pooled with Rubin’s rules.

``` r

head(pooled_estimates(fit), 3)
#>            est         se         t       df            p  ci.lower ci.upper
#> f=~y2 1.154987 0.07868703 14.678239 8.461538 2.590283e-07 0.9752428 1.334731
#> f=~y3 1.192796 0.15386787  7.752081 8.461538 4.020284e-05 0.8413170 1.544275
#> f=~y4 1.405607 0.27205119  5.166700 8.461538 7.174873e-04 0.7841623 2.027051
#>       riv lambda       fmi m
#> f=~y2   0      0 0.1744966 2
#> f=~y3   0      0 0.1744966 2
#> f=~y4   0      0 0.1744966 2
pooled_wald(fit, "z == 0")
#>          F df1      df2      pvalue          riv
#> 1 18.45062   1 8.461538 0.002311303 1.457622e-29
fitMeasures(fit, c("chisq.scaled", "cfi.robust.mean", "cfi.robust.min", "cfi.robust.max"))
#>    chisq.scaled cfi.robust.mean  cfi.robust.min  cfi.robust.max 
#>       0.5671328       1.0000000       1.0000000       1.0000000
```

The status table has one row per imputation. `on_fail = "omit"` drops a
failed imputation and records the message, as long as two fits remain.
The default is `"stop"`.

``` r

fit@status
#>   no.conv no.se no.npd message
#> 1   FALSE FALSE  FALSE        
#> 2   FALSE FALSE  FALSE
```

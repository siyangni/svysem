# Design-based SEM with svysem

A design-based structural equation model takes a lavaan model and a
survey design. Sampling weights, strata, and primary sampling units
enter the variance. The point estimates are the weighted estimates.

``` r

library(svysem)
design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 1)
```

## Estimators

`estimator = "MLM"` is maximum likelihood with a Satorra-Bentler
adjustment. It is for complete continuous data. A regression with a mean
structure and free covariates matches the design used by
[`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html).

``` r

regression <- svysem(
  "y ~ x1 + x2", design,
  fixed.x = FALSE, meanstructure = TRUE
)
coef(regression)
#>   y~x1   y~x2   y~~y x1~~x1 x1~~x2 x2~~x2    y~1   x1~1   x2~1 
#>  0.441 -0.483  0.880  0.726 -0.085  0.958  0.902  0.088  0.104
```

`estimator = "MLR"` is robust maximum likelihood with a Yuan-Bentler
adjustment. `missing = "ml"` scores each row on the indicators that row
observed. The sampling weight of an incomplete row stays in the design.

``` r

incomplete <- design
incomplete$variables$y2[c(1, 8, 15)] <- NA
mlr <- svycfa(
  "f =~ y1 + y2 + y3 + y4", incomplete,
  estimator = "MLR", missing = "ml"
)
mlr@estimator
#> [1] "MLR"
lavaan::lavInspect(mlr, "nobs")
#> [1] 48
```

`estimator = "WLSMV"` is for ordered indicators. The point estimates use
lavaan’s weighted diagonally weighted least squares matrix. The standard
errors use replicate weights. With strata, the default replicate type is
a stratified jackknife.

``` r

frame <- design$variables
frame$a <- ordered(findInterval(frame$y1, c(-0.5, 0.5)))
frame$b <- ordered(findInterval(frame$y2, c(0)))
frame$c <- ordered(findInterval(frame$y3, c(-0.4, 0.6)))
ordinal <- survey::svydesign(
  ids = ~psu, strata = ~strat, weights = ~w, data = frame
)
ordered_fit <- svycfa(
  "f =~ a + b + c", ordinal,
  estimator = "WLSMV", ordered = c("a", "b", "c")
)
ordered_fit@variance_method
#> [1] "replication"
```

## Domains and replicate-weight designs

[`subset()`](https://rdrr.io/r/base/subset.html) keeps the primary
sampling units of the full design, which is the variance for a domain. A
design that already carries replicate weights is used as it stands.
`variance_method` is then `"replication"` for every estimator.

``` r

domain <- subset(design, strat == 1)
domain_fit <- svysem(
  "y ~ x1", domain, fixed.x = FALSE, meanstructure = TRUE
)
domain_fit@df
#> [1] 5

replicated <- survey::as.svrepdesign(design, type = "auto")
rep_fit <- svycfa("f =~ y1 + y2 + y3 + y4", replicated, estimator = "MLR")
rep_fit@variance_method
#> [1] "replication"
```

## Factor scores

[`predict()`](https://rdrr.io/r/stats/predict.html) returns posterior
means for continuous outcomes.
[`lavaan::lavPredict()`](https://rdrr.io/pkg/lavaan/man/lavPredict.html)
on a weighted fit returns a different number. Draws from the posterior,
when requested, are plausible values for one completed data set.

``` r

predict(mlr)[1, , drop = FALSE]
#>             f
#> 1 -0.05815887
```

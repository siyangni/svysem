# Design-based structural equation models.

`svysem()`, `svycfa()` and `svygrowth()` fit a lavaan model to a
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)
object. Sampling weights, strata and primary sampling units enter the
standard errors through the design-based covariance of the sample
moments (Gamma), which is passed to lavaan as `NACOV` together with the
raw data and the sampling weights. The fit is not built from sample
moments. `svysem()` uses the defaults of
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html), `svycfa()`
those of [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html), and
`svygrowth()` those of
[`lavaan::growth()`](https://rdrr.io/pkg/lavaan/man/growth.html).
Arguments in `...` are passed on to that fitter.

## Usage

``` r
svysem(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

svycfa(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

svygrowth(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# S3 method for class 'survey.design2'
svysem(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# S3 method for class 'survey.design'
svysem(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# S3 method for class 'svyrep.design'
svysem(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# Default S3 method
svysem(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# S3 method for class 'survey.design2'
svycfa(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# S3 method for class 'survey.design'
svycfa(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# S3 method for class 'svyrep.design'
svycfa(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# Default S3 method
svycfa(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# S3 method for class 'survey.design2'
svygrowth(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# S3 method for class 'survey.design'
svygrowth(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# S3 method for class 'svyrep.design'
svygrowth(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)

# Default S3 method
svygrowth(
  model,
  design,
  estimator = c("MLM", "WLSMV", "MLR"),
  replicates = NULL,
  replicate.type = NULL,
  parallel = NULL,
  ...
)
```

## Arguments

- model:

  Model syntax, as for
  [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html),
  [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html) or
  [`lavaan::growth()`](https://rdrr.io/pkg/lavaan/man/growth.html).

- design:

  A `survey.design`, `survey.design2`, or `svyrep.design` object,
  including the result of
  [`subset()`](https://rdrr.io/r/base/subset.html) or
  [`survey::as.svrepdesign()`](https://rdrr.io/pkg/survey/man/as.svrepdesign.html).

- estimator:

  Variance estimator. `"MLM"` (Satorra-Bentler for one group, replicate
  weights when `group=` is set; complete continuous data), `"WLSMV"`
  (ordinal weighted least squares with replicate-weight standard
  errors), or `"MLR"` (Yuan-Bentler, including `missing = "ml"`).

- replicates:

  Number of bootstrap replicate weights. Used when
  `replicate.type = "bootstrap"` for WLSMV or for MLM with `group=`.
  Jackknife and BRR counts are fixed by the primary sampling units.

- replicate.type:

  How
  [`survey::as.svrepdesign()`](https://rdrr.io/pkg/survey/man/as.svrepdesign.html)
  builds weights for WLSMV and for MLM with `group=`. One of `"auto"`,
  `"JK1"`, `"JKn"`, `"BRR"`, or `"bootstrap"`. `NULL` is `"auto"`.

- parallel:

  `NULL`, or a function of `(X, FUN)` with the same contract as
  [`lapply()`](https://rdrr.io/r/base/lapply.html), used to loop WLSMV
  replicates and multigroup MLM replicates. On a `svyimputationList`,
  `parallel` loops imputations instead, and the replicates inside each
  imputation stay sequential.

- ...:

  Passed to [`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html),
  [`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html) or
  [`lavaan::growth()`](https://rdrr.io/pkg/lavaan/man/growth.html).
  `data`, `sampling.weights`, sample moments, `NACOV`, and `WLS.V` are
  taken from the design. `ov.order` is set internally so Gamma matches
  the model.

## Value

An S4 object of class `svysem` that contains the lavaan fit. Slots are
`df` (design degrees of freedom), `estimator` (`"MLM"`, `"WLSMV"`, or
`"MLR"`), `variance_method` (`"linearization"` or `"replication"`) and
`svysem_call`.

## Details

The returned object extends lavaan, so
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`lavaan::lavTestWald()`](https://rdrr.io/pkg/lavaan/man/lavTestWald.html)
and
[`lavaan::fitMeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html)
report design-based results.
[parameterEstimates()](https://siyangni.github.io/svysem/reference/parameterEstimates.md)
and
[standardizedSolution()](https://siyangni.github.io/svysem/reference/standardizedSolution.md)
include the design-based standard errors and confidence intervals.
lavaan's own functions omit those columns unless `se = TRUE`, because
the object's class is `svysem` rather than `lavaan`.
[`summary()`](https://rdrr.io/r/base/summary.html) and the printed
display add the design degrees of freedom and the estimator.

`"MLM"` is maximum likelihood with a Satorra-Bentler adjustment for
complete continuous data. `missing = "ml"` is an error and names MLR as
the estimator that handles it. Rows that are incomplete on the model
variables are not dropped: their sampling weight is set to 0 inside the
design used for Gamma, which is how
[`subset()`](https://rdrr.io/r/base/subset.html) represents a domain,
and the primary sampling unit structure is kept. Zero weights are legal.
A negative weight, or a design with no positive weight, is an error.

lavaan's default `sampling.weights.normalization = "total"` rescales
weights so that they sum to the number of rows. Means and covariances of
products are invariant to that constant, so Gamma is computed from the
design's sampling weights on the same rows. The products use the maximum
likelihood divisor (the weighted residual products, with no
`n / (n - 1)` factor).
[`survey::svyvar()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
is not used. Means and second-moment products are stacked in one
linearization, so the covariance between them is kept; that term is what
makes a regression intercept match
[`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html).

Gamma follows lavaan's model order (`ov.order = "force.model"` on the
fit that receives `NACOV`). Without that argument lavaan reorders
variables to the data order and pairs Gamma with the wrong columns.
`fixed.x = TRUE`, the default for
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html), is supported
when those columns are the observed second moments, plus the means when
a mean structure is estimated. `conditional.x = TRUE` is not supported
for MLM or MLR. MLR fits one group of continuous variables. Single-group
MLM does too, by linearization. `group=` on MLM refits the weighted
model on replicate weights, because lavaan keeps one `NACOV` per group
and that drops the covariance between groups. The multigroup chi-square
is the maximum-likelihood statistic, not a Satorra-Bentler statistic.
`fitMeasures()` does not report a scaled chi-square for that fit: the
scaled statistic lavaan adds from the sampling weights ignores primary
sampling units, so it is dropped. Clustering is taken from the design;
do not pass `cluster=`, which would double-count primary sampling units.
`survey.lonely.psu` is read by
[`survey::svymean()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
and is never set here.

`"MLR"` is robust maximum likelihood with a Yuan-Bentler adjustment. The
variance is a sandwich: lavaan's observed information, and a meat from
the design variance of the casewise scores
([`survey::svytotal()`](https://rdrr.io/pkg/survey/man/surveysummary.html)).
Equality constraints, including a shared label, `==`, and a nonlinear
`==`, are included. The sandwich is solved on lavaan's constraint basis
and expanded back to the free parameters. `std.lv` and `effect.coding`
are accepted. `missing = "ml"` scores the observed indicators; it does
not fall back to complete cases. On a `survey.design` the variance
method is linearization. On a `svyrep.design` it is replication, because
[`survey::svytotal()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
uses the replicate weights.

`"WLSMV"` is lavaan's weighted least squares with the diagonally
weighted least squares matrix taken from the weighted sample. Point
estimates match that weighted fit. Standard errors and the scaled test
use the replicate-weight covariance of lavaan's sample statistics
(thresholds, continuous means and variances, and correlations), with the
scale from
[`survey::as.svrepdesign()`](https://rdrr.io/pkg/survey/man/as.svrepdesign.html)
and [`survey::svrVar()`](https://rdrr.io/pkg/survey/man/svrVar.html).
`replicate.type = "auto"` is a stratified jackknife when the design has
strata and JK1 otherwise. `replicates` is only the bootstrap count.
Ordered indicators, `group=`, and `missing = "pairwise"` or `"listwise"`
are accepted. lavaan stores one `NACOV` per group, so the covariance of
sample statistics between groups is not passed on. The weight matrix
Mplus uses for `TYPE = COMPLEX` has not been compared with lavaan's. For
WLSMV, MLR, and multigroup MLM, `se` and `test` are set by the design
and are not taken from `...`.

A `svyrep.design` is accepted for every estimator. `variance_method` is
then `"replication"`. `replicates` and `replicate.type` are an error on
a design that already has replicate weights.

[`predict()`](https://rdrr.io/r/stats/predict.html) returns factor
scores for continuous outcomes. See
[`predict.svysem()`](https://siyangni.github.io/svysem/reference/predict.svysem.md).

## References

Oberski D (2014). "lavaan.survey: An R Package for Complex Survey
Analysis of Structural Equation Models." Journal of Statistical
Software, 57(1), 1–27.
[doi:10.18637/jss.v057.i01](https://doi.org/10.18637/jss.v057.i01)

Yuan KH, Bentler PM (2000). "Three Likelihood-Based Methods for Mean and
Covariance Structure Analysis with Nonnormal Missing Data." Sociological
Methodology, 30(1), 165–200.
[doi:10.1111/0081-1750.00078](https://doi.org/10.1111/0081-1750.00078)

## See also

[svysem](https://siyangni.github.io/svysem/reference/svysem-class.md),
[`predict.svysem()`](https://siyangni.github.io/svysem/reference/predict.svysem.md),
[`simulate_design()`](https://siyangni.github.io/svysem/reference/simulate_design.md),
[`lavaan::sem()`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html),
[`survey::svyglm()`](https://rdrr.io/pkg/survey/man/svyglm.html)

## Examples

``` r
design <- simulate_design(n_psu = 6, n_per = 4, n_strata = 2, seed = 1)
fit <- svysem(
  "y ~ x1 + x2", design,
  fixed.x = FALSE, meanstructure = TRUE
)
fit
#> svysem fit
#> Estimator: MLM
#> Variance method: linearization
#> Design df: 10
#> lavaan 0.7-3 ended normally after 1 iteration
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         9
#> 
#>   Number of observations                                  48
#>   Sampling weights variable                   .svysem_weight
#> 
#> Model Test User Model:
#>                                               Standard      Scaled
#>   Test Statistic                                 0.000       0.000
#>   Degrees of freedom                                 0           0
coef(fit)
#>   y~x1   y~x2   y~~y x1~~x1 x1~~x2 x2~~x2    y~1   x1~1   x2~1 
#>  0.441 -0.483  0.880  0.726 -0.085  0.958  0.902  0.088  0.104 
vcov(fit)
#>          y~x1   y~x2   y~~y x1~~x1 x1~~x2 x2~~x2    y~1   x1~1   x2~1
#> y~x1    0.048                                                        
#> y~x2    0.007  0.020                                                 
#> y~~y   -0.011 -0.008  0.030                                          
#> x1~~x1  0.022  0.004 -0.004  0.039                                   
#> x1~~x2  0.000 -0.001 -0.001 -0.006  0.005                            
#> x2~~x2 -0.016  0.002  0.011 -0.004 -0.005  0.021                     
#> y~1    -0.006 -0.010  0.007  0.002  0.001 -0.006  0.021              
#> x1~1    0.020  0.000 -0.005  0.009  0.001 -0.010  0.004  0.015       
#> x2~1   -0.026 -0.008 -0.009  0.001  0.000  0.004  0.001 -0.014  0.046
summary(fit)
#> Estimator: MLM
#> Variance method: linearization
#> Design degrees of freedom: 10
#> lavaan 0.7-3 ended normally after 1 iteration
#> 
#>   Estimator                                         ML
#>   Optimization method                           NLMINB
#>   Number of model parameters                         9
#> 
#>   Number of observations                                  48
#>   Sampling weights variable                   .svysem_weight
#> 
#> Model Test User Model:
#>                                               Standard      Scaled
#>   Test Statistic                                 0.000       0.000
#>   Degrees of freedom                                 0           0
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                           Robust.sem
#>   Information                                 Expected
#>   Information saturated (h1) model          Structured
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   y ~                                                 
#>     x1                0.441    0.218    2.020    0.043
#>     x2               -0.483    0.143   -3.384    0.001
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   x1 ~~                                               
#>     x2               -0.085    0.070   -1.210    0.226
#> 
#> Intercepts:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .y                 0.902    0.144    6.270    0.000
#>     x1                0.088    0.123    0.715    0.474
#>     x2                0.104    0.214    0.485    0.627
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .y                 0.880    0.173    5.076    0.000
#>     x1                0.726    0.197    3.690    0.000
#>     x2                0.958    0.145    6.595    0.000
#> 
cfa_fit <- svycfa("f =~ y1 + y2 + y3 + y4", design)
growth_fit <- svygrowth(
  "i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
   s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4",
  design
)
coef(cfa_fit)
#>  f=~y2  f=~y3  f=~y4 y1~~y1 y2~~y2 y3~~y3 y4~~y4   f~~f 
#>  1.044  1.071  1.219  0.334  0.118  0.216  0.383  1.053 
coef(growth_fit)
#> y1~~y1 y2~~y2 y3~~y3 y4~~y4   i~~i   s~~s   i~~s    i~1    s~1 
#>  0.266  0.117  0.205  0.295  1.067  0.032  0.031  0.249 -0.039 
predict(cfa_fit)[1, , drop = FALSE]
#>           f
#> 1 0.1295933
mlr_fit <- svycfa("f =~ y1 + y2 + y3 + y4", design, estimator = "MLR")
coef(mlr_fit)
#>  f=~y2  f=~y3  f=~y4 y1~~y1 y2~~y2 y3~~y3 y4~~y4   f~~f 
#>  1.044  1.071  1.219  0.334  0.118  0.216  0.383  1.053 
```

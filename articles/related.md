# How svysem relates to Mplus, lavaan and lavaan.survey

svysem follows the design-based procedure in Oberski (2014,
<doi:10.18637/jss.v057.i01>): a lavaan model, and a variance that uses
the survey design. No code is taken from lavaan.survey. The agreements
that the package tests are with `survey`, with stored lavaan.survey
results, and with replicate weights from the same design. Agreement with
Mplus output has not been checked. In particular, the weight matrix
Mplus uses for `TYPE = COMPLEX` with WLSMV has not been compared with
lavaan’s diagonally weighted least squares matrix.

``` r

library(svysem)
design <- simulate_design(n_psu = 6, n_per = 5, n_strata = 2, seed = 4)
model <- "f =~ y1 + y2 + y3 + y4"
fit <- svycfa(model, design)
weighted <- lavaan::cfa(model, data = design$variables, sampling.weights = "w")
max(abs(coef(fit) - lavaan::coef(weighted)[names(coef(fit))]))
#> [1] 1.332268e-15
```

The coefficients above agree because both fits use the sampling weights.
The standard errors do not. lavaan treats the rows as independent.
svysem also uses the strata and the primary sampling units.

lavaan rescales the sampling weights so that they sum to the number of
rows. Means and maximum-likelihood covariances are unchanged by that
rescaling, so the same rows are passed to the design with their original
weights.

`survey.lonely.psu` is read when the design has a stratum with one
primary sampling unit. The fit does not change that option.

MLR fits one group. Single-group MLM is a linearization. With `group=`,
MLM uses replicate weights: lavaan keeps one covariance matrix of sample
statistics per group, and that matrix drops the covariance between
groups. The omission is large enough to move a between-group contrast
when primary sampling units contain both groups, so the replicate
standard errors are the ones returned. The multigroup MLM chi-square
remains the maximum-likelihood statistic. WLSMV accepts `group=` and
still passes one matrix per group. Clustering belongs in the design.
Passing `cluster=` would count the primary sampling units twice, and
svysem rejects it.

MLM uses a Satorra-Bentler adjustment for complete continuous data. MLR
uses a Yuan-Bentler adjustment (Yuan and Bentler 2000,
<doi:10.1111/0081-1750.00078>), including when `missing = "ml"`. The MLR
variance is a sandwich: lavaan’s observed information, and the design
variance of the casewise scores. WLSMV standard errors are the
replicate-weight variance of lavaan’s own sample statistics. With
strata, those replicates are a stratified jackknife: one primary
sampling unit is removed from its stratum at a time. For a nonlinear
statistic from a stratified multistage sample in which the primary
sampling units are drawn with replacement, that jackknife and balanced
repeated replication are consistent variance estimators (Krewski and Rao
1981, <doi:10.1214/aos/1176345580>).

For a multiply imputed design, each imputation is fitted and the
parameters are pooled with Rubin’s rules. The within-imputation
variances are the design-based ones. When the imputation model and the
analysis design do not match, that pooled variance can be wrong. Kott
(1995) describes this for weighted survey data, and Kim, Brick, Fuller
and Kalton (2006, <doi:10.1111/j.1467-9868.2006.00546.x>) show that the
bias can be large for a domain mean when many values are imputed. svysem
does not correct for it.


<!-- README.md is generated from README.Rmd. Please edit that file -->

# svysem

Design-based structural equation models for complex survey data. The
model is fitted with [lavaan](https://lavaan.ugent.be/). Sampling
weights, clustering and stratification come from a
[survey](https://cran.r-project.org/package=survey) design, the
procedure of Mplus `TYPE = COMPLEX`.

Continuous outcomes use MLM (Satorra-Bentler) or MLR (Yuan-Bentler),
including `missing = "ml"`. A model with `group=` uses replicate weights
for MLM, so a between-group contrast keeps the covariance from primary
sampling units that contain both groups. Ordered indicators use WLSMV,
with replicate-weight standard errors. `predict()` returns factor scores
for continuous outcomes. A `svyrep.design` is accepted.

Completed data sets from a multiple imputation are fitted one by one and
pooled with Rubin’s rules. The design degrees of freedom are the
complete-data degrees of freedom.

The method for a single design follows Oberski (2014,
<doi:10.18637/jss.v057.i01>). No code is taken from lavaan.survey.

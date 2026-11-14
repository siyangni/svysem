# svysem

Design-based structural equation models for complex survey data. The model is fitted with [lavaan](https://lavaan.ugent.be/). Sampling weights, clustering and stratification come from a [survey](https://cran.r-project.org/package=survey) design, the procedure of Mplus `TYPE = COMPLEX`.

Completed data sets from a multiple imputation are fitted one by one and pooled with Rubin's rules. The design degrees of freedom are the complete-data degrees of freedom.

The method for a single design follows Oberski (2014, <doi:10.18637/jss.v057.i01>). No code is taken from lavaan.survey.

This repository is the package specified in `SVYSEM_PLAN.md`.

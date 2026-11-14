# svysem: design-based structural equation models for complex survey data

Revised 2026-11-13 (fourth version, after a verification review). Replaces all earlier versions of this file. Section "What the review changed" lists the corrections and the evidence for each.

## Goal and scope

One CRAN-ready R package, `svysem`, at `/home/siyang/work/svysem`, public at `github.com/siyangni/svysem`. It has one methodological scope: **design-based estimation of structural equation models fitted with lavaan**. Design-based means sampling weights, clustering (PSU) and stratification, the procedure of Mplus `TYPE = COMPLEX`. The scope covers:
- design-based standard errors and test statistics for MLM, WLSMV and MLR (with FIML);
- subpopulations and replicate-weight designs;
- factor scores that are correct for weighted fits;
- fitting each completed data set of a multiple imputation and pooling the SEM results.

**Scope rule.** An exported function belongs in svysem only if (a) it is part of design-based SEM estimation and (b) no maintained CRAN package provides it for design-based fits. Design objects, linearization and replicate variances are reused from `survey`. Model fitting is reused from `lavaan`. Everything that is not SEM (GLMs, hazard models, descriptive statistics) stays with `survey` and the private dissertation repo.

## What the review changed (all checked by running code on 2026-11-13)

| # | Finding | Evidence | Change |
|---|---|---|---|
| 1 | The MLM design is sound | `growth(model, data, sampling.weights =, estimator = "MLM", NACOV = design Gamma)` reproduced `lavaan.survey` on a stratified cluster sample: estimates to 1e-15, SEs to 8e-15 relative, identical scaled χ² | Kept. Raw data are passed to lavaan, so no fitting from moments |
| 2 | Extending the lavaan class works | An S4 class with `contains = "lavaan"` passed through `parameterEstimates()`, `lavTestWald()` and `mitml::testEstimates()` | Kept |
| 3 | **Plan error:** the MI "methods" cannot exist | `parameterEstimates`, `standardizedSolution`, `lavTestWald`, `lavTestLRT` and `lavPredict` are plain functions in lavaan, not generics. Only `coef`, `vcov`, `summary`, `nobs`, `fitMeasures`, `predict` are generic | MI results get svysem's own functions `pooled_estimates()` and `pooled_wald()`, plus methods for the true generics |
| 4 | **Plan error:** factor scores from weighted fits are wrong in lavaan 0.6-21 | `lavPredict()` on a fit with `sampling.weights` differed from the posterior mean α + ΨΛ′Σ⁻¹(y − μ) by up to 1.7 (unweighted fit: 5e-15). `newdata` and other methods do not help. Clearing the weights slot of a copy restores 4e-15. `semTools::plausibleValues()` calls `lavPredict()` and inherits the problem | svysem defines a `predict()` method that returns correct scores and plausible values. The old test 8a would have failed. Reported upstream; rechecked on lavaan 0.7 |
| 5 | **The gap is narrower than stated** | `Grexxi/lavaan.survey-modernized` on GitHub (created 2026-05, unofficial fork, v1.2.0, 0 stars, not on CRAN) updates `lavaan.survey` for current lavaan and adds ordinal SEM through replicate weights and MI by pooling moments. `daob/lavaan.survey` had one commit in 2025-10 (a CI workflow) | Value statement rewritten. The fork is cited as prior art |
| 6 | lavaan is changing in this area | The development version is 0.7-3 and has new options `sampling.weights.type = c("design", "frequency")` and `se = "robust.cluster"`, `"robust.cluster.sem"`. It still has no strata argument, and the "categorical + clustered" refusal is still in the source | CI tests against lavaan's development version from phase 1. The slot-injection part (MLR) moves to the last phase |
| 7 | A complete-data shortcut for MLR does not hold | MLM with observed information differed from MLR SEs by 2.2% | MLR always uses the casewise-score sandwich |
| 8 | Importing `mitml` is heavy | `mitml` adds 39 packages beyond lavaan + survey (among them lme4, ordinal, haven, readr, vroom, jomo, pan). `mitools` adds none but its small-sample df differs: 5.985 against 5.824 from `mitml`, `mice::pool.scalar` and Paper 1's `pool_scalar()` on the same input | Rubin's rules, Barnard–Rubin df, D1 and D2 become about 60 lines of internal code, tested against `mitml`, `mice` and `lavaan.mi` as oracles in Suggests. See "Design decision A" |
| 9 | `mitml::testConstraints` cannot take lavaan labels or `:=` names | `"b - c"`: "object 'b' not found". A `:=` name: coercion error. `` `y~x` - `y~z` `` works | Another reason for finding 8's change: `pooled_wald()` evaluates lavaan constraint syntax itself |
| 10 | Writing polychoric influence functions by hand is the most error-prone piece, and lavaan's are not exported | WLSMV accepts a user `NACOV` with raw ordered data (verified: SEs scale as √1.5 when Gamma is multiplied by 1.5). The fork above uses replicate weights for the same purpose | WLSMV Gamma comes from `survey::as.svrepdesign()` + `withReplicates()` over lavaan's own sample statistics. Only public functions are used. Analytic linearization becomes a later option |
| 11 | Multigroup models with PSUs that span groups | lavaan takes `NACOV` as one matrix per group, so covariance between groups' sample statistics is dropped (reasoned, not yet run). `lavaan.survey` has the same property | Phase 1 measures the size of the effect against the replicate-weight oracle. Documented as a limitation of the MLM route if confirmed |
| 12 | Test 1's tolerance was too strict for variances | `svyvar` and a saturated SEM use different denominators (n − 1 vs n) | The test states the convention or compares against a hand-computed influence function |
| 13 | A dissertation inconsistency, outside the package | Paper 3 pools with `mitools::MIcombine`, whose df differs from the Barnard–Rubin df used in Papers 1 and 2 (finding 8) | Listed under decisions. svysem's pooling matches Papers 1 and 2 |

## Does the package address a real gap?

Yes, with a narrower claim than before.

**Demand.** `lavaan.survey` was the standard R tool for this problem and was archived from CRAN on 2023-06-20 ("issues were not corrected despite reminders"). Since then no CRAN package fits SEMs with strata. lavaan itself offers `sampling.weights` and `cluster`, without strata or subpopulations, and refuses clusters with ordered indicators.

**What exists off CRAN.** The unofficial fork `Grexxi/lavaan.survey-modernized` covers continuous MLM, ordinal WLSMV through replicate weights, and MI by pooling moments. It keeps `lavaan.survey`'s interface (fit naively first, then correct).

**What svysem would add that nothing has, on or off CRAN:**
1. A maintained CRAN package with a `survey`-style interface (`svysem(model, design)`), tested against lavaan's development version.
2. MLR with strata, including FIML (the default Mplus analysis under `TYPE = COMPLEX`).
3. Multiple imputation the Mplus way: the model is fitted to each completed data set and the parameters are pooled with Rubin's rules and design degrees of freedom, with pooled standardized estimates, Wald tests in lavaan syntax and fit summaries. The fork and `lavaan.survey` pool moments and fit once, which is a different estimator.
4. Correct factor scores and plausible values for weighted fits (finding 4).
5. Subpopulations handled as domains, and replicate-weight designs for every estimator.

**Position toward `lavaan.survey` (decided 2026-11-13).** The package is not on CRAN, so svysem fills the gap as a new package without waiting on its author's plans. Oberski (2014) is cited for the method, and no code is copied.

**Optional upstream reports** (they help other users; svysem does not depend on them): finding 4 to the lavaan maintainers, finding 9 to the mitml maintainer.

## Existing tools and what is reused

| Need | Provided by | Use in svysem |
|---|---|---|
| Design object, `nest`, lonely PSU option, domains via `subset()` | `survey::svydesign` | Input class `survey.design2` |
| Linearization variance of totals (PSU totals, stratum centering, n_h/(n_h − 1)) | `survey::svytotal` / `svyrecvar` | Applied to casewise contributions that svysem computes |
| Replicate weights and replicate variances | `survey::as.svrepdesign`, `withReplicates` | WLSMV Gamma; `svyrep.design` input; test oracle for every estimator |
| Design degrees of freedom | `survey::degf` | Reference df for tests and pooling |
| Imputed data with a design | `svydesign(data = mitools::imputationList(...))` | Input class `svyimputationList` |
| Model fitting, robust SEs and scaled tests given a Gamma | `lavaan` with `NACOV` | MLM and WLSMV |
| Sample statistics for ordered data (thresholds, polychoric correlations), pairwise deletion | `lavaan` | Recomputed per replicate |
| Fisher z for user-defined correlations | lavaan `:=` with `atanh()` (verified) | Vignette example |
| Invariance model syntax | `semTools::measEq.syntax(return.fit = FALSE)` | Syntax passed to `svycfa()` |
| Rubin's rules, Barnard–Rubin df, D1, D2 | `mitml`, `mice`, `lavaan.mi` | **Test oracles** (Suggests). Formulas are internal: design decision A |

### Design decision A: pooling formulas are internal (final, 2026-11-13)

Earlier versions imported `mitml` so that no formula was re-implemented. The review found three costs: 39 extra packages for every user, an archival risk (a package is archived when a hard dependency is archived, which is how `lavaan.survey` users lost it), and the label failure in finding 9, which needs a translation layer anyway. The pooling formulas are short and published. svysem therefore:
- implements Rubin's rules, Barnard–Rubin df, D1 and D2 in one internal file, not exported as general-purpose pooling functions;
- tests them to 1e-10 against `mitml::testEstimates`, `mitml::testConstraints`, `mice::pool.scalar` and `lavaan.mi::calculate.D2`, all in Suggests;
- exports only SEM-specific functions that have no equivalent for design-based fits.

The decision keeps the package at three non-base dependencies, all of which lavaan and survey already require.

## Name

**`svysem`**:
- It follows the naming of design-based model functions around the `survey` package: `svyglm`, `svycoxph`, `svyolr`, and the package `svylme`. `svysem()` reads as the SEM member of that family.
- It is free: not on CRAN, never archived there, and no GitHub repository uses it (checked 2026-11-13).
- Rejected: `lavaan.complex` / `lavaan.design` (the `lavaan.*` prefix suggests affiliation with the lavaan project); `complexsem` (an unrelated GitHub repository uses it).

## Public API

```r
svysem(model, design, estimator = c("MLM", "WLSMV", "MLR"), ...)   # lavaan::sem defaults
svycfa(model, design, estimator = ..., ...)                        # lavaan::cfa defaults
svygrowth(model, design, estimator = ..., ...)                     # lavaan::growth defaults
```

S3 generics that dispatch on `design`, the convention of `survey::svyglm`.

| `design` class | Returns | Variance method |
|---|---|---|
| `survey.design2` (also after `subset()`) | S4 `svysem`, `contains = "lavaan"` | Taylor linearization (MLM, MLR); replicate weights built from the design (WLSMV) |
| `svyrep.design` | S4 `svysem` | The design's replicate weights, for every estimator |
| `svyimputationList` | S4 `svysem_mi` | Each completed design fitted separately |

**A `svysem` object** is a lavaan object whose covariance matrix and test statistic are design-based. Every lavaan function therefore reports design-based results without further steps: `coef()`, `vcov()`, `summary()`, `parameterEstimates()` (with `:=` rows), `standardizedSolution()`, `lavTestWald()`, `lavTestLRT()`, `fitMeasures()`. Extra slots: design df, estimator, variance method, and the call. One method is added:
- `predict(object, newdata = NULL, draws = 0L)`: posterior means of the latent variables at the design-based estimates. With `draws > 0` it returns that many plausible-value sets drawn from the posterior. It exists because of finding 4. The help page warns that `lavaan::lavPredict()` on a weighted fit does not return these values in lavaan 0.6-21.

**A `svysem_mi` object** holds the list of `svysem` fits, the design df, and a status table per imputation (converged, standard errors available, latent covariance matrix positive definite, error message; column names follow `lavaan.mi`'s `omit.imps` vocabulary). Arguments of the `svyimputationList` method:
- `on_fail = c("stop", "omit")`, default `"stop"`;
- `parallel = NULL`: any lapply-like function, for example `function(X, FUN) parallel::mclapply(X, FUN, mc.cores = 20)`.

Functions for `svysem_mi`:

| Function | Kind | Result |
|---|---|---|
| `summary()`, `print()`, `coef()`, `vcov()`, `nobs()` | methods of true generics | Pooled estimates; total covariance matrix T = Ū + (1 + 1/m)B |
| `pooled_estimates(object, standardized = FALSE, level = 0.95)` | new function | Data frame: est, se, t, df (Barnard–Rubin, complete-data df = design df), p, CI, riv, lambda, fmi, m. With `standardized = TRUE`, each fit's design-based standardized estimates are pooled; correlations are pooled on the Fisher z scale and back-transformed |
| `pooled_wald(object, constraints)` | new function | D1 test for lavaan constraint syntax: labels (`"b2h == b2c"`), several constraints, and nonlinear constraints on labeled or `:=` parameters (delta method within each imputation) |
| `fitMeasures()` | method of lavaan's S4 generic | D2-pooled scaled χ²; mean, minimum and maximum of CFI, TLI, RMSEA and SRMR (lavaan's robust versions). The indices are not passed through Rubin's rules |
| `predict()` | method | Factor scores averaged over imputations, or per imputation |

### Estimators

- **MLM** (default; reproduces the current papers). svysem computes the design-based covariance matrix of the sample means and covariances (Gamma) by linearization and passes it to lavaan as `NACOV`, with the raw data and `sampling.weights`. lavaan computes robust SEs and the Satorra–Bentler scaled test, so `lavTestLRT()` difference tests are design-based too. Complete data only: `missing = "ml"` with MLM is an error that names MLR.
- **WLSMV**. Gamma for thresholds and polychoric correlations is the replicate-weight covariance of lavaan's own sample statistics: `as.svrepdesign(design)` (the type and number of replicates are arguments; default JKn for stratified cluster designs) and `withReplicates()`. Each replicate recomputes the statistics with the replicate weight as `sampling.weights`, in parallel when requested.
  - Supports `missing = "listwise"` and `"pairwise"`, mixed ordered and continuous indicators, and several groups, because lavaan computes the statistics.
  - The DWLS weight matrix stays lavaan's default from the weighted data, so point estimates equal lavaan's weighted WLSMV estimates. The design changes SEs and tests only. The help page documents this choice and that Mplus's choice of weight matrix under `TYPE = COMPLEX` has not been compared.
- **MLR**. Observed-information bread, and a meat from linearization of casewise scores. Scores come from `lavScores()` on complete data, and from numerical differentiation of `lavInspect(fit, "loglik.casewise")` under weighted FIML (lavaan's `lavScores()` errors on weighted FIML fits). The covariance matrix, the `:=` standard errors and the Yuan–Bentler scaling are then written into the object's slots. This is the Mplus `TYPE = COMPLEX; ESTIMATOR = MLR` procedure and the only estimator that needs slot injection.

Implementation rules:
- No `:::` calls.
- Input checks at the top of every exported function with informative, classed errors: design class, variables of the model present in the design, weights positive, estimator compatible with the data and the `missing` setting.
- The `survey.lonely.psu` option is read, never set.

## Internal files

```
R/
  svysem-package.R    "_PACKAGE" documentation, imports
  generics.R          svysem(), svycfa(), svygrowth() and their methods
  class.R             S4 classes svysem and svysem_mi; show, summary
  gamma_linear.R      casewise moment contributions -> linearization -> NACOV (MLM)
  gamma_replicate.R   replicate-weight Gamma of lavaan sample statistics (WLSMV, svyrep designs)
  sandwich_mlr.R      casewise scores (analytic, numerical), observed information, slot injection
  predict.R           posterior means and plausible values
  mi.R                svyimputationList method, status table, parallel
  pooled.R            pooled_estimates(), pooled_wald(), fitMeasures method
  constraints.R       parse lavaan constraint syntax; expand labels and := definitions
  rubin.R             internal: Rubin's rules, Barnard–Rubin df, D1, D2
  checks.R            input validation and classed conditions
  simulate.R          internal DGPs for tests and examples
data-raw/             fixture generation (Rbuildignored)
tests/testthat/fixtures/  stored oracle values (.rds)
```

## Tests that define correctness

All data are simulated. `lavaan.survey` is archived and cannot be in Suggests, so its results are stored fixtures made by `data-raw/make_fixtures.R` with the archived version 1.1.3.1. The script records package and R versions.

1. **Saturated models vs `survey`.** A regression through `svysem()` with `fixed.x = FALSE` matches `svyglm` SEs to 1e-6 relative; a saturated mean model matches `svymean`. Variances are compared with a hand-computed influence function, because `svyvar` uses n − 1. Cases: no strata; strata; `nest = TRUE` with reused PSU labels; each `survey.lonely.psu` setting; a domain via `subset()` (its SEs differ from a fit on the subset rows alone).
2. **Weight scale invariance.** Multiplying every weight by a constant leaves SEs and tests unchanged to 1e-10.
3. **MLM vs the `lavaan.survey` fixture.** Composite piecewise growth, bivariate growth, RI-CLPM. `coef` to 1e-6; `vcov`, `:=` SEs, standardized SEs and the scaled χ² to 1e-6 relative (the review obtained 1e-14).
4. **Replicate-weight oracle, every estimator.** Linearization SEs (MLM, MLR, FIML) agree with `withReplicates()` on the JKn version of the same design within 2% on a design with 60 or more PSUs. This is the only check available for MLR and FIML that needs no archived or commercial software.
5. **Multigroup.** A two-group model with PSUs spanning both groups: the MLM route is compared with the replicate oracle, and the size of any gap in between-group contrasts is recorded (finding 11). The result decides the documentation and whether multigroup MLM switches to the replicate route.
6. **WLSMV.** With PSU = case and no strata, replicate SEs agree with lavaan's weighted WLSMV SEs within 2%. With clusters, SEs and the scaled test change. A configural → metric → scalar sequence from `measEq.syntax(return.fit = FALSE)` fits through `svycfa()` and `lavTestLRT()` runs. Both `missing` settings are covered.
7. **MLR.** The vcov equals a hand-built `solve(H) %*% meat %*% solve(H)`. Numerical scores match `lavScores()` to 1e-6 on complete data. After injection, `parameterEstimates()`, `standardizedSolution()`, `lavTestWald()` and `fitMeasures()` report the design-based values, including `:=` rows.
8. **FIML.** n is unchanged (no listwise deletion); the vcov is finite and symmetric; the replicate oracle of test 4 holds. If scores are unstable, `missing = "ml"` errors with the reason and never falls back to complete cases.
9. **Factor scores.** `predict()` equals α + ΨΛ′Σ⁻¹(y − μ) at the estimates to 1e-8, with `newdata` and under FIML (observed indicators per case). Plausible values have the posterior mean and covariance Ψ − ΨΛ′Σ⁻¹ΛΨ within simulation error. A test documents that `lavPredict()` on the same object does not match; when lavaan fixes this, the test fails and prompts removal of the warning.
10. **Pooling.** `pooled_estimates()` equals `mitml::testEstimates(fits, df.com =)` and `mice::pool.scalar()` to 1e-10 (est, se, df, riv, fmi). `pooled_wald()` equals `mitml::testConstraints` on the same constraints written with lhs-op-rhs names, and handles labels and `:=`. D2 equals `lavaan.mi::calculate.D2`. Pooled standardized estimates equal a hand computation on the Fisher z scale.
11. **MI behavior.** Parallel and sequential fits are identical. One failing imputation errors under `"stop"` and is dropped and recorded under `"omit"`.
12. **Input validation.** Each classed error is triggered once (snapshot tests of the messages).
13. **Hygiene.** No file contains `dissertation`, `mcs` or `MCS`. `R CMD check --as-cran`: 0 errors, 0 warnings, no notes other than "New submission". At most 2 cores in examples, tests and vignettes. Examples under 5 seconds. Tests over 10 seconds use `skip_on_cran()`.

## Acceptance gates (private dissertation repo, `experiments/svysem_gates/`)

Each gate reruns the current analyses through svysem with MLM and compares with the locked CSVs: point estimates to 1e-6, SEs and df to 1e-4 relative.

| Gate | Reruns | Compares with |
|---|---|---|
| P1 | `shared/10` blocks A–C, `shared/13` SAT | `p1_mi_pooled.csv`, `p1_mi_tests.csv`, `p1_mi_fit.csv`, `p1_cc.csv` |
| P2 | `shared/20`, `paper2_parenting/14, 16` | `p2_mi_pooled.csv`, `p2_mi_tests.csv`, `p2_mi_fit.csv`, `domain_riclpm_*.csv`, `*_std_ci.csv` |
| P3 | `paper3/33` A and B, `21`, `R/survey_factor_scores.R`; later `47`, `53` | `onestep_sem_inc17*.csv`, `plausible_values_inc17*.csv`, `childhood_growth_*.csv`; later `model5_factor.csv`, `model5_trait_score.csv` |

Expected differences, reported by the gate and not treated as failures:
- `fmi`: Paper 1's `pool_scalar()` labels λ as `fmi`. svysem reports both `lambda` and `fmi`, so the old column is matched against `lambda`.
- Plausible values use a different random number stream from `survey_plausible_values()`; the gate compares posterior moments and the pooled coefficients within Monte Carlo error.

After the gates, the code that stays in the dissertation repo for Papers 1 and 2 is model syntax, data preparation, the fallback ladder (a loop over the status table) and the proportionality test. Outside svysem by the scope rule:
- Paper 1 blocks D and E (`svymean`, `svyvar`);
- Paper 3's hazard models, predicted risks, attrition analyses, fixed-effects check, regression calibration and the bootstrap in `34`.

## Decisions that change reported numbers (user decides; none is made here)

1. **MLM → MLR** for SEM in Papers 1–3.
2. **Paper 1 measurement models: weights only → full-design WLSMV**, after test 6 passes.
3. **Paper 3 trait models (`47`, `53`): weights only → full design.** Loadings and trait scores stay the same; CFI, RMSEA and SRMR change.
4. **Paper 3 pooling df (deferred by the user, 2026-11-13; to be handled later).** `33` and `41` pool with `mitools::MIcombine`, whose small-sample df differs from the Barnard–Rubin df of Papers 1 and 2 (finding 13). Switching to one rule changes confidence limits slightly where t-based intervals are used.
5. **Paper 1 proportionality test:** keep the minimum-distance test as analysis code, or replace it with `pooled_wald()` on the 3 nonlinear constraints (asymptotically equivalent, numerically different).
6. **FIML instead of MI:** available after phase 5, not proposed.

Each decision is a rerun plus a NUMBERS.md re-lock under the verify-every-claim rule.

## CRAN and current package standards

Part of the definition of done for 0.1.0. *Writing R Extensions* and the *CRAN Repository Policy* are re-read before submission.

**DESCRIPTION**
- `Title: Design-Based Structural Equation Models for Complex Surveys` (title case, under 65 characters, no period).
- `Description:` one paragraph, not starting with the package name; software names in single quotes ('lavaan', 'survey', 'Mplus'); method references as `Author (year) <doi:...>`, each DOI checked against the publisher record.
- `Authors@R:` `person("Siyang", "Ni", email = ..., role = c("aut", "cre", "cph"), comment = c(ORCID = ...))`. Oberski's method is cited; no code from `lavaan.survey` or the fork is copied. If any code is adapted later, its authors are added with role `ctb` or `cph`, as the GPL and CRAN policy require.
- `License: GPL (>= 2)`, no LICENSE file.
- `URL`, `BugReports`; `Depends: R (>= 4.1.0)`.
- `Imports: lavaan (>= 0.6-21), survey (>= 4.5), methods, stats, numDeriv`. Three non-base packages; `numDeriv` is already required by lavaan and survey.
- `Suggests: mitools, mitml, mice, lavaan.mi, semTools, parallel, testthat (>= 3.0.0), knitr, rmarkdown`, each used conditionally (`requireNamespace()`, `skip_if_not_installed()`).
- `VignetteBuilder: knitr`; `Config/testthat/edition: 3`; `Roxygen: list(markdown = TRUE)`; `Encoding: UTF-8`; `Language: en-US`.

**Code and documentation**
- roxygen2 with markdown. Every exported function, method and class has a title, description, `@param` for every argument, `@return` describing the class and its components, `@references`, `@seealso`, and `@examples` that run in under 5 seconds on simulated data. No `\dontrun{}`.
- S4 classes documented with `@slot`. A package-level help page.
- Errors, warnings and messages only through `stop()`, `warning()`, `message()` with condition classes (`svysem_error_*`), so users and tests can catch them. No `print()` or `cat()` outside print and show methods.
- No change to `options()`, `par()` or the working directory. No files written outside `tempdir()`. Seeds are set only when the user passes one.
- testthat 3rd edition; snapshot tests for printed output; `skip_on_cran()` for slow tests; coverage of 80% or more.
- Formatting and linting with `styler`/`lintr` defaults. Functions under 50 lines where practical.

**Vignettes** (knitr, each under 1 minute, simulated data)
1. *Design-based SEM with svysem*: `svydesign` → `svysem()`; estimators; domains; replicate designs; factor scores.
2. *Multiply imputed survey data*: `imputationList` → `svydesign` → `svysem()` → `pooled_estimates()`, `pooled_wald()`, `fitMeasures()`; Fisher z through `:=`; the status table.
3. *How svysem relates to Mplus, lavaan and lavaan.survey*: what each estimator computes; lavaan's weight normalization; lonely PSUs; multigroup models and clusters; the Kott (1995) / Kim et al. (2006) caveat on Rubin's variance with weights. It states which agreements are tested (`lavaan.survey`, `survey`, replicate weights) and that agreement with Mplus output is not yet tested.

**Other files**: `README.Rmd`; `NEWS.md`; `inst/CITATION`; `inst/WORDLIST`; `cran-comments.md`; `.Rbuildignore`; `_pkgdown.yml`; GitHub Actions from `r-lib/actions` v2 (`check-standard` on Linux, macOS and Windows with R release, devel and oldrel; a job with lavaan's development version; `test-coverage`; `pkgdown`).

**Pre-submission checklist**: `devtools::check(cran = TRUE)`; `check_win_devel()`, `check_win_release()`, `check_mac_release()`; `rhub::rhub_check()`; `urlchecker::url_check()`; `spelling::spell_check_package()`; check time under 10 minutes.

## Phases and versions

Ordered so that the dissertation is served first and the riskiest piece comes last.

- **Phase 0, serial.** Package skeleton with `usethis`. `git init -b main`, `gh repo create svysem --public`. First tasks: rerun the finding 4 check on lavaan's development version and record the result; draft the upstream reports for the user to send.
- **Phase 1, two agents in parallel worktrees.**
  - *MLM agent*: `gamma_linear.R`, `class.R`, `generics.R`, `checks.R`, `simulate.R`, fixtures. Tests 1–5, 12.
  - *MI agent*: `rubin.R`, `constraints.R`, `pooled.R`, `mi.R` against a stub fit. Tests 10–11.
  - Merge. Run **Gates P1 and P2**.
- **Phase 2, serial.** `predict.R` (test 9). Run **Gate P3** for `33`, `21` and the factor scores. Help pages, vignettes 1–2, README. Release **0.1.0** on GitHub: MLM, MI, factor scores.
- **Phase 3, serial.** `gamma_replicate.R`: WLSMV and `svyrep.design` (test 6). Vignette 3. Release **0.2.0**. The first CRAN submission is possible from here, as a user decision.
- **Phase 4, serial.** `sandwich_mlr.R`: MLR on complete data, then FIML (tests 7–8). Before coding, check whether lavaan 0.7 offers a public way to supply a user meat or information matrix; if it does, use it and skip slot injection. Gate P3 for `47`, `53`. Release **0.3.0**.

## References (each checked against the full text before it enters the package)

Asparouhov & Muthén (2005); Barnard & Rubin (1999); Binder (1983); Kim, Brick, Fuller & Kalton (2006); Kott (1995); Li, Raghunathan & Rubin (1991); Li, Meng, Raghunathan & Rubin (1991); Lumley (2010); Muthén & Muthén, Mplus User's Guide ch. 9; Oberski (2014); Rosseel (2012); Rubin (1987); Satorra & Bentler (1994); Yuan & Bentler (2000); a reference for replicate-weight variance estimation under stratified cluster designs, to be selected and verified.

## Non-goals

- General-purpose pooling functions (the internal formulas are not exported).
- Wrappers for `svyglm`, conditional logit or predicted risks.
- A design specification object.
- Imputation; mixture models; logit or probit outcomes inside SEM; a new likelihood.
- Mplus, Mplus license files or UK Data Service files in the repo.

## Risks

- **lavaan 0.7.** A major version is in development and touches weights and clusters. MLM and WLSMV rely only on documented arguments (`NACOV`, `sampling.weights`). MLR relies on slot injection and is scheduled last, with a CI job on lavaan's development version from phase 1.
- **Replicate-weight cost for WLSMV.** Each replicate recomputes polychoric correlations. For a model with 50 ordered items and several hundred PSUs this can take a long time. The number and type of replicates are arguments, the loop runs in parallel, and phase 3 starts by timing a model of Paper 1's size.
- **Overlap with the unofficial fork.** If it reaches CRAN first, svysem's distinct contributions are items 2–5 of the gap section. The README states the relation plainly.
- **Multigroup MLM** (finding 11) until test 5 quantifies it.
- **Numerical FIML scores** near boundary solutions; checked by the replicate oracle.
- **Fixture provenance** recorded by `data-raw/make_fixtures.R`.
- **Time.** The package is a separate project from the dissertation. Phases 1–2 serve Papers 1–3's current analyses. Phases 3–4 serve optional decisions 2 and 3 and other users.

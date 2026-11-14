# Build svysem from the revised plan

Source of truth: `/home/siyang/work/svysem/SVYSEM_PLAN.md` (fourth version, 2026-11-13). This file records the feasibility check and the corrections the build follows. It does not replace that plan.

## Where the work happens

Every file an agent creates, edits, or commits is under `/home/siyang/work/svysem`. That includes the package, tests, fixtures, `data-raw/`, and any git worktree. Worktrees live at `/home/siyang/work/svysem/.worktrees/` and are listed in `.gitignore`. Nothing is written under `/home/siyang/work/dissertation`, under other directories in `/home/siyang/work`, or into the home directory.

Agents do not read or modify the dissertation repository. The specification they implement is `SVYSEM_PLAN.md` plus the corrections in this file. Dissertation gates, locked CSVs, and `experiments/svysem_gates/` are out of this build. They are not a later step of this session.

## Feasibility

The revised plan is feasible. The statistical design is the right one: `survey` for the design, `lavaan` for the model, linearization or replicate weights for the variance, and fit-then-pool for imputations. The scope rule (no GLM wrappers, no imputation, no mixtures) is what keeps the package small enough to finish.

Checked on this machine, lavaan 0.6-21 and survey 4.5, on 2026-11-14:

- `parameterEstimates`, `standardizedSolution`, `lavTestWald`, `lavTestLRT` and `lavPredict` are not generics. Finding 3 is correct. `pooled_estimates()` and `pooled_wald()` are required.
- `fitMeasures` is an S4 generic inside lavaan. `coef`, `vcov`, `summary`, `predict`, `anova`, `nobs` and `logLik` are generics. `nobs` is not defined in the lavaan namespace; a method is set on `stats::nobs`.
- `NACOV` is an argument of `lavaan()`. It is not a name in `lavOptions()`. Finding 1 can be implemented.
- `numDeriv` is already imported by lavaan and survey. The three-dependency claim holds.
- `gh search repos svysem` returned no public repository.

No blocking error. The items below are corrections the agents follow so the first implementation does not miss the fixtures.

## Corrections the build follows

1. **No invented ORCID.** `Authors@R` uses Siyang Ni, `johnni.nj@gmail.com`, roles `aut`, `cre`, `cph`, matching `surveymixr`. An ORCID is added only when supplied.
2. **Weight normalization is part of Gamma, not a vignette aside.** lavaan exposes `sampling.weights.normalization`. The MLM Gamma is computed from the weights lavaan uses after that normalization. Raw `svydesign` weights are not assumed to match. Test 3 (the lavaan.survey fixture) is the check.
3. **Test 2 is about means and model parameters, not totals.** Multiplying all weights by a constant leaves a `svymean` and a regression coefficient's standard error unchanged after normalization. It does not leave `svytotal` unchanged. The test must not use a total.
4. **Register both design classes.** S3 methods for `survey.design2` and `survey.design`. The returned object is the S4 class `svysem`. The function and the class may share the name.
5. **`fitMeasures` for `svysem_mi` uses lavaan's existing generic.** The method signature is copied from `getGeneric("fitMeasures")`. No new generic.
6. **`lavaan.survey` is not a dependency.** Fixtures are built by `data-raw/make_fixtures.R` with the installed archived version 1.1.3.1 and stored under `tests/testthat/fixtures/`. The script records package versions.
7. **The package never sets `survey.lonely.psu`, `options()`, or the working directory.**
8. **No dissertation work.** Do not create `experiments/svysem_gates/`, do not open the dissertation repo, and do not compare against its CSVs in this build. Correctness is the tests inside `/home/siyang/work/svysem`.
9. **Phase 4 stays last.** MLR slot injection is not started in this wave. MLM and the pooling functions are.

## What the agents build now

**Step 0, serial, in `/home/siyang/work/svysem`.** Skeleton only: `DESCRIPTION`, `NAMESPACE`, `R/svysem-package.R`, `NEWS.md`, `README.md`, `.gitignore`, `.Rbuildignore`, `tests/testthat.R`. `git init -b main`, one commit, `gh repo create svysem --public --source=. --remote=origin --push`. Description: "Design-based structural equation models: Mplus TYPE=COMPLEX and Rubin's rules for imputed data sets." Do not overwrite `SVYSEM_PLAN.md`.

**Step 1, two agents, separate git worktrees inside the package.** Worktrees are `/home/siyang/work/svysem/.worktrees/mlm` and `/home/siyang/work/svysem/.worktrees/mi`. No worktree is created anywhere else. Neither agent edits `DESCRIPTION` except by appending a roxygen import that the merge will reconcile. Neither runs `devtools::document()` on the shared checkout. Both agents are started with working directory `/home/siyang/work/svysem`.

- **MLM agent** owns `R/gamma_linear.R`, `R/class.R`, `R/generics.R`, `R/checks.R`, `R/simulate.R`, `data-raw/make_fixtures.R`, and tests 1, 2, 3 and 12 from the revised plan. `svysem()` for class `survey.design2`, estimator `"MLM"` only. Other estimators error with class `svysem_error_estimator`. Complete data only: `missing = "ml"` errors with class `svysem_error_missing` and names MLR as the later estimator. Gamma goes to lavaan as `NACOV` together with the raw data and `sampling.weights`. Return value is an S4 object containing the lavaan object, so `coef`, `vcov`, `parameterEstimates`, `lavTestWald` and `fitMeasures` work without new methods.
- **MI agent** owns `R/rubin.R`, `R/constraints.R`, `R/pooled.R`, `R/mi.R`, and tests 10 and 11. Pooling accepts a list of objects that have `coef()` and `vcov()`, plus `df.com`. It does not require the MLM fitter. `pooled_estimates()` and `pooled_wald()` match `mice::pool.scalar` and, when mitml is installed, `mitml::testEstimates` / `testConstraints`, to 1e-8 or tighter. Barnard–Rubin df uses the supplied design df. Labels and `:=` constraints are parsed by svysem, not handed to mitml as text. The `svyimputationList` method calls `svysem()` and is tested with a local stub if the MLM worktree is not merged.

**Step 2, serial, after both agents finish.** Merge the two worktrees back into `/home/siyang/work/svysem`. `roxygen2::roxygenise()`. `testthat::test_local()`. Fix until the owned tests pass. Commit and push `main`. Remove `.worktrees/` when the merge is done. Do not start WLSMV, MLR, factor scores, or vignettes in this wave. Do not leave the package directory.

## Agent rules

- The working directory and every path an agent writes is under `/home/siyang/work/svysem`.
- Simulated data only. No path containing dissertation, MCS, or UKDS. No reads of the dissertation repository.
- No `:::`.
- Errors through `stop()` with classes `svysem_error_*`.
- Examples and tests use at most 2 cores.
- A failing comparison against `lavaan.survey` or `svyglm` is fixed in the Gamma code. The tolerance is not widened to hide it.

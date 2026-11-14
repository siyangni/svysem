# CRAN comments for svysem 0.0.1

New submission. There are no reverse dependencies.

Checked file: sha256 f1947a1e9a0ad30b8ef6075d1e8576c4baccd81b6cb6f18ec46c46b678502777, commit eb3fb1435d0d566db7600a8dd744946d0646331e, built before this comments file.

## Test environments

* Local Ubuntu 24.04, R 4.5.3 (2026-03-11), old-release, lavaan 0.7-3. NOT_CRAN=true and _R_CHECK_CRAN_INCOMING_USE_ASPELL_=true. R CMD check --as-cran of this tarball.
* The same source tree on lavaan 0.6-21, NOT_CRAN=true, testthat::test_local(). Finished with status 0.
* Docker rocker/r-devel. The check command was RD CMD check --as-cran. R Under development (unstable) (2026-11-08 r90859). The PDF manual was built.
* Minimum R. On R 4.1.0, current CRAN does not install this dependency set: Matrix 1.7-6 requires R >= 4.4, and the Matrix bundled with R 4.1 is older than the Matrix version lme4 links to. Depends is R (>= 4.2.0). On rocker/r-ver:4.2, R 4.2.3 (2023-03-15), this tarball installed from current CRAN with lavaan 0.7-3, survey 4.5, mice 3.19.0, mitml 0.4-5, and testthat 3.3.2. testthat::test_local() under NOT_CRAN=true reported FAILED 0, SKIPPED 0, PASSED 518, WARNINGS 0.
* Win-builder. Uploaded once at 2026-11-14 16:46:45 UTC to R-devel, R-oldrelease, and R-release. Mail from before that time was ignored. The R version in each line is copied from that result's 00check.log.
  * R-oldrelease: R version 4.5.3 (2026-03-11 ucrt). https://win-builder.r-project.org/5ecSkSvQar1p/00check.log
  * R-release: R version 4.6.2 (2026-10-29 ucrt). https://win-builder.r-project.org/Aso1hKVGb621/00check.log
  * R-devel: R Under development (unstable) (2026-11-13 r90888 ucrt). https://win-builder.r-project.org/7ZNnQ34NE2C8/00check.log
* macOS builder: results page not ready. https://mac.r-project.org/macbuilder/v1/submit returned HTTP 502 from 2026-11-14 16:38 UTC through 17:08 UTC, and on a later check. One r-release post came back in under two seconds as nginx's page-unavailable document and contained no job link, so it was not an accepted submission. r-devel was not posted after the endpoint stayed down. https://mac.r-project.org/reports/status.html lists the MacBuilder submission service as FAILED.

## R CMD check results

0 errors and 0 warnings on the local check, on all three Win-builder logs, and on the Docker R-devel check.

The local check and all three Win-builder logs are Status: 1 NOTE.

    New submission

    Possibly misspelled words in DESCRIPTION:
      Bentler (22:13)
      Oberski (8:64)
      Satorra (22:5)
      subpopulation (6:67)

Satorra, Bentler, and Oberski are surnames in the citations. subpopulation is the survey term for a domain taken with subset(). The Description spells out structural equation models, maximum likelihood, multiple imputation, confirmatory factor analysis, primary sampling units, and weighted least squares. These four words are kept. The New submission note needs no code change.

Docker R-devel is Status: 3 NOTEs. One is the new submission. The other two come from the image, not from the package:

* Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.
* Skipping checking HTML validation: no command 'tidy' found.

checking PDF version of manual ... OK. checking tests ... OK.

On Win-builder, examples, tests, the PDF manual, and the HTML manual were OK for R 4.5.3, R 4.6.2, and R-devel r90888.

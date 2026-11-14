## Test environments

* Local Ubuntu 24.04.5, R 4.5.3, lavaan 0.7-3. `R CMD check --as-cran` with `_R_CHECK_CRAN_INCOMING_USE_ASPELL_=TRUE`.
* R 4.2.3, lavaan 0.7-3, testthat 3.3.2. This is the oldest R allowed by Depends.
* R 4.6.2, lavaan 0.7-3.
* R-devel r90859 (2026-11-08), lavaan 0.7-3. The PDF manual was built.
* Win-builder R 4.6.2 and R 4.5.3, checked 2026-11-14. Each reported 1 NOTE: New submission.
* Win-builder R-devel was uploaded on 2026-11-14. The result mail had not arrived when these comments were written.
* The macOS builder returned HTTP 502 on 2026-11-14 for both R-release and R-devel, so there is no macOS result.

## R CMD check results

There were no ERRORs or WARNINGs.

On R 4.5.3 there was 1 NOTE:

* New submission.

The description spell check was clean. `.aspell` is included so the checker accepts Oberski, Satorra, Bentler, and subpopulation.

Win-builder R 4.6.2 and R 4.5.3 reported only that same NOTE.

R 4.6.2 and R-devel r90859 reported the same New submission NOTE, plus two notes from the container:

* `NEWS.md` was not checked because pandoc is not installed.
* HTML validation was skipped because HTML Tidy is not installed.

The local R 4.5.3 check has pandoc and HTML Tidy. It checked `NEWS.md` and the HTML manual.

R 4.2.3 reported the New submission NOTE, plus:

* Possibly invalid DOIs `10.1002/9780470316696` (Rubin 1987) and `10.1093/biomet/86.4.948` (Barnard and Rubin 1999). Both publishers return HTTP 403 to the checker. R 4.5.3, R 4.6.2, and R-devel did not report them.
* Unable to verify the current time, so future timestamps were not checked.
* pandoc is not installed, so `NEWS.md` was not checked there.

Examples, tests, and the PDF manual passed on every completed platform. All of these checks used lavaan 0.7-3. The package depends on lavaan >= 0.6-21.

## Test environments

* Local Ubuntu 24.04, R 4.5.3 (old-release), lavaan 0.6-21 and 0.7-3
* R 4.2.3, the oldest version on which the current suggested packages install
* Win-builder R-devel (2026-11-13 r90888), R 4.6.2, and R 4.5.3, checked 2026-11-14
* R-devel r90859 in a Debian container, with the PDF manual built
* The macOS builder was unavailable on 2026-11-14 (the submission service returned HTTP 502), so there is no macOS result yet

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* New submission.

Possibly misspelled words in DESCRIPTION:

* Bentler
* Oberski
* Satorra
* subpopulation

Satorra, Bentler, and Oberski are surnames. subpopulation is the usual survey term.

The Debian R-devel container also notes that pandoc and HTML Tidy are not installed, so NEWS.md and HTML validation were not checked there. The PDF manual was built.

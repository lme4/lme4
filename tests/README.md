Testing protocols
================

At present `lme4` uses two different testing protocols:

* "vanilla" R testing -- every `.R` file gets run during `R CMD check` and compared with its corresponding `.Rout.save` file (if it exists)
* `testthat` testing: tests in `inst/tests` can be run in the framework of the `testthat` package, specifically via `tests/test-all.R`

There are some unique challenges in testing a package like `lme4` that does numerical computations:

* because of floating-point precision issues, detailed numerical results will differ across platforms, compilers, etc.: this can be handled e.g. by setting an appropriate tolerance in `all.equal`, but we do have to decide on an appropriate tolerance.  Also, this means that we should *not* typically print out full results of fits to the output file, because these will lead to `Rout`/`Rout.save` comparisons being flagged.
* detailed numerical results will also change if we change the details (order etc.) of the internal algorithms (including changes to the underlying computational algorithms in `RcppEigen`); the same tolerance issues arise
* in addition to well-targeted unit tests that test specific aspects of the interface, we really need a battery of tests that may be time-consuming, which will be a problem both for rapid development and for CRAN.  *Solution*:
 * users can set an environment variable `LME4_TEST_LEVEL` which controls which tests will be run. Each `.R` file can contain use
```
if(lme4:::testLevel() >= *) { ... }
```
where that gets the value of an environment variable `LME4_TEST_LEVEL` if
defined as a number and '1' otherwise, and so specify the default testing level.  The general idea would be that 1=default=quick (for CRAN compliance, tests running in a few seconds or less); 2=fairly quick (tests running in <10 seconds); 3=moderate (tests running in <30 seconds); 4=long/exhaustive.  `LME4_TEST_LEVEL` could be set less than 1, or zero, or negative, if we wanted really quick tests.
 * at present this is implemented at the level of each `.R` file, and indeed even at a finer level, by putting `if (testLevel-condition) {}` blocks in the file.  It might be nice to have something more global, but this granularity is convenient too.
 * the `Rout.save` results will differ according to test level: we can either make sure that results are saved at the CRAN default level, or use `.Rbuildignore` to make sure the `Rout.save` files stay off CRAN entirely.



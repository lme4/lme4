## `flexSigmaMinimum`

* `grep -e "FIXME flexSigmaMinimum"`
* check that existing examples work and existing tests pass
  - `git diff master tests` looks OK now: a few tests fail but those
    particular failures are expected; these tests are skipped or adapted
  - `git diff master man` shows that `profile` is broken when `theta`
    is not composed of segments of length `nc*(nc+1)/2`
    * hence TODO: adapt usage of (or generalize) the `*_to_*` functions
      in `vcconv.R`; see, e.g., `devfun2` in `R/profile.R`
  - otherwise looking good ... !
* check that reverse dependencies pass *their* checks
* new tests
  - unit tests for stuff in `R/covariance.R`
  - integration tests for `lmer`, `glmer`, `nlmer`, and
    methods for class `"merMod"`
* new documentation
  - update `vignette("lmer")` (or is that static ... ?)
  - write a `vignette("covariance")` (or whatever)


## not `flexSigmaMinimum`

- "toast-scraping": improve and/or get rid of post-hoc convergence testing
- build reliable downstream-package-testing infrastructure
- mixed models task view, MM package comparison grid
- improve lmList
- stress-test with tibbles
- more examples of using modular structure, e.g. for correlation models
- bring `flexLambda` branch up to the present?
- `glmer`: return `devfun` without re-evaluating
- implement various kinds of marginal predictions
- parallelized likelihood calcs for simple (single/nested groups) models?
- implement Kenward-Roger for GLMMs? (after Stroup)
- AGHQ
    - implement for random-slopes models (at least)
	- write up AGHQ simulation paper
	- diagnostics?
	- importance sampling?
- revisit *post hoc* MCMC sampling
- finish `glmer` paper
- think about sampling weights
- data-sharing models
- back-port improved formula processing from `glmmTMB`
- multinomial and ordinal models; multitype models???


-----
## Old "to do" list (pre-2009)

*items here have uncertain status*

- Change the `terms` object to be the terms for the fixed-effects only so
that the `drop1` method doesn't try to drop the grouping factors for the
random effects.
- Modify the effect of the verbose setting in the Nelder-Mead optimizer. In particular, it should count evaluations but define an "iteration" as a change in the best value encountered so far.
- The paper by Sophia Rabe-Hesketh et al describes a spherical form of the Gauss-Hermite quadrature formula.  Look that up and use it.
- Because the Gauss-Hermite quadrature is formed as a sum, it is necessary to divide the contributions to the deviance according to the levels of the random effects.  This means that it is only practical to use AGQ when the response vector can be split into sections that are conditionally independent. As far as I can see this will mean a single grouping factor only.
- Allow for a matrix of responses in lmer so multiple fits can be performed without needing to regenerate the model matrices.
- Determine what a `coef` function should do for multiple, possibly non-nested, grouping  factors.
- add nicer (more realistic?) pedigree examples and tests
- document print(<mer>) including an example  print(<lmer>, corr = FALSE) and one with many fixed effects (*) and print(<lmer>, symbolic.cor = TRUE)


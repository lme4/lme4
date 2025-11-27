## `flexSigmaMinimum`

* document that theta != par != profPar (getVC)
   * theta = unique scaled cholesky components
   * par = optimization parameters
   * profPar = profiling parameters
   * getVC() gets scaled SDs and correlation parameters
* document clearly that scale = "varcor" is only allowed for us() matrices
* DOCUMENT that internal order of profiling parameters has changed; was previously in Lambda-matrix order
      (i.e. sd1, cor12, sd2 for a 2x2), now in [more sensible, glmmTMB-matching] (sd1, sd2, cor12) [this is going to be a "user-visible change"]
* check that reverse dependencies pass *their* checks (!!)
* new tests
  - unit tests for stuff in `R/covariance.R`
  - integration tests for `lmer`, `glmer`, `nlmer`, and
    methods for class `"merMod"`
* new documentation
  - write a `vignette("covariance")` (or whatever)
* adapt `formatVC` to handle new covariance structures smoothly
* make `VarCorr()` structure-aware: in particular, be careful *not* to generate full cov/cor matrix, unless explicitly requested, for structured cov matrices (to avoid huge and redundant structures)


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


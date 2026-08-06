# Log-link Gaussian mini parameter-recovery survey (GH #936)

Much more limited sibling of `../paramsurvey/`: investigates
[GH #936](https://github.com/lme4/lme4/issues/936) ("matching up
gaussian(link="log") results across packages") specifically, rather than
the broader Gamma GLMM bias question `../paramsurvey/` covers. One
dataset (`nlme::Rail`, the exact example from the issue), one family
(`gaussian(link="log")`), a single scalar random intercept, and four
methods (not six -- no R-level joint-phi/PIRLS-phi devfun arms this time):

1. **glmmTMB** -- reference/gold standard.
2. **glmer** -- current dev build (this branch's fixes).
3. **mgcv::gam** (`s(Rail, bs="re")`, `method="ML"`).
4. **MASS::glmmPQL** -- penalized quasi-likelihood, not full ML; its
   `negll` is always `NA` (no comparable marginal likelihood).

Per the (partial) resolution already recorded in
`../README_Gamma_GLMMs.md` ("Backward-compatibility switch, and nAGQ>1
scoping"/GH #936 section): this branch's dispersion-profiling fix
improved `glmer`'s logLik on the real Rail data dramatically (-72.3 to
-65.1, much closer to glmmTMB's -64.6) but left `sigma`/`theta` still
visibly off from glmmTMB/mgcv/glmmPQL, on this one (tiny, n=18) real
dataset. This mini-survey exists to check whether that residual sigma
gap is systematic (a real remaining bug) or just small-sample noise, by
simulating repeatedly instead of relying on one dataset.

## Scripts (run in this order)

| script | purpose |
|---|---|
| `toolkit.R` | shared library: per-method fit wrappers, result-row helpers. Sourced by everything below. |
| `01_prep_rail.R` | fit a glmmTMB reference to the real `Rail` data, round to "pretty" true parameters (beta=4.1, sd=0.4, sigma=4.0), simulate `B` new datasets from those parameters (real Rail design, 6 rails x 3 obs) via `glmmTMB::simulate_new()`. Takes `B` as `commandArgs()[[1]]` (default 10). |
| `02_fit_glmmTMB.R` | fit glmmTMB to each simulated dataset. |
| `03_fit_glmer.R` | fit current lme4 dev build (`glmer()`) to each simulated dataset. |
| `04_fit_mgcv.R` | fit `mgcv::gam()` to each simulated dataset. |
| `05_fit_pql.R` | fit `MASS::glmmPQL()` to each simulated dataset. |
| `06_analysis.R` | combine the four per-method result files, print status/singular counts, timing, and parameter-recovery summaries; save `rail_results_combined.rds`. |

## Outputs

Per-method fits are saved as `rail_results_<method>.rds`;
`06_analysis.R` combines them into `rail_results_combined.rds`. `.rds`
files are gitignored (see `../paramsurvey/`'s convention) -- only the
scripts and this README are meant to be committed.

## Status

Set up but not yet run at scale (B=10 to start, per this session's
request) -- see the parent conversation for results once the sweep has
been run.

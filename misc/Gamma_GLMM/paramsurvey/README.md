# Parameter-recovery survey toolkit

Broader successor to the multistart diagnostic (see `../README_Gamma_GLMMs.md`
§11-12 for full background/motivation). For each of four real datasets, fit
**six** methods to **B** freshly-simulated datasets drawn from "pretty"
(rounded) parameters near a real reference fit, and compare parameter
recovery, reliability (clean/warning/error/singular), timing, and
-2\*logLik across methods.

## Six methods

1. **glmmTMB** — reference/gold standard, full joint TMB optimization.
2. **joint-phi (R)** — phi as a first-class outer `bobyqa` parameter
   alongside theta/beta (R-level devfun).
3. **PIRLS/digamma (R)** — nested fixed-point: PIRLS for beta/theta given
   phi, phi re-derived each outer iteration as the true Gamma conditional
   MLE (digamma equation).
4. **PIRLS/moment (R)** — same nested fixed-point, but phi = deviance/n
   (crude moment plug-in) instead of the conditional MLE.
5. **PIRLS/moment (C++)** — the compiled lme4 fix (`glmerLaplace()`,
   current dev version), same moment-phi algorithm as #4, isolates
   R-vs-compiled overhead.
6. **PIRLS/fixed-phi (CRAN)** — unmodified lme4 2.0-6: the pre-fix,
   buggy method that never actually profiles phi during fitting (the bug
   this whole investigation exists to fix). Included as the "what users
   get today from CRAN" baseline, not a fair phi-estimator alternative.

## Four datasets

- **epil2 (simple)**: `y ~ trt + (1|subject)`, single random intercept.
- **epil2 (complex)**: `y ~ Base*trt + Age + Visit + (Visit|subject)`,
  correlated random intercept+slope, uncentered covariates (Base, Age).
- **Report4BB**: `crate ~ (1|location) + (1|fyear)`, two independent
  small-group random-intercept terms, no correlation parameter.
- **schizophrenia**: `imps79 ~ TxDrug*Week + (1|id)`, single random
  intercept, much larger dataset (1603 obs, 437 subjects).

## Scripts (run in this order)

| script | purpose |
|---|---|
| `toolkit.R` | shared library: joint-phi/PIRLS-phi devfuns, RE sd/corr extraction, per-method fit wrappers, result-row helpers. Sourced by everything below. |
| `01_prep_epil2_simple.R`, `01_prep_epil2_complex.R`, `01_prep_report4bb.R`, `01_prep_schizophrenia.R` | one per dataset: fit a glmmTMB reference (with a regularizing ranef prior where needed) to real data, round to "pretty" true parameters, simulate B new datasets from those parameters using the real design. Takes `B` as `commandArgs()[[1]]` (default 10). |
| `02_fit_glmmTMB.R` | fit glmmTMB (no prior) to each simulated dataset. |
| `03_fit_jointphi.R` | fit the joint-phi devfun to each simulated dataset. |
| `04_fit_lme4current.R` | fit current lme4 dev build (`glmer()`) to each simulated dataset. |
| `05_fit_lme4old.R` | fit unmodified lme4 2.0-6 (isolated library) to each simulated dataset. Must load the old lme4 before `toolkit.R`'s `library(lme4)` call, since R won't reload an already-attached package from a different `lib.loc`. |
| `06_fit_pirls_phi.R` | fit the R-level nested-fixed-point PIRLS devfun for a given `phiType` (`digamma` or `moment`) to each simulated dataset. |
| `07_analysis.R` | combine one dataset's six per-method result files, print status/singular counts, timing, and parameter-recovery summaries; save `<example>_results_combined.rds`. |
| `08_summary_plots.R` | build all four summary plots (see below) from all four datasets' per-method result files. |

Fitting scripts 02-06 all take `example` (dataset name) as their first
`commandArgs()` argument and `MC_CORES` as their last (parallelized via
`mclapply`); `06_fit_pirls_phi.R` additionally takes `phiType` as its
second argument.

## Outputs

Per-dataset, per-method fits are saved as
`<example>_results_<method>.rds`; `07_analysis.R` combines them into
`<example>_results_combined.rds`. `08_summary_plots.R` produces:

- **`param_summary_distrib.png`** — per-parameter estimates across all
  six methods and four datasets: full per-replicate distribution
  (violin + boxplot + points), one facet row per parameter, one panel
  per dataset (patchwork).
- **`param_summary_stderr.png`** — same parameter summary, but mean ±
  2 SE (`geom_pointrange`) instead of the full distribution.
- **`param_summary_stderr_newonly.png`** — same as
  `param_summary_stderr.png`, excluding PIRLS/fixed-phi (CRAN)/lme4
  2.0-6, so the five current methods' much smaller differences aren't
  swamped by that method's much larger bias.
- **`time-negll_summary.png`** — two panels stacked vertically via
  patchwork (not `facet_grid`, since the value ranges differ too much
  across datasets for one shared axis per row): elapsed time (all six
  methods) on top, and paired per-replicate Δ(-2\*logLik) vs glmmTMB
  (same simulated dataset, method minus glmmTMB) on the bottom, which
  cancels between-replicate variation and isolates each method's
  systematic gap from glmmTMB. The bottom panel excludes glmmTMB itself
  and PIRLS/fixed-phi (CRAN)/lme4 2.0-6 (see script comments for why).

## Caveats

- All scripts hardcode this session's scratch-directory path
  (`wd <- "/tmp/claude-.../scratchpad/param_survey"`) — update before
  rerunning outside this session.
- `05_fit_lme4old.R` expects an isolated lme4 2.0-6 library at
  `../lme4_206_lib` (sibling of this directory, in the scratch tree).

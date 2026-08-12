# Parameter-recovery survey toolkit

Broader successor to the multistart diagnostic (see `../README_Gamma_GLMMs.md`
§11-12 for full background/motivation). For each of four real datasets, fit
**up to nine** methods to **B** freshly-simulated datasets drawn from "pretty"
(rounded) parameters near a real reference fit, and compare parameter
recovery, reliability (clean/warning/error/singular), timing, and
-2\*logLik across methods.

## Nine methods

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
7. **MixedModels.jl (pa/dispersion-again)** — Julia's
   [JuliaStats/MixedModels.jl](https://github.com/JuliaStats/MixedModels.jl)
   branch `pa/dispersion-again`, which implements essentially the same
   fix as #5: profiles phi-hat = pwrss/n at the converged PIRLS mode,
   matching lme4's `sigma()`. Run via a separate, isolated Julia project
   in `julia/` (see below) -- not part of the R `toolkit.R` machinery.
8. **SAS GLIMMIX (RSPL/PQL)** — `PROC GLIMMIX`'s default estimation
   method, residual/restricted pseudo-likelihood via linearization
   (SAS's implementation of PQL, the same family as `MASS::glmmPQL`).
   A genuinely different bias mechanism from #1-7 (Breslow & Clayton
   1993's classic PQL bias, not the phi-profiling issue this whole
   investigation is about) -- included as another point of comparison,
   not a redundant one. Doesn't produce a real marginal -2\*logLik (its
   fit statistics are on the linearized pseudo-data scale), so `negll`
   is left missing for this method.
9. **SAS GLIMMIX (Laplace)** — same `PROC GLIMMIX`, `METHOD=LAPLACE`: a
   proper marginal-likelihood Laplace approximation, the most directly
   comparable SAS method to lme4's/Julia's default (`nAGQ=1`-equivalent).
   Run via a separate `sas/` directory (see below) -- not part of the R
   `toolkit.R` machinery. **`METHOD=QUAD` (adaptive Gauss-Hermite
   quadrature) deliberately not included**: like lme4's `nAGQ>1`, it's
   restricted to a single scalar random effect per subject, and none of
   methods #1-8 use anything beyond `nAGQ=1`/Laplace either, so it
   wouldn't be a like-for-like addition without first extending the
   comparison to quadrature everywhere.

## Four datasets

- **epil2 (simple)**: `y ~ trt + (1|subject)`, single random intercept.
- **epil2 (complex)**: `y ~ Base*trt + Age + Visit + (Visit|subject)`,
  correlated random intercept+slope, uncentered covariates (Base, Age).
- **Report4BB**: `crate ~ (1|location) + (1|fyear)`, two independent
  small-group random-intercept terms, no correlation parameter. Data
  file (`report4bb_data/data_4_article_clickrates_deep_dive.rda`, ~50KB)
  copied directly from
  [TiagoAMarques/Report4BB](https://github.com/TiagoAMarques/Report4BB)
  into this directory -- self-contained, no repo clone needed.
- **schizophrenia**: `imps79 ~ TxDrug*Week + (1|id)`, single random
  intercept, much larger dataset (1603 obs, 437 subjects).

## Scripts (run in this order)

| script | purpose |
|---|---|
| `toolkit.R` | shared R library: joint-phi/PIRLS-phi devfuns, RE sd/corr extraction, per-method fit wrappers, result-row helpers. Sourced by scripts 02-06. |
| `01_prep_epil2_simple.R`, `01_prep_epil2_complex.R`, `01_prep_report4bb.R`, `01_prep_schizophrenia.R` | one per dataset: fit a glmmTMB reference (with a regularizing ranef prior where needed) to real data, round to "pretty" true parameters, simulate B new datasets from those parameters using the real design. Takes `B` as `commandArgs()[[1]]` (default 10). |
| `02_fit_glmmTMB.R` | fit glmmTMB (no prior) to each simulated dataset. |
| `03_fit_jointphi.R` | fit the joint-phi devfun to each simulated dataset. |
| `04_fit_lme4current.R` | fit current lme4 dev build (`glmer()`) to each simulated dataset. |
| `05_fit_lme4old.R` | fit unmodified lme4 2.0-6 (isolated library) to each simulated dataset. Must load the old lme4 before `toolkit.R`'s `library(lme4)` call, since R won't reload an already-attached package from a different `lib.loc`. |
| `06_fit_pirls_phi.R` | fit the R-level nested-fixed-point PIRLS devfun for a given `phiType` (`digamma` or `moment`) to each simulated dataset. |
| `07_analysis.R` | combine one dataset's per-method result files (whichever of the seven exist), print status/singular counts, timing, and parameter-recovery summaries; save `<example>_results_combined.rds`. |
| `08_summary_plots.R` | original 6-method (methods #1-6)/5-dataset (including `epil2_phigt1`) summary plots -- thin wrapper around `plotlib.R::make_summary_plots()`, unchanged output filenames. |
| `09_export_csv_for_julia.R` | export one example's B simulated replicate datasets to `julia/data/<example>/rep_NNN.csv` (only the formula's columns), so the Julia side can fit them without RCall/RDS. |
| `10_ingest_julia_results.R` | read `julia/results_<example>.csv` back into R, reshape into the same `resultsToDF()` column shape as the other methods, save as `<example>_results_juliaMixedModels.rds`. |
| `11_summary_plots_julia.R` | 4-method (glmmTMB/joint-phi/lme4current/Julia) comparison plots via `plotlib.R::make_summary_plots()`, output prefixed `julia_compare_` so it never collides with `08_summary_plots.R`'s files. |
| `12_export_csv_for_sas.R` | export one example's B simulated replicate datasets to a single stacked `sas/data/<example>.csv` (all B replicates, keyed by a `rep` column) -- SAS's `PROC GLIMMIX BY rep;` fits all replicates in one step, unlike Julia's per-file loop. |
| `13_ingest_sas_results.R` | read `sas/results_<example>_<rspl|laplace>.csv` back into R, reshape into the same `resultsToDF()` column shape as the other methods, save as `<example>_results_sasRSPL.rds` / `<example>_results_sasLaplace.rds`. Takes `example` and `rspl`/`laplace` as its two `commandArgs()`. |
| `plotlib.R` | shared, parameterized plotting library (factored out of the original `08_summary_plots.R`): `make_summary_plots(examples, methods, out_prefix, ...)` builds all four plot types for any subset of the registered methods/datasets; `make_stderr_plot(examples, methods, outfile, ...)` builds just the mean±2SE pointrange plot for one-off subset comparisons that don't need the full three-plot set (e.g. `08_summary_plots.R`'s `_newonly` plot). Master method/colour/dataset-label registries live at the top of this file -- extend them there when adding a new method or dataset, never reassign an existing method's colour. |
| `run_remaining_datasets.sh`, `run_full_survey.sh <B>` | bash orchestration: per dataset, run prep -> CSV export -> R fits (glmmTMB/joint-phi/lme4current concurrently via mclapply) -> Julia fit (via the warm DaemonMode server) -> ingest -> `07_analysis.R`. `run_full_survey.sh` is the general one (takes B as an argument, covers all four datasets incl. epil2_simple); `run_remaining_datasets.sh` was the one-off B=50 calibration version, kept for reference. |

Fitting scripts 02-06 all take `example` (dataset name) as their first
`commandArgs()` argument and `MC_CORES` as their last (parallelized via
`mclapply`); `06_fit_pirls_phi.R` additionally takes `phiType` as its
second argument. `09_export_csv_for_julia.R` and
`10_ingest_julia_results.R` take only `example`.

## Julia side (`julia/`)

Isolated Julia project (own `Project.toml`/`Manifest.toml`, *not* the
global Julia environment) with MixedModels.jl installed from
`JuliaStats/MixedModels.jl@pa/dispersion-again` via
`Pkg.add(PackageSpec(url=..., rev=...))`, plus CSV/DataFrames/GLM/
StatsModels/DaemonMode.

| file | purpose |
|---|---|
| `setup.jl` | one-time environment setup: installs the branch + deps, precompiles. Rerun to pick up new commits on the branch. |
| `fitlib.jl` | shared fitting/extraction routine, `run_survey(example, form, family, link, sdcorr_spec, nbeta)` -- handles both single-group and two-independent-group RE structures; writes a CSV matching the R side's `resultsToDF()` columns. |
| `fit_epil2_simple.jl`, `fit_epil2_complex.jl`, `fit_report4bb.jl`, `fit_schizophrenia.jl` | thin per-dataset scripts, each just declares its `@formula` and calls `run_survey()`. |
| `daemon_start.sh [port]` | starts a persistent DaemonMode.jl server (default port 3141) with this project activated, so MixedModels.jl loads/precompiles once and stays warm across script runs -- avoids paying Julia's JIT cost on every invocation. |
| `run_via_daemon.sh <script.jl> [args]` | runs a script against the warm daemon instead of a fresh `julia` process. ~9x faster than cold start once warm (measured: 27.5s cold vs 2.97s warm for a 50-replicate fit). |

Extending to a new dataset: add a `fit_<example>.jl` following the
existing pattern (formula + `run_survey()` call with the right
`sdcorr_spec` -- a single grouping-factor `Symbol` for one RE term (with
1 or 2 correlated components), or a `Tuple{Symbol,Symbol}` for two
independent single-term grouping factors).

## SAS side (`sas/`)

**Written but not yet run against real SAS** -- there's no SAS install
in this environment, so this arm was authored from `PROC GLIMMIX`/`ODS
OUTPUT` documentation rather than verified interactively. Expect the
first real run to need debugging; `sas/fitlib.sas`'s header and inline
comments flag the specific things most likely to need adjustment (exact
`CovParm` label text, `ConvergenceStatus` column names, Fit Statistics
row labels -- these can vary a little by SAS release/options).

| file | purpose |
|---|---|
| `fitlib.sas` | shared macro, `%fit_glimmix(example=, method=, classvars=, modelrhs=, nbeta=, randomstmt=, ...)`: fits `PROC GLIMMIX` to all B replicates of one example via `BY rep;` processing under one `METHOD=` (`rspl` or `laplace`), reshapes `ParameterEstimates`/`CovParms`/`ConvergenceStatus` into the same column shape the other methods use, writes `results_<example>_<method>.csv`. Its `indir=`/`outdir=` defaults (and the plain `%include "fitlib.sas"` in each `fit_<example>.sas`) assume SAS is launched **from inside `sas/`** as the working directory -- `run_sas_survey.sh` (below) does this via `cd sas && ...`. |
| `fit_epil2_simple.sas`, `fit_epil2_complex.sas`, `fit_report4bb.sas`, `fit_schizophrenia.sas` | thin per-dataset scripts, each just supplies the model/class/random-effect specifics and calls `%fit_glimmix` twice (rspl + laplace). |

`../run_sas_survey.sh` orchestrates all four datasets end to end (export
CSV -> `PROC GLIMMIX` RSPL+Laplace -> ingest -> combined analysis),
mirroring `run_full_survey.sh`'s structure. Point it at your actual SAS
batch executable first (install-dependent, not guessable from here):
`SAS_EXE=/path/to/sas ./run_sas_survey.sh`. SAS's own process exit code
isn't a reliable success signal, so the script greps each `.log` for
`^ERROR` instead of trusting `$?`, and keeps going through the remaining
datasets even if one fails.

Two things worth knowing before running these:
- **Fixed-effect term order matters and isn't auto-corrected.** R's
  `terms()` moves interaction terms to the end of the model regardless
  of where they're written in the formula (e.g. `Base*trt + Age + Visit`
  becomes `Base, trt, Age, Visit, Base:trt` in `beta_pretty`'s order);
  SAS's `MODEL` statement does not do this. The `modelrhs=` arguments in
  each `fit_<example>.sas` are written in R's final term order by hand
  (e.g. `Base trt Age Visit Base*trt`, not the more natural-looking
  `Base*trt Age Visit`) so the positional beta-matching in
  `13_ingest_sas_results.R` lines up -- same positional-matching
  approach the Julia arm uses, since SAS's own effect-name strings don't
  string-match R's factor-contrast names either.
- **`negll` is only filled in for `method=laplace`.** RSPL's fit
  statistics are on the linearized pseudo-data scale, not a real
  marginal -2\*logLik, so reporting a number there would look
  comparable to the other methods' `negll` without actually being so.
- **`time_sec` is a whole-batch average, not a per-replicate
  measurement.** `BY rep;` processing fits all B replicates inside one
  `PROC GLIMMIX` step, so there's no per-replicate wall-clock split the
  way R's `system.time()`-per-fit or Julia's per-file timing gives --
  every row for a given (dataset, method) gets the same value
  (total elapsed / B).

## Outputs

Per-dataset, per-method fits are saved as
`<example>_results_<method>.rds`; `07_analysis.R` combines them into
`<example>_results_combined.rds`. `08_summary_plots.R`/
`11_summary_plots_julia.R` (via `plotlib.R`) produce, per call:

- **`<prefix>param_summary_distrib.png`** — per-parameter estimates
  across the chosen methods/datasets: full per-replicate distribution
  (violin + boxplot + points), one facet row per parameter, one panel
  per dataset (patchwork).
- **`<prefix>param_summary_stderr.png`** — same parameter summary, but
  mean ± 2 SE (`geom_pointrange`) instead of the full distribution.
- **`<prefix>time-negll_summary.png`** — two panels stacked vertically
  via patchwork (not `facet_grid`, since the value ranges differ too
  much across datasets for one shared axis per row): elapsed time (all
  chosen methods) on top, and paired per-replicate Δ(-2\*logLik) vs a
  reference method (default glmmTMB; same simulated dataset, method
  minus reference) on the bottom, which cancels between-replicate
  variation and isolates each method's systematic gap from the
  reference.
- `08_summary_plots.R` additionally produces
  **`param_summary_stderr_newonly.png`** (excludes PIRLS/fixed-phi
  (CRAN)/lme4 2.0-6) via an explicit third `make_plot_pointrange()`
  call in that script, not part of `make_summary_plots()` itself.

`<prefix>` is `""` for `08_summary_plots.R` (original filenames,
unchanged) and `"julia_compare_"` for `11_summary_plots_julia.R` --
different prefixes never collide, so re-running either never overwrites
the other's output. Any new comparison (different method/dataset
subset) should pick its own `out_prefix` when calling
`make_summary_plots()` directly, for the same reason.

## Real-data findings so far (Julia vs lme4/R methods, B=50 and B=500)

- **epil2 (simple)**: Julia's random-intercept SD estimate runs
  noticeably low vs every R method (e.g. B=500 median 0.640 vs
  glmmTMB/joint-phi/lme4current's ~0.77, true 0.80) -- consistent at
  both B=50 and B=500, so it's a real, reproducible gap, not noise.
- **Report4BB** (two small independent grouping factors, true SDs tiny
  relative to n): Julia has a *much lower* singular-fit rate than every
  R method (e.g. B=500: 66/500 vs 276-323/500 for the R methods) --
  opposite direction from the epil2-simple finding.
- **schizophrenia** (larger n, single RE): only a small Julia-vs-lme4
  gap; Julia's phi estimate was actually the more accurate of the two.
- **epil2 (complex)** (correlated intercept+slope RE) -- the most
  substantial finding: **both** lme4's C++ moment-phi fix and Julia's
  port badly overestimate the random-**slope** variance (B=50: sd2
  medians 0.82 (lme4) / 0.71 (Julia) vs glmmTMB 0.27 / joint-phi 0.22,
  true 0.15), while reporting *zero* singular fits (unlike glmmTMB/
  joint-phi's ~35/50) -- i.e. a real, confident, non-boundary bias
  shared by both moment-phi implementations, not something unique to
  either codebase. Not yet root-caused; see the project memory file for
  the current working hypothesis and suggested follow-up.

## Caveats

- `05_fit_lme4old.R` expects an isolated lme4 2.0-6 library at
  `lme4_206_lib` (sibling of this directory) -- not committed, rebuild
  before rerunning (see comment at the top of that script).
- The `sas/` scripts have never been run against real SAS -- see the
  "SAS side" section above before trusting their output.

# dispsim: dispersion estimators for Gamma GLMMs

Small simulation study comparing estimators of the dispersion phi (and the
other parameters) for free-dispersion GLMMs, prompted by the deviance
offset in `misc/TODO_AGQ.md` (the deviance uses two different phi values).

## Files

- `dispsim_funs.R`: simulation (`simdata()`), fitting (`fit_all()`), and
  post-hoc reporting-only phi estimators (`phi_alt()`).
- `run_dispsim.R`: `Rscript run_dispsim.R nsim ngrp nper shape sd_re workers`
  (defaults 100, 20, 20, 0.5, 0.5, 12). Runs replicates in parallel with
  furrr (`multicore`), one BLAS/OpenMP thread per worker. Replicate `i`
  uses `set.seed(1000 + i)`, so results don't depend on scheduling. Writes
  `dispsim_<setup>.rds`: one row per replicate x method, with the true
  values and setup (including lme4/glmmTMB versions and wall time) as
  attributes.
- `summarize_dispsim.R`: `Rscript summarize_dispsim.R dispsim_<setup>.rds`:
  failure counts, timing, and bias / SE(bias) / RMSE per method and parameter.
- `dispsim_<setup>_summary.txt`: saved output of `summarize_dispsim.R`.

## Design

Gamma(log), `y ~ x + (1|g)`, x ~ N(0,1), beta = (1, 0.5). Methods:

- `glmer_n`, `glmer_n-p`, `glmer_n-q`: `glmerControl(disp_dof_correction =
  "n" / "n-p" / "n-q")`, i.e. phi = dev/n, dev/(n - p) or dev/(n - q) inside
  PIRLS, q = rank of [X Z]. The corrected fits are warm-started from `glmer_n`
  (this saves little time: about 4.1 vs 4.3 s per fit).
- `glmmTMB_ML`, `glmmTMB_REML`.
- Reporting-only phi alternatives, computed after each glmer fit at its
  conditional means (theta and beta unchanged): `digamma` (conditional ML
  given mu, `log(nu) - digamma(nu) = D/(2n)`, phi = 1/nu), `pearson_p`
  (X^2/(n - p)) and `pearson_q` (X^2/(n - q)).

## Results

### g20_n20_shape0.5_sd0.5 (20 clusters x 20 obs, phi = 2, RE sd = 0.5), 2026-10-26

100 replicates, 130 s wall on 12 workers (about 15 s CPU per replicate;
glmer about 4 s per fit, glmmTMB under 1 s). No errors, warnings, or
singular fits. Full output in `dispsim_g20_n20_shape0.5_sd0.5_summary.txt`.

| phi estimator | bias (SE about 0.01-0.02) | RMSE |
|---|---|---|
| glmer dev/n | +0.48 | 0.51 |
| glmer dev/(n-p) | +0.49 | 0.52 |
| glmer dev/(n-q) | +0.62 | 0.64 |
| glmmTMB ML | +0.004 | 0.11 |
| glmmTMB REML | +0.008 | 0.11 |
| glmer + digamma (reporting only) | -0.04 | 0.12 |
| glmer + Pearson X^2/(n-p) | -0.19 | 0.27 |
| glmer + Pearson X^2/(n-q) | -0.09 | 0.22 |

- At large dispersion the deviance-based moment estimator is about 24% high,
  in the direction predicted by D/n ~ phi(1 + phi/6) (a small-phi expansion,
  which overpredicts at phi = 2: 2.67 vs 2.48 observed); the n-p and n-q corrections make it
  worse. The choice of denominator matters much less than the choice of
  estimator.
- Computing the digamma (conditional ML) phi after the fit, from the same
  conditional means, nearly matches glmmTMB ML (-2% bias, RMSE 0.116 vs
  0.111). The result is the same whichever denominator was used inside PIRLS.
  Pearson estimators are biased low and noisier.
- RE sd: all glmer variants have bias +0.008 to +0.017 and RMSE about 0.096;
  glmmTMB ML -0.034 (RMSE 0.106), REML -0.015 (0.104).
- beta0: about -0.04 bias (about 2.8 SE) for all glmer variants and glmmTMB
  REML, but not glmmTMB ML. Not investigated.
- This differs from the earlier finding (README_Gamma_GLMMs.md, paramsurvey)
  that moment estimators were 15-25% *low* relative to glmmTMB, which
  motivated the n-q correction (partly on gaussian(log), where dev/n is the
  exact conditional ML). Reconciled by the factorial below.
- Earlier history of the digamma estimator: README_Gamma_GLMMs.md §6, §7
  (in lme4pureR, digamma *inside* PIRLS gave -16% RE-variance bias at
  shape 0.5, B=30, vs +4.6% for moment) and §9 (compiled code uses moment
  only). Using digamma only for reporting avoids that issue.

### Factorial: phi x cluster size, 2026-10-26

`run_factorial.sh`: phi in {0.2, 1, 2} (shape 5, 1, 0.5) x cluster size in
{4, 20}, 20 clusters, RE sd 0.5, 100 replicates per cell (about 9 minutes
in total on 12 workers). q/n = 21/80 = 0.26 for cluster size 4 (close to epil2's
0.28) and 21/400 = 0.05 for 20. Adds `phi_digamma_q` (digamma with D/(2(n - q))).
Combined table: `Rscript combine_dispsim.R` -> `dispsim_combined_summary.csv`
(long format, relative bias / SE / relative RMSE). Per-cell output in
`dispsim_g20_n*_summary.txt`. Singular fits (included in all summaries):
up to 11/100 for glmer and 22/100 for glmmTMB ML in the cluster size 4, phi = 2
cell, 0 for cluster size 20.

phi, relative bias (SE about 0.01-0.03):

| | 0.2 / 4 | 1 / 4 | 2 / 4 | 0.2 / 20 | 1 / 20 | 2 / 20 |
|---|---|---|---|---|---|---|
| glmer dev/n | -0.17 | +0.06 | +0.18 | -0.02 | +0.11 | +0.24 |
| glmer dev/(n-q) | +0.15 | +0.41 | +0.57 | +0.04 | +0.17 | +0.31 |
| digamma (cond. ML) | -0.19 | -0.07 | -0.06 | -0.05 | -0.03 | -0.02 |
| digamma, n-q | +0.09 | +0.21 | +0.22 | +0.01 | +0.01 | +0.03 |
| Pearson, n-q | +0.05 | +0.12 | +0.09 | -0.01 | -0.03 | -0.05 |
| glmmTMB ML | +0.00 | +0.00 | -0.00 | -0.00 | -0.00 | +0.00 |

(Reporting-only rows shown for the `glmer_n` fit; they hardly depend on
which glmer arm they come from.)

**Explains the paramsurvey discrepancy.** dev/n has two competing biases:
an upward one from the Gamma deviance's nonlinearity (about phi/6 for small
phi, growing with phi) and a downward one from fitting the conditional
modes (growing with q/n). The paramsurvey's "15-25% low" (R-level
PIRLS/moment, e.g. epil2_simple -22%) is the small-phi, large-q/n corner
(here -17%). The +24% in the first run above is the large-phi, small-q/n
corner. epil2_phigt1 (phi = 1.2, q/n = 0.28, PIRLS/moment +0.7%) is the
near-cancelling case (here phi = 1, cluster size 4: +6%).

Other observations:
- **n - q over-corrects everywhere**, even at small phi (+15% at 0.2 / 4;
  cf. lme4current +9% on epil2_simple). Shrunken conditional modes use fewer
  effective df than rank [X Z]; an effective-df (hat-matrix trace)
  correction might be closer. Not tried.
- Digamma (conditional ML) removes the nonlinearity bias but keeps the
  overfitting bias (-19% at 0.2 / 4); its n - q version over-corrects
  (+9 to +22%). Only glmmTMB ML/REML (integrating over the random effects)
  is unbiased in every cell, and it has the lowest phi RMSE everywhere.
- **RE sd** (relative bias / RMSE at 2 / 4): glmer dev/n -0.02 / 0.54,
  n - q +0.31 / 0.62, glmmTMB ML -0.26 / 0.65, REML -0.18 / 0.62. In
  small clusters glmer dev/n has the best RMSE for the RE sd despite its
  biased phi; n - q inflates the RE sd. With cluster size 20 all biases are within
  about 7%.
- **beta0** is biased low, increasingly with phi and small clusters (2 / 4:
  glmer -0.14 to -0.19, glmmTMB REML -0.11, ML -0.04). Not investigated.

## Next step: effective-df (hat-matrix trace) correction (not yet tried)

phi = dev/(n - edf), with edf = tr(H) for the penalized weighted
least-squares problem at the converged mode:

edf = ||C_L||_F^2 + ||C_R||_F^2,
C_L = L^{-1} P Lambda' Z' W^{1/2},  C_R = R_X^{-T} (X' W^{1/2} - R_ZX' C_L)

(the same computation as `hatvalues.merMod()`, `R/lmer.R`). edf runs from p
(theta -> 0) to rank [X Z] = q (no shrinkage), so it lies between the n - p
and n - q corrections. For GLMMs W must be the *working* weights
(`fit@pp$Xwts`, which on this branch already include 1/phi; `getME(fit, "L")`
is factored with them), not the prior weights `hatvalues.merMod()` uses
(hence its "may not make sense for GLMMs" warning). edf depends on phi
through W, so it needs a fixed point.

Implementation options, cheapest first:

1. **Reporting only** (about 20 lines of R in `phi_alt()`): compute edf from
   the fitted model's `L`, `RX`, `RZX`, `Xwts` and report dev/(n - edf),
   iterating phi/edf with u fixed. Changes only the reported phi, not theta.
2. **Approximate in-fit** (R, no C++): refit with `qEff` set to the current
   edf (via the existing numeric `qEff` hook, `mkGlmerDevfun()` ->
   `rho$qEff`), recompute edf at the new estimates, repeat until it settles
   (warm starts; probably 2-4 refits). Affects theta. Approximation: edf is
   fixed at its value at the final estimates instead of varying with theta
   during the optimization (not checked).
3. **Proper C++** (a day or more): compute edf inside `profilePhi()` at each
   outer phi iteration. C_L needs a solve with L against a q x n right-hand
   side: trivial when L is diagonal (single scalar RE), but for large or
   crossed designs possibly much more than PIRLS itself, on every devfun
   call. Options: compute it once per devfun call after PIRLS converges
   (edf lags one phi iteration), or a stochastic (Hutchinson) trace
   estimator. Needs benchmarking on large models (the same scaling concern
   that moved `computeQEff()` to sparse QR).

Plan: try 1 and 2 in this study first. If the edf correction removes the
bias in the phi = 0.2, cluster size 4 cell without inflating the others, option 3 is
worth doing.

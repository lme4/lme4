# Caching more `glmer.Rnw` chunks: notes (2026-10-22)

Decision for 2.1-0: **do nothing**. The check time (excluding the
incoming-feasibility URL timeouts) is already below 2.0-6's, and the
size cost of storing results seemed not worth it. Notes kept in case we
revisit.

## Context

win-builder R-devel check of 2.1-0 (2026-10-22) took 643 s, over the
nominal 10-minute target. Of that, 136 s was the incoming-feasibility URL
check timing out on jstatsoft.org / lmdvr.R-Forge (URLs in
`inst/doc/lmer.pdf`). Remaining big items: tests 178 s (testthat 120 s,
since reduced by gating two slow new tests at `testLevel > 1`, commit
79aaae17), vignette re-building 111 s. For comparison, CRAN's own checks of
2.0-6 report 702 s (r-devel-windows) and 610 s (r-release-windows) check
time, status OK.

Locally, `glmer.Rnw` takes ~62 s to build (49 s R code, rest LaTeX);
`lmer.Rnw` ~24 s; the other vignettes 1-4 s each.

## Slowest `glmer.Rnw` chunks (local timings, lme4 built 2026-10-21)

| chunk | line | time | contents | shown? |
|---|---|---|---|---|
| `cm1` | 1371 | 12.0 s | 6 Contraception `glmer` fits (`cm1`-`cm6`) | code echoed |
| `cm_varcorr` | 1511 | 5.5 s | 3 fits with `diag()`/`cs()` covariance structures | hidden |
| `diagnose` | 980 | 4.8 s | `allFit(gm2)` for convergence messages | hidden |
| `allfit_cbpp1` | 1239 | 4.3 s | `allFit(gm1)` | code echoed |
| `preliminaries` | 123 | 3.3 s | package loading | - |
| `toenailplot` | 810 | 3.1 s | toenail fit + `zeta()` plot | hidden |
| `perf_chk_mod_cbpp` | 1098 | 2.9 s | `performance::check_model()` figure | hidden |
| `cm_cs` | 1521 | 2.4 s | `cs()` Contraception fit | code echoed |

Remaining 35 chunks are ~1 s or less each. (Timed by `knitr::purl()` and
evaluating each chunk in sequence.)

## Cost of storing results (xz-compressed `.rda`, full `merMod` objects)

| chunk(s) | objects | size | time saved |
|---|---|---|---|
| `cm1` | `cm1`-`cm6` | 627 KB | 12.0 s |
| `cm_varcorr` | `cm.diag`, `cm.homdiag`, `cm.homcs` | 339 KB | 5.5 s |
| `cm_cs` | `cm.cs` | 149 KB | 2.4 s |
| `diagnose` | `aa` | 107 KB | 4.8 s |
| `allfit_cbpp1` | `gm_all` | 123 KB | 4.3 s |
| all | | ~1.1 MB | ~29 s |

The source tarball is currently ~4.5 MB, so storing everything would take
it to ~5.6 MB.

## How it would be done

There is already a mechanism: `inst/vignette_data/` has batch scripts
(`glmer_Contraception_batch.R`, `glmer_cbpp_batch.R`) and a `Makefile`
producing `Contraception_batch.rda` / `cbpp_batch.rda`, loaded by the
`loadbatch` chunk of `glmer.Rnw`. Move the fits into the batch scripts,
make the echoed chunks `eval=FALSE` display chunks, and load the stored
objects in hidden chunks. Stored results must be regenerated whenever
fitting changes; the `allFit` results also depend on which optimizer
packages are installed.

Cheaper alternative: store only what the vignette prints (`anova()`
tables, `VarCorr()`, `logLik()`, ...) rather than whole fits, which carry
the model frame, data and factorization -- probably tens of KB instead of
~1 MB, but need to check which later chunks use `cm1`-`cm6` (predictions,
`pframe`, ...).

## CRAN size limits

From the CRAN Repository Policy, as tracked in
https://github.com/eddelbuettel/crp:

- Since at least svn r2908 (2013-09): "As a general rule, neither data nor
  documentation should exceed 5MB" (and documentation trimmed to a maximum
  of 5MB). No explicit limit on the tarball.
- Added in r5711 (archived 2023-07-25), still current: "Source package
  tarballs should if possible not exceed 10MB."

There is no (and was no) 5 MB limit on the whole tarball. `R CMD check`'s
"installed size is X Mb" NOTE triggers above 5 MB installed by default;
lme4 already exceeds that on Linux because of `libs/`.

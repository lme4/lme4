# `tolPwrss` default (`1e-7` vs `1e-9`): timing/correctness impact

(by Claude, 8 August 2026)

## Background

[GH #289](https://github.com/lme4/lme4/issues/289) ("minor Laplace
discrepancies with toenail data") reports that `glmer()`'s fit to the toenail
dataset disagrees with `glmmML`, `glmmADMB`, and SAS by more than expected for
a canonical-link (binomial) model, where lme4 is supposed to compute the exact
Laplace approximation rather than an approximation to it. The issue thread
traces this to the default PIRLS convergence tolerance, `tolPwrss = 1e-7`: the
default fit's log-likelihood is measurably below (about 0.008 units) what it
reaches under a tighter tolerance, and no convergence warning is raised.
Tightening `tolPwrss` (e.g. to `1e-9` or `1e-13`) reproduces the other
packages' values to several significant figures.

This was investigated as a side question while reviewing the `Gamma_GLMM`
branch: that branch does *not* touch `tolPwrss` or PIRLS convergence (its
commits are entirely about dispersion/sigma-dof estimation for free-dispersion
families, which don't apply to binomial's fixed dispersion), so a rebuilt
`Gamma_GLMM` reproduces the GH #289 discrepancy byte-for-byte. This note
records what changing the default `tolPwrss` from `1e-7` to `1e-9` would cost,
as a separate, orthogonal fix.

## Method

Built lme4 from source twice (`R CMD INSTALL --preclean`), once with the
default in `R/lmerControl.R`'s `glmerControl` (`tolPwrss=1e-7`, unchanged) and
once patched to `tolPwrss=1e-9`. Timed every GLMM-relevant file under both
builds:

- all of `tests/*.R` and `tests/testthat/test-*.R` that call `glmer()` or
  `glmer.nb()` (58 files)
- every `man/*.Rd` example that fits a non-Gaussian `glmer`/`glmer.nb` model,
  extracted via `tools::Rd2ex()` (19 files)

run with `LME4_TEST_LEVEL=4` so that every level-gated block in the suite
executes (the highest threshold actually used anywhere in the tree is
`testLevel() >= 3`; nothing is gated higher than that, so level 4 exercises
everything). Runs were sequential (tests, then examples, one setting at a
time) after an initial pass showed that letting the two sweeps overlap
introduced enough CPU contention to swamp the effect being measured. The
`R/lmerControl.R` edit was reverted and the package reinstalled at the
original default afterward.

## Timing result: no measurable effect

| | tests (58 files) + examples (19 files), 79 total |
|---|---|
| `tolPwrss=1e-7` | 770.8s |
| `tolPwrss=1e-9` | 773.9s |
| **delta** | **+3.2s (+0.4%)** |

Individual files swing up to ±9s in either direction (`VerbAgg` example
-9.3s, `test-gamma_glmm_bias.R` +8.0s, `test-nbinom.R` -4.1s, ...) with no
consistent sign. PIRLS iteration count is a small fraction of total runtime
for essentially every fit in the suite, dwarfed by simulation/bootstrap/
profiling loops that dominate wall time. **Tightening `tolPwrss` by two
orders of magnitude has no measurable aggregate timing cost.**

## Correctness result: not free

Three files that pass cleanly at `tolPwrss=1e-7` fail at `1e-9`:

- **`tests/testcrab.R`** — a genuine behavioral change, not just a tolerance
  break: a fit that converges cleanly at `1e-7` now throws `Model failed to
  converge with max|grad| = 0.0745227 (tol = 0.002, component 1)` at `1e-9`,
  and the resulting predictions differ from the `1e-7` fit by ~3.3%. This is
  the one that would need real attention before changing the default.
- **`tests/glmer-1.R`** — a `stopifnot(all.equal(...))` comparing two
  equivalent parameterizations of the same model (weights given as a matrix
  vs. as proportions) that agreed to <4e-5 at `1e-7` now differ by 3.2e-4 at
  `1e-9` — a reference-tolerance break that would need the test's tolerance
  loosened, not obviously a sign of a real problem.
- **`tests/predict_basis.R`** — trivial: two `predict()` calls differ by
  ~4e-6, tripping a tight default `all.equal`/`identical()` tolerance. Not
  meaningful.

## Conclusion

Moving the `glmerControl()` default `tolPwrss` from `1e-7` to `1e-9` would fix
the GH #289 discrepancy (canonical-link binomial fits reaching the true
Laplace optimum instead of stopping short of it with no warning) at
essentially zero timing cost across the existing GLMM test/example suite.
It is not risk-free, though: `testcrab.R` shows the tighter tolerance can
also destabilize a *different* fit's convergence, and at least one existing
test (`glmer-1.R`) would need its tolerance loosened to match. Before
changing the default it's worth deciding between:

- a more conservative target than `1e-9` (e.g. `1e-8`, which was already
  enough to reproduce the other packages' toenail values to 8 figures per the
  GH #289 thread), and/or
- pairing the tolerance change with an explicit PIRLS convergence check/
  warning, rather than silently changing the default and hoping affected fits
  still converge.

# Gamma GLMM random-effect variance bias (GH #643)

Investigation notes, 2026-07-28, comparing `lme4pureR`, `lme4` (dev HEAD,
v2.1-0), and `glmmTMB`. Originating issue:
<https://github.com/lme4/lme4/issues/643> ("Reliability of glmer to fit
Gamma-distributed model - Strange Ranef Variance"), most recently revived by
<https://github.com/TiagoAMarques/Report4BB>. Scratch scripts referenced below
live in the session scratchpad (not committed); `misc/GH643.R` in this repo
has the original exploratory session.

## tl;dr

`glmer` systematically **overestimates** random-effect variances for Gamma
(and presumably other estimated-dispersion-family) GLMMs, relative to
`glmmTMB`. Fixed effects and the dispersion parameter itself are essentially
unbiased. The bias grows as the true dispersion moves away from 1 and is
worst in the regime real Gamma GLMMs usually live in (small dispersion,
large shape). Root cause is understood precisely (see below): the PIRLS
working weights that drive both the mode-finding step and the Laplace
log-determinant are computed disp-free (as if phi=1 always), while only the
final reported deviance/logLik term is disp-aware. Two naive
post-hoc-formula patches (§5, §6) failed to fix it — one even made it much
worse. A working fix (§7) was found and validated by simulation: reweight
the PIRLS working weights by `1/phi` and profile phi via a nested
fixed-point loop inside PIRLS (Gamma-specific digamma-MLE equation) — a
real change to the inner fitting loop, not a small patch to the final
deviance formula, but it does *not* require any change to the outer
optimizer interface (theta/beta stay the only externally-visible
parameters). The theoretical vignette (`vignettes/glmer.Rnw`) has the same
gap in its derivation and needs a corresponding update — see §8.

## 1. Empirical confirmation

### Real data (Report4BB's `crglmer` model)

`crate ~ (1|location) + (1|fyear)`, Gamma(log), n=103:

| | glmer | glmmTMB |
|---|---|---|
| sd(location) | 0.308 | 0.112 |
| sd(fyear) | 0.250 | 0.063 |
| dispersion | 0.278 | 0.318 |
| logLik | -21.1 | -23.0 |

### Known-truth simulation (B=150, same design, true sdL=0.2, sdY=0.1, disp=0.12)

| | truth | glmer (% bias) | glmmTMB (% bias) |
|---|---|---|---|
| mean | 1.00 | -1.5% | -0.04% |
| dispersion | 0.12 | +0.3% | -0.4% |
| sd(location) | 0.20 | **+75%** | -28% |
| sd(year) | 0.10 | **+158%** | -27% |

glmmTMB's bias is small and *uniform* across the two REs (~-27% for both) —
the textbook downward bias of an ML (as opposed to REML) variance-component
estimator with few groups, unrelated to family. glmer's bias is much larger,
*asymmetric* between the two REs, and in the *opposite direction* — a
distinct extra mechanism specific to glmer.

### `glmmTMB::epil2` (58 subjects — enough groups to rule out small-sample ML
variance bias as the main driver)

`y ~ trt + (1|subject)`, Gamma(log): glmer sd(RE) = 1.18, glmmTMB sd(RE) =
0.79. AIC differs by ~59 points despite near-identical fixed effects.

## 2. nAGQ test: not an optimizer artifact

`glmer` only allows `nAGQ>1` for single-scalar-RE models. On `epil2`:

| nAGQ | sd(RE) | logLik |
|---|---|---|
| 1 (Laplace, default) | 1.179 | -553.8 |
| 5 | 0.670 | -72.6 |
| 15, 25 | 0.670 (same) | -72.6 (same) |

The logLik jump (-553.8 → -72.6) is silent — no convergence warning at
either end. Scanning the **Laplace (nAGQ=1) deviance surface directly** as a
function of theta confirms it has its own genuine minimum at theta≈1.18 —
this is not the optimizer getting stuck at a bad start; the
Laplace-approximated objective itself prefers the inflated value.

Caveat: repeating this with only 7 real grouping levels (Report4BB's
`location` factor) makes `nAGQ=15` collapse to singular (zero-variance) fits
in 88% of replicates — a well-known, unrelated small-number-of-groups ML
instability. `epil2`'s 58 groups avoid that confound, which is why it's the
cleaner diagnostic case.

## 3. Dispersion scan: pins the mechanism to phi, specifically

Single-RE simulation (30 groups × 20 reps, true sigma=0.3), `glmer`'s
default fit, across true dispersion values:

| disp | shape (1/disp) | bias in RE sd |
|---|---|---|
| 0.05 | 20 | +114% |
| 0.20 | 5 | +49% |
| 1.00 | 1 | **+1.6%** (essentially unbiased) |
| 3.00 | 0.33 | -76% (see update below — collapse regime, not noise) |
| 8.00 | 0.125 | +60% (small B=5-ish scan; see update below) |

Bias vanishes almost exactly at disp=1. Real Gamma GLMMs are almost always
disp≪1 (Report4BB: disp≈0.12-0.32; `epil2`: disp≈0.27-0.65) — exactly where
this bias is worst.

**Update, 2026-07-29:** the disp>1 (shape<1) regime is no longer
unexplained. A cleaner repeat at disp=2 (shape=0.5), same single-RE 30×20
design, B=100: mean bias **-67.2%**, with the *median* estimate exactly
0 — a majority of fits hit a singular (zero-variance) boundary solution.
So disp>1 isn't "erratic," it's a **different, consistent failure mode**:
instead of inflating theta (disp<1 regime), `glmer` systematically deflates
it, often collapsing to a singular fit entirely. The -76%/+60% single-run
numbers above were just noisy draws from this same collapse-prone regime,
not a separate unexplained phenomenon. Locked in as a characterization test
in `tests/testthat/test-gamma_glmm_bias.R` (gated behind
`LME4_TEST_LEVEL > 1`), alongside the disp=0.05 inflation case, so both
directions are covered.

**Further update, 2026-07-29 (later the same day):** turns out the §7 fix
*does* also resolve this side, even though it was only derived/validated
for disp<1 — see §7 for the full before/after comparison now that the fix
is implemented (not just prototyped) in `lme4pureR::pirls()`.

A first attempt to reproduce this in a 2-RE random-slope model (`y ~ 1 + x +
(1+x|f)`, disp=2, 20 groups) did *not* show a clean bias: at 500 obs/group
it looked essentially unbiased (enough per-group information apparently
swamps the effect), and at 50 obs/group it showed large bias but confounded
by frequent singular fits from having only 20 groups to estimate a 2×2
unstructured RE covariance (the same well-known few-groups instability
noted in §2, not obviously the same mechanism as the single-RE collapse
above). Isolating the single-RE case (this update) was needed to get a
clean, unconfounded read; revisiting the random-slope design is a
reasonable next step.

## 4. Mechanism, traced through code

Working from `~/R/pkgs/lme4pureR/R/pirls.R` (much easier to read than the
compiled internals; `lme4`'s C++ implements the same structure in
`RglmerWrkIter`/`glmerLaplace` in `R/lmer.R` and `src/predModule.cpp`).

The Laplace deviance actually optimized over theta (and beta) is:

```
Lm2ll <- aic(y, mu, wt, dev) + sum(u^2) + ldL2
```

(`pirls.R:122-125`; C++ analog `resp$resDev() + pp$sqrL(1)`, `R/lmer.R:423`).

- `aic()` (`src/glmFamily.cpp::gammaDist::aic`, `glmFamily.cpp:239-246`) computes
  `disp <- dev/nn` internally and plugs it into `dgamma(...)` — this term
  *is* dispersion-aware (using a moment-based `dev/n` plug-in for phi, not
  the true conditional MLE — see §5).
- `sum(u^2)` — the random-effect prior penalty — has **no** phi anywhere
  (`predModule.cpp:164`: `sqrL(f) = u(f).squaredNorm()`, literally `‖u‖²`).
  This is actually *correct* as-is: the standard GLMM spec (used by both
  `glmer` and `glmmTMB`) is `U ~ N(0,I)` fixed, not scaled by phi (unlike
  the LMM case — see below).
- `ldL2`, the log-determinant/curvature term that governs how much theta
  gets shrunk, comes from a Cholesky factor built from working weights
  `Whalf = sqrt(weights / variance(mu))`, and
  `gammaDist::variance(mu) = mu^2` (`glmFamily.cpp:269`) — **no phi**. This
  *should* have a `1/phi` factor (see §7): the true curvature of the
  Gamma conditional log-density in u scales with 1/phi, but lme4's working
  weights (both for the mode-finding Newton step and for the log-determinant)
  never include it.

glmmTMB has no such split: one joint TMB AD tape differentiates the full
log-density (beta, theta, phi together), so 1/phi correctly enters every
derivative, including the innermost Laplace curvature over u — this is why
"glmmTMB uses Laplace too" doesn't contradict the finding above; the gap
isn't Laplace vs. not-Laplace, it's phi-blind curvature vs. phi-aware
curvature.

### Vignette theory has the identical gap

`vignettes/glmer.Rnw` §3.3 ("The Generalized Linear Mixed Model") derives
the PSUD/Laplace deviance (eq. `condmodeGLMM`, `LaplaceDev`) using a bare
`sum(u^2)`/`||u||^2`, with no phi anywhere, and states explicitly (line
316-317) that "the scale parameter sigma is not present in the most common
GLMMs (binomial and Poisson families)" — i.e. the derivation is scoped to
disp-free families from the outset and never generalized.

Contrast with §3.2 ("The Linear Mixed Model", lines 244-296): there, `sigma^2`
*is* present, and critically it's **shared** — `U ~ N(0, sigma^2 I)`
(line 249), the *same* sigma^2 that scales the residual noise. That shared
scale is what makes closed-form exact profiling of sigma^2 possible (line
293: "the minimum doesn't depend on sigma^2"). The GLMM section copies the
"spherical random effect" structure from the LMM section but silently drops
this shared-scale linkage — fine for binomial/Poisson (no phi to share),
never revisited for Gamma/inverse Gaussian (which do have a phi, but it's
*not* meant to be shared with the RE prior the way LMM's sigma^2 is — see
§7; §5-6 chased that wrong analogy before §7 sorted it out).

(The commented-out footnote about `glmer.nb`/negative binomial near line
339-345 is unrelated: NB's dispersion/size parameter is not a standard
exponential-dispersion-family scale parameter and is handled by a
completely separate outer 1-D profile loop, not the PIRLS/Laplace machinery
discussed here.)

## 5. Attempted fix #1: naive phi-scaling of the u penalty — FALSIFIED

Direct analogy to the LMM's shared-scale trick suggests replacing
`sum(u^2)` with `sum(u^2)/phi_hat` (and originally, incorrectly, an
additional `-q*log(phi_hat)` correction to `ldL2` — algebraically this
term cancels against the RE-prior's own normalizing constant if `U ~
N(0,phi*I)` is assumed, so was dropped in the second attempt).

Tested in `pirls_dispfix.R` (reimplementing lme4's PIRLS/Laplace loop from
scratch via `glFormula`, validated to reproduce real `glmer` output almost
exactly at the same parameter values before trusting it — this reimplementation
also surfaced a **latent bug in lme4pureR's own reference `pirls()`**: it
calls `aic(y, rep.int(1,n), mu, weights, NULL)`, passing `dev=NULL`; this is
silently broken for Gamma (whose `aic()` needs `dev` to compute `disp`) even
though it happens to work for Poisson/binomial, whose `aic()` ignores the
`dev` argument entirely).

Result: **made the bias much worse**, not better, at every tested dispersion:

| disp | orig bias | "fixed" (u²/phi) bias |
|---|---|---|
| 0.05 | +117% | +359% |
| 0.20 | +43% | +105% |
| 1.00 | +1% | -2.5% |

## 6. Attempted fix #2: replace the dev/n plug-in with the true Gamma MLE — ALSO FALSIFIED

Hypothesis: maybe attempt #1's failure was because `phi_hat = dev/n` (a
crude moment estimator) is a poor stand-in for the true conditional MLE of
phi, and using a better phi_hat would fix it. The Gamma dispersion MLE has a
simple closed form: solve for shape `nu = 1/phi` in

```
log(nu) - digamma(nu) = deviance / (2n)
```

(standard result, e.g. underlying `MASS::gamma.shape`; LHS is monotonic
decreasing from +inf to 0, so `uniroot` finds a unique root). Verified this
recovers the true dispersion accurately on a large single sample (true
disp=0.05: moment estimate 0.0525, MLE estimate 0.0521).

Tested (`pirls_dispfix_mle.R`) substituting this MLE phi_hat everywhere the
moment-based one had been used (in both the `aic()` term and the `u^2/phi`
penalty):

| disp | orig bias | moment-plugin "fix" | **MLE-plugin "fix"** |
|---|---|---|---|
| 0.05 | +117% | +359% | +360% |
| 0.20 | +43% | +105% | +108% |
| 1.00 | +1% | -2.5% | -1.0% |

**No meaningful difference from attempt #1.** This rules out "phi_hat
wasn't accurate enough" as the explanation — the correction formula itself
is structurally wrong, independent of how phi is estimated. This is a
Gamma-specific test (the digamma equation above is derived from the Gamma
log-density); it would need separate re-derivation for other
estimated-dispersion families such as inverse Gaussian.

## 7. The working fix

Re-deriving once more: attempts #1-2 assumed (by direct analogy with the
LMM's shared-sigma^2 trick) that the RE prior should become `U ~ N(0,
phi*I)`. But that's the *wrong* target model — the standard GLMM
specification (used by both `glmer` and `glmmTMB`) is `U ~ N(0,I)` **fixed**,
independent of phi; theta alone parameterizes the actual RE covariance. Only
the *deviance/data-fit* term should be reweighted by `1/phi`, not the `‖u‖²`
prior penalty. That distinction matters a lot mechanically: a differential
reweighting (deviance term only) does **not** cancel out of the location of
the PIRLS mode the way a uniform reweighting of the whole bracket would —
meaning `u_hat` itself depends on phi under the correctly-specified model,
and the fix has to reach into the **inner PIRLS mode-finding loop and working
weights** (reweight `Whalf` by `1/phi`, leave the `Imult=1` prior-penalty
contribution to the precision matrix unscaled), not just patch the final
`Lm2ll`/`ldL2` formula.

**First prototype (`pirls_joint_phi.R`):** treat `phi` as a genuine free
parameter and optimize `(theta, log(phi), beta)` jointly via the same outer
`optim()` call (closer in spirit to glmmTMB's single joint AD-based ML).
This directly confirmed the math is right — bias essentially disappears:

| disp | theta bias (orig) | theta bias (joint theta/phi/beta ML) | phi bias |
|---|---|---|---|
| 0.05 | +117% | **-0.3%** | -1.4% |
| 0.20 | +43% | **-8.5%** | +0.8% |
| 1.00 | +1% | **-1.5%** | 0.0% |

**Second prototype, less invasive (`pirls_innerloop_phi.R`):** promoting phi
to a full outer-optimizer parameter is not actually necessary — it can be
profiled out via a fixed-point loop *nested inside* PIRLS instead, exactly
the same pattern `glmer` already uses to fold beta into the inner loop for
`nAGQ=0`. Concretely: for the current (theta, beta) from the (unchanged)
outer optimizer, alternate (a) run PIRLS to convergence for `u` given the
current phi estimate (using the phi-reweighted working weights), (b)
re-estimate phi via the digamma-MLE from the resulting fit, repeat to
convergence — all inside the devfun evaluation. **The outer optimizer
interface is completely unchanged** (still just theta, or theta+beta,
exactly as today); only the internals of one devfun evaluation grow an
extra nested convergence loop, no different in kind from how e.g. `nAGQ=0`
already nests beta/u profiling. Result is essentially identical to the
joint-optimization version:

| disp | theta bias (orig) | theta bias (inner-loop-profiled phi) |
|---|---|---|
| 0.05 | +117% | **-0.4%** |
| 0.20 | +43% | **-8.8%** |
| 1.00 | +1% | **-2.2%** |

So: the fix does *not* require promoting phi to an outer ML parameter or any
other "radical" change to the optimization architecture. It requires (a)
reweighting the PIRLS working weights by `1/phi` (deviance curvature only,
not the RE-prior penalty) and (b) profiling phi via a Gamma-specific
digamma-MLE fixed point nested inside the existing PIRLS convergence loop,
alongside the existing u-finding (and, for `nAGQ=0`, beta-finding)
iteration. This is a real change to `PIRLS`/`glmerLaplace`/`glmFamily.cpp`'s
Gamma path (not a one-line patch), but it's contained to the inner loop and
doesn't touch the outer optimizer, the `theta`/`beta` parameter vector, or
downstream code (summary, profiling, confint, etc.) that consumes the
current interface.

**Caveat / scope:** the digamma-MLE step (`log(nu) - digamma(nu) =
dev/(2n)`) is Gamma-specific, derived from the Gamma log-density. Other
estimated-dispersion families available in `glmer` (chiefly inverse
Gaussian) would need the analogous MLE equation re-derived from their own
density before the same nested-fixed-point structure could be reused; this
has not been attempted here. The `1/phi`-reweighting-of-working-weights
half of the fix, in contrast, is generic (follows directly from the general
exponential-dispersion-family deviance decomposition), only the specific
root-finding equation for phi is family-specific.

All four prototypes above were validated against real `glmer` output before
being trusted: the unmodified reimplementation matches real `glmer`'s fitted
theta/beta closely on a spot-checked dataset (0.6351 vs. 0.6346), confirming
the reimplementation faithfully mirrors lme4's actual PIRLS/Laplace
formula and isn't just a coincidentally-different computation.

### Implemented in `lme4pureR::pirls()`, and fixes both bias regimes (2026-07-29)

The standalone-script prototypes above were superseded by an actual
implementation in `lme4pureR`'s own `pirls()` (not just a from-scratch
reimplementation): a new `phiType` argument, `c("none", "moment",
"digamma")`. `"none"` (the default) reproduces the original, unmodified
disp-blind behaviour exactly (verified: all pre-existing `lme4pureR` tests
still pass byte-for-byte unchanged with this default). `"moment"` and
`"digamma"` both apply the nested-fixed-point fix from above, differing
only in how phi is re-estimated each outer iteration (`dev/n` vs. the
digamma-MLE).

Testing this real implementation (not just the standalone prototype)
against *both* known bias regimes — including the disp>1/collapse regime
from the update above, which the fix had never actually been tried against
— gave a bigger result than expected: **it fixes both**, not just the
disp<1/inflation side it was derived for.

| | shape=20 (disp=0.05, inflation) | shape=0.5 (disp=2, collapse) |
|---|---|---|
| `phiType="none"` | +117.0% bias | -94.5% bias, 97% of fits collapse to ~0 |
| `phiType="moment"` | **-0.1%** bias | **+4.6%** bias, 0% collapse |
| `phiType="digamma"` | **-0.1%** bias | -16.0% bias, 0% collapse |

(B=30 replicates per cell, single-RE 30×20 design as in §3.) `moment` and
`digamma` are essentially tied in the inflation regime, but `moment`
noticeably *beats* `digamma` in the collapse regime — the crude `dev/n`
plug-in appears to be fully adequate for correcting the RE-variance bias in
both directions, making the more complex Gamma-specific digamma-MLE
possibly unnecessary for this purpose (it may still be worth keeping
around for more accurate dispersion *reporting*, just not required for the
theta-bias fix itself). This weakens the scope caveat below: the
family-specific piece of the fix (the digamma equation) turns out not to be
the load-bearing part.

Regression tests locking in this comparison:
`lme4pureR`'s `inst/tinytest/test_pirls_phiType.R`.

## 8. Corresponding changes needed in `vignettes/glmer.Rnw`

The theoretical derivation has the identical gap as the code (§4), and
would need updating alongside any code fix:

- **§3.2 "The Linear Mixed Model" (lines 244-296)** is fine as-is and is the
  right model for the fix: it already shows the *correct* pattern for a
  family with a free scale parameter — sigma^2 is shared between the
  residual noise and the "spherical" random-effect prior,
  `U ~ N(0, sigma^2 I)` (line 249), which is exactly why sigma^2 profiles
  out in closed form (line 293). Nothing to change here; it's the reference
  case.

- **§3.3 "The Generalized Linear Mixed Model" (lines 310-366)**: line
  316-317 currently states "the scale parameter, sigma, is not present in
  the most common GLMMs (binomial and Poisson families)" and then never
  revisits the estimated-dispersion case. Needs a new sentence/paragraph
  noting that for families *with* an estimated dispersion phi (Gamma,
  inverse Gaussian), phi enters differently than LMM's sigma^2 does: it
  reweights the unit-deviance term but, unlike sigma^2, does *not* rescale
  the random-effect prior (`U ~ N(0,I)` stays fixed regardless of phi —
  see finding in §7).

- **Equation `condmodeGLMM` (line 324)**: currently
  `sum_i d(y_i,mu_i) + ||u||^2`. For families with an estimated dispersion,
  this needs to become `sum_i d(y_i,mu_i)/phi + ||u||^2` (deviance term
  reweighted by `1/phi`, `||u||^2` unchanged) — matching §7's finding that
  `u_hat` genuinely depends on phi once this is written correctly, unlike
  the current (implicitly phi=1) formula.

- **Equation `LaplaceDev` (lines 504-508)**, the main formula actually
  implemented (`Lm2ll` in the code): currently
  `sum_i d(y_i,mu_i) + ||u~||^2 + log(|L|^2) + (q/2)log(2*pi)`. Needs the
  same `/phi` on the deviance term as above, i.e.
  `sum_i d(y_i,mu_i)/phi + ||u~||^2 + log(|L|^2) + (q/2)log(2*pi)`
  (`||u~||^2` and the log-determinant term are unaffected — see §7's
  derivation for why the naive "also scale the log-determinant by phi"
  correction turns out to be wrong/unnecessary once the u-prior is
  correctly treated as phi-independent). Should also note, since this
  equation is now phi-dependent, that phi has to be determined
  concurrently with `u~` (a nested profiling step, §7), rather than only
  after the fact as `sigmaML`/`disp = dev/n` currently is.

- A new remark (could go near the commented-out `glmer.nb` footnote around
  line 339-345, or as a new paragraph after eq. `LaplaceDev`) should
  explain that the *current* code and vignette both implicitly assume
  phi=1 in these two equations, which is exact for binomial/Poisson
  (hence no bug there) but is the source of the systematic random-effect
  variance inflation documented in this file for Gamma GLMMs.

None of this has been applied to `glmer.Rnw` yet — recorded here as the
spec for that edit once/if the code fix (§7) is actually implemented and
merged, so the vignette and the code change together rather than the
vignette silently going stale again.

## 9. Implemented in the compiled path (branch `Gamma_GLMM`), 2026-07-29

Per user request, the "moment" approach (only — no digamma, no opt-out) is
now implemented directly in lme4's compiled fitting path, unconditionally
for the Gamma family:

- `src/respModule.h`/`.cpp`: `glmResp` gets a `d_phi` member (default 1.0,
  reproduces original behavior when unused), `phi()`/`setPhi()` accessors,
  and `sqrtWrkWt()` divides the deviance/curvature term by `d_phi`.
- `src/external.cpp`: `internal_glmerWrkIter`'s PWRSS convergence criterion
  divides `resDev()` by `phi()` (keeps step-halving consistent with the
  reweighted objective; `sqrL`/`||u||^2` stays unscaled, per §7).
  `glmerLaplace()` — the entry point used for nAGQ=0/1 fits, i.e. nearly
  all real GLMM fits (`glmerAGQ`, the separate nAGQ>1 single-scalar-RE
  path, was deliberately left untouched) — wraps the single `pwrssUpdate()`
  call in an outer fixed-point loop over phi, gated on
  `rp->family() == "Gamma"` only.

**Validated working**, matching the lme4pureR R-level `"moment"` prototype
closely on real `glmer()` calls (single-RE, B=30 batches):

| | shape=20 (disp=0.05, inflation) | shape=0.5 (disp=2, collapse) |
|---|---|---|
| before | +117.0% bias | -94.5% bias, 97% singular |
| after | -0.1% bias | +4.4% bias, 0% singular |

Full test suite (`devtools::test()`, `LME4_TEST_LEVEL=2` — note
`testthat::test_dir()` after a plain `library(lme4)`/`load_all()` gives
false-positive "object not found" errors for internal/non-exported
objects like `glmFamily`; `devtools::test()` is the correct invocation)
passes cleanly except three now-mechanical hardcoded-reference-value
updates (`test-glmFamily.R:155`, `test-methods.R:659`,
`test-gamma_glmm_bias.R` — the latter rewritten from "characterize the
known bug" to "verify the fix stays working," since it's now unconditional)
and one open item, below.

### Robustness investigation: `test-isSingular.R`'s "checking singular fit for merMod"

This test (a genuine illustration of singularity detection, two crossed
random-slope terms, one deliberately singular in truth, Gamma shape=0.25)
started throwing `"(maxstephalfit) PIRLS step-halvings failed..."` under
the fix — previously it converged (to a biased-but-valid fit). Tried,
in order: 50/50 and then gentler (90/10, log-space) damping of the phi
update; a try/catch fallback resetting `delu`/`delb` to zero and retrying
at phi=1; retrying at `lastGoodPhi` instead of phi=1; returning the
partial (non-converged) Laplace value instead of a flat sentinel on
double failure (this one made things *worse* — a misleadingly low partial
value pulled the optimizer toward a degenerate boundary rather than away).
None of these change the final outcome — confirmed by explicit tracing,
not just the end result: **245 of 246 internal failures occur at
`outer=0`, `phi=1`** — i.e. at the exact computation the original,
unmodified code always did. So this isn't new fragility introduced by the
phi fix.

Two further diagnostics (reprex scripts `misc/GH643_reprex_pirls_fragility.R`
and `misc/GH643_reliability_comparison.R`, run against lme4 2.0-6 installed
into an isolated library via `remotes::install_version` alongside the
current dev build, so both can be compared without clobbering either):

1. **Near the true parameters, both versions are equally robust.** 200
   evaluations of the joint (theta,beta) Laplace deviance at points
   perturbed ±10% around the true (deliberately singular) structure: 0
   errors for *either* version, nearly identical deviance values.
2. **The actual failing candidates are not near the truth.** Tracing the
   real fit showed candidates like `theta=(1,0,1,1,0,1)` (literal starting
   values) recurring for nearly the whole optimization, and — critically —
   `beta` itself was `NaN` in many failing evaluations. Directly
   re-evaluating one such (valid, non-NaN) point fresh against lme4 2.0-6
   **succeeds** (deviance 1813.0), which only makes sense if the failures
   in the live fit are driven by *warm-started internal state*, not the
   candidate point itself.

**Mechanism:** the fix's altered deviance surface routes `bobyqa` through
a few early failing evaluations. Once several evaluations in a row return
the same flat sentinel (`1e10`), `bobyqa`'s own trust-region/quadratic
interpolation model — a known failure mode for derivative-free optimizers
fed a locally-flat objective — degrades and starts proposing **NaN
candidates**. lme4's R-level optimizer wrapper commits whatever candidate
it's given as `pp`'s new baseline (`beta0`) without checking finiteness
first (every devfun call still "succeeds" from its point of view, since it
gets back *some* finite number). Once `beta0` itself is NaN, every
subsequent evaluation inherits a NaN-poisoned baseline and fails too,
regardless of what `(theta,beta)` `bobyqa` asks for next — a
self-reinforcing cascade, matching the near-total (245/246) failure rate.

This is a **pre-existing gap in lme4's optimizer/devfun contract**
(candidates are never validated for finiteness before being committed as
the persistent baseline) that the original, unfixed code never had to
confront, because it had no try/catch at all — a single PIRLS failure
propagated immediately as an R error and aborted the fit, rather than
continuing and accumulating corruption. The Gamma dispersion fix is simply
the first thing to route the optimizer somewhere this latent gap gets
triggered. Fixing *that* gap (e.g. validating/rejecting non-finite
candidates before they're committed) is a separate, pre-existing lme4
robustness issue, out of scope for this fix.

**Reliability comparison, completed** (`misc/GH643_reliability_comparison_params.R`,
run overnight 2026-07-29→30, parallelized with `parallel::mclapply` at
`mc.cores=14` — machine has 32 cores, not 4 as first assumed; verified via
`/proc/cpuinfo` and `parallel::detectCores()`): B=100 *resimulated*
datasets (same generative process as `test-isSingular.R`'s near-singular
scenario — 20×20 balanced design, two crossed random-slope terms, group1
truly singular (`theta=0,0,0`), group2 not (`theta=1.0,0.5,0.3`),
Gamma(link=log), shape=0.25 — fresh random draws each time, shared seed
sequence so both lme4 versions see identical datasets), with the actual
fitted `theta`/`beta`/`sigma` captured for every non-error replicate (not
just pass/fail), so "clean convergence" can be checked for whether it's
actually *correct*, not just unflagged:

| | clean | warnings | **error** | runtime (100 fits, 14 cores) |
|---|---|---|---|---|
| lme4 2.0-6 (unmodified) | 85/100 | 3/100 | **12/100 (12%)** | 0.21 min |
| current (Gamma_GLMM fix) | 28/100 | **71/100** | **1/100 (1%)** | 4.0 min |

The fix essentially eliminates outright crashes (12%→1%) but converts most
of that into flagged convergence warnings (3%→71%) rather than clean fits
— consistent with the mechanism traced above (the nested phi loop hitting
real optimizer difficulty far more often, but now surfacing it via
warnings instead of an R error).

The more consequential finding is what happens *inside* "clean
convergence" — fits that report no problem at all. Group1's diagonal theta
(should be ~0, true singular structure) among clean fits:

| | old, clean (n=85) | current, clean (n=28) |
|---|---|---|
| group1 θ mean / median | 0.109 / 0.000 | **0.466 / 0.495** |
| clean fits with group1 θ>0.5 (silently wrong) | 12/85 (14%) | **13/28 (46%)** |

So among fits the *current* version reports as fully clean (no warning at
all), nearly half have silently converged to a clearly non-singular
estimate for the block that is truly singular by construction — a
materially higher silent-failure rate than old lme4's 14% (95% CIs
roughly 8–23% vs 28–65% — not just small-sample noise given n=28 vs 85).
Fixed effects (`beta`, `sigma`) are comparable between versions and
reasonably close to truth in both (medians ~1.8-2.0 and ~1.8, true 2 and
2, for old and current respectively).

**Net read:** the fix trades crashes for warnings (good — a warned,
inspectable fit beats an aborted one) but appears to somewhat *increase*
the rate of confidently-wrong, unflagged fits specifically for detecting
a truly-zero random-slope variance in one block while another block in
the same model is genuinely non-zero. This wasn't visible in the
single-random-effect bias-correction validation earlier in this document
(§3, §9) because that used a single, non-singular random effect throughout
— this is a distinct effect specific to mixed singular/non-singular
multi-term structures, not yet root-caused. Worth a closer look at whether
it connects to the same phi-fallback/optimizer-state mechanism diagnosed
above, or is a separate phenomenon, before deciding whether/how it should
change the fix.

## 10. "Clean but wrong" is genuine multimodality, not an optimizer miss

Follow-up to §9's finding above (`misc/GH643_multimodal_check.R`): is the
elevated "clean but wrong" rate under the fix (group1's truly-zero
variance silently converging to a substantial non-zero estimate) an
optimizer-convergence problem — `bobyqa` starting somewhere unlucky and
failing to find the true-adjacent optimum — or a genuine second mode in
the fixed objective? Tested on three of the affected datasets (reps 14,
26, 48 from the B=100 comparison in §9) two ways:

1. **Refit starting exactly at the true generating parameters**
   (`start = list(theta = c(0,0,0,1,0.5,0.3), fixef = 2)`) instead of
   lme4's default starting values. Result: converges back to essentially
   the *same* non-singular solution every time — bit-for-bit identical
   theta for reps 26 and 48, and within noise for rep 14. So this is not
   "the optimizer wandered away from a good start" — even starting
   exactly at the truth, it moves away from it.

2. **Directly evaluate the (unoptimized) Laplace deviance** at the
   default-converged point vs. at the true parameters exactly, via
   `updateGlmerDevfun(..., nAGQ=1)`. The non-singular point has **lower**
   deviance (better fit) in all three cases (e.g. rep 14: 1359.9 vs
   1363.5; rep 26: 1447.2 vs 1450.5; rep 48: 1540.2 vs 1545.9).

**Conclusion: this is genuine multimodality, not an optimizer failure.**
Under the fix's phi-reweighted objective, the non-singular solution
really is the better (likely global) optimum for these datasets — the fix
has shifted where the mode of the approximate likelihood sits, and not in
a way that only affects hard-to-reach corners of parameter space.

**Mechanism, corrected:** an initial hypothesis here was that phi might
somehow be mis-scoped between the two random-effect terms (e.g. "a phi
calibrated mostly by group2's scale doesn't correctly balance group1's
boundary case"). That doesn't survive scrutiny, and shouldn't have been
proposed: **the model as specified has exactly one conditional response
distribution and one shape/dispersion parameter, shared by every
observation regardless of which random-effect term applies to it** —
there is no such thing as a separate "phi for group1" vs "phi for
group2" to miscalibrate; phi is not decomposable by term at all. Checked
directly: at both the non-singular ("wrong") mode and starting from the
true parameters, the internally-profiled phi converges to essentially the
identical value (3.3152 vs 3.3149) — confirming the fix's single, jointly
consistent phi is doing exactly what it's supposed to at both points, and
ruling out the per-term-miscalibration idea entirely.

So the real, still-open question is why the (correctly single-phi)
profile-over-theta objective has two modes at all for this design — one
at the true, singular structure and one substantially away from it, with
the non-singular one deeper. This is a property of the *profile
likelihood surface itself*, not of anything term-specific about phi.
Plausible framings, not yet distinguished: (a) this is a genuine feature
of Gamma GLMM likelihoods near a variance-component boundary when another
term in the same model carries substantial signal — a well-known general
difficulty (profile likelihoods for near-zero variance components are
often irregular near the boundary), which the fix's more accurate
handling of phi could be *revealing* rather than *causing*; or (b) it's
specific to the Laplace approximation itself (with phi now handled
correctly) being a poor approximation to the true/AGQ marginal likelihood
in this regime, in a way the previous disp-blind version's own bias
happened to mask (old lme4's own "clean but wrong" rate was 14%, not
0% — so some version of this already existed, just less often). This is
a distinct effect from the NaN-cascade/optimizer-robustness issue in §9;
it doesn't involve any failed `pwrssUpdate` call or fallback path at all,
just the successfully-converged objective having a genuine second optimum.

**Still open:** distinguish (a) from (b) above — e.g. compare against
`glmmTMB` (a genuinely joint single-phi ML fit via full automatic
differentiation, not lme4's profile/PIRLS-based approach, though it also
uses a Laplace approximation for the random-effect integral) on the same
datasets: if glmmTMB also finds the non-singular point as its optimum,
that points to (a) — a real feature of this design's likelihood, not an
lme4-specific artifact; if it reliably recovers the near-zero solution,
that points to (b) — something more specific to lme4's implementation.
Also still open: decide how to adjust `test-isSingular.R`'s expectation
given both §9 and §10 — the user does not consider that test
"adversarial"/contrived (it illustrates genuine singularity), so any
adjustment should be a real accepted-limitation note, not a brushed-off
tolerance loosening; and consider whether §10's finding changes the
recommended scope of the fix itself (e.g. whether it should be restricted
to single-random-effect-term models until the mixed-term case is better
understood). The C++ changes (`src/respModule.h`, `src/respModule.cpp`,
`src/external.cpp`) are committed and pushed to `Gamma_GLMM` (commit
`99e015f2`).

## 11. Multistart diagnostic (Raue et al. 2013): real but modest multimodality, and a correction along the way

Follow-up to §10, using the method from Raue, A., Schilling, M., Bachmann,
J., et al. (2013), "Lessons Learned from Quantitative Dynamical Modeling in
Systems Biology," *PLOS ONE* 8(9):e74335: refit the *same* dataset (the
rep=14 target from §10, seed 9014) from B=200 random starting points drawn
from a broad hypercube, and plot the ECDF of the achieved -2*logLik at
convergence. A single dominant plateau means the optimizer reliably finds
one (presumably global) optimum regardless of starting point; multiple
separated plateaus mean multiple genuinely distinct local optima exist —
real multimodality, not an artifact of any one starting-value choice.

Run for **lme4 2.0-6**, **current (dev, fix)**, **glmmTMB**, and (added
below) **lme4 joint-phi**, all on the identical dataset, same random seed
(20260730) generating the starting-point ensemble in each script (not
numerically identical starting values across engines, which have different
native theta parameterizations). Hypercubes: lme4 theta diagonal ~U(0,2),
off-diagonal ~U(-1,1) (its native raw-Cholesky units); glmmTMB theta ~U(-2,2)
(its native log-SD/unconstrained-corr units, default-initialized at 0);
intercept `beta` ~U(0,4) for all (directly comparable — same units, the
log-link intercept). Parallelized with `parallel::mclapply`, 10 cores each
for the first three (30 of 32 cores simultaneously), 28 cores for the
joint-phi arm (run afterward, alone). Scripts: `misc/GH643_multistart_lme4.R`,
`misc/GH643_multistart_glmmTMB.R`, `misc/GH643_multistart_jointphi.R`,
`misc/GH643_multistart_analysis.R`; plots: `misc/GH643_multistart_ecdf.png`,
`misc/GH643_multistart_dispersion_hist.png`.

### Correction: `deviance()` is not `-2*logLik()` for these fits

The first pass of this analysis used `deviance(fit)` and found dramatic,
widely-separated plateaus (hundreds of units apart) for both lme4 versions.
That turned out to be measuring largely the wrong thing. Flagged directly:
for a fixed dataset, `deviance()` and `-2*logLik()` should differ by at
most a parameter-independent constant — but the observed gap varied by
over 1000 across fits of the *same* data, which isn't possible if that were
true. Traced to ground: **`deviance(fit)` for these Gamma GLMM fits is
exactly the raw, phi-free unit-deviance sum**, `sum(Gamma()$dev.resids(y,
mu, wt))` — verified identical to machine precision — with no `1/phi`
scaling and no phi-dependent normalizing constant, and no contribution from
the random-effects integration term (`‖u‖² + log|L|`) that `-2*logLik()`
correctly includes. Since phi is *estimated* and differs across fits (each
local optimum has its own fitted dispersion), the two quantities aren't
simply offset — `deviance()` just isn't a comparable likelihood-based
quantity across different fits of a dispersion-estimated GLMM.

Independently confirmed this is the conventional, expected role of the
`dev.resids()`/`wt` argument (not something specific to lme4): surveyed
every call to `dev.resids(...)` with a Gamma or gaussian family across all
R packages under `~/R/pkgs` (44 files). `wt` is *always* ordinary prior/case
weights, matching base R's `glm.fit()` convention — never used to carry
dispersion. `mgcv/R/efam.r:22` is the clearest confirmation: it computes
`sum(family$dev.resids(y, mu, wt, theta))/scale`, applying the dispersion
scaling as an explicit separate step after the raw call — exactly the
`/phi` scaling `deviance()` is missing here. This isn't a bug introduced by
the fix; both lme4 2.0-6 and current show the same `deviance()`/`-2*logLik()`
divergence, so it looks like a longstanding property of how `merMod`
defines `deviance()` for GLMMs with an estimated dispersion generally.

All results below use the corrected quantity, `-2*logLik(fit)` for lme4
(field `neg2ll` in the updated scripts) and `2*fit$obj$fn(fit$fit$par)` for
glmmTMB (already correct — matches `-2*logLik` exactly whenever the Hessian
is PD, and per glmmTMB's troubleshooting vignette still gives a usable NLL
when it isn't).

### Also found while rerunning: a small silent-degenerate-fit issue in the fix

2 of 200 lme4-current fits hit the C++ fix's flat `1e10` sentinel (§9 —
returned when even the safe phi=1 retry inside `glmerLaplace()` fails) but
did *not* trigger an R-level error, only an ordinary convergence warning
indistinguishable in kind from the other 37 warned-but-fine fits. Both are
recoded as "degenerate" and excluded from the analysis below; this is a
loose end in the §9 fallback (a fit that internally gave up should probably
be more clearly flagged than a generic convergence warning) worth tightening
separately, not folded into the multimodality question here.

### Results (corrected)

| | clean | warning | degenerate | error | dominant -2logLik plateau (count) |
|---|---|---|---|---|---|
| lme4 2.0-6 | 156 | 15 | — | 29 | **1393 (118)**, 1396 (32), 1385 (16), + small tail to ~1500 |
| lme4 current (fix) | 161 | 37 | 2 | 0 | **1360 (131)**, 1393 (18), 1396 (17), + longer tail to ~1500 |
| lme4 joint-phi | 200 | 0 | — | 0 | **1329 (158)**, 1332 (24), 1330 (2), + short tail to ~2745 |
| glmmTMB | 0 | 200 | — | 0 | **1326 (190)**, 1328 (5), 1330 (2) |

(lme4 joint-phi is discussed separately below — see "Fourth arm" — since it
uses a different devfun/optimizer setup than the other two lme4 rows, not
just a different profiling rule.)

The picture is much less dramatic than the `deviance()`-based first pass,
but a real, qualitative difference survives: **all three still show more
than one plateau**, but the plateaus for both lme4 versions are now only
tens of units apart (not hundreds) — still likely-meaningful separation on
a -2logLik scale (roughly matching a chi-square-type threshold, not just
optimizer tolerance noise), but nowhere near as severe as first reported.
glmmTMB's spread (1326/1328/1330, ~4 units total) is small enough to be
ordinary numerical tolerance around one optimum, not a second mode.

Visually (`misc/GH643_multistart_ecdf.png`): lme4-old rises in one fairly
sharp jump near its own dominant plateau with a short tail. **lme4-current
also jumps early to its (better/lower) dominant plateau, but then has a
visibly longer, more gradual staircase tail out to ~1500** — i.e. the fix's
dominant mode is objectively a better fit than old lme4's dominant mode,
but the fix's optimization landscape has *more* reachable alternative local
optima than old lme4's, not fewer. That's a genuine, previously-unstated
cost of the fix worth being aware of, distinct from the boundary-collapse
question in §9-§10.

### Dispersion (phi = sigma²), true value 4

`misc/GH643_multistart_dispersion_hist.png`:

| | median phi | range | shape |
|---|---|---|---|
| lme4 2.0-6 | 3.510 | 2.90–5.65 | dominant cluster ~3.5 (below truth), small secondary cluster ~5.6 |
| lme4 current (fix) | 3.315 | 2.82–5.65 | dominant cluster ~3.3-3.4 (below truth), longer spread-out secondary cluster ~4.5-5.0 |
| glmmTMB | 4.009 | 4.01–4.08 | essentially a single spike at the true value |

glmmTMB recovers the true dispersion almost exactly regardless of starting
point. Both lme4 versions underestimate phi somewhat at their dominant
mode (~3.3-3.5 vs true 4), with the fix's secondary cluster more spread
out than old lme4's — consistent with the longer ECDF tail above. This is
a different, more complex scenario (two random-effect terms, one singular)
than the single-RE cases validated as unbiased earlier in this document
(§3, §9), so it doesn't contradict that validation, but it's a reminder
the earlier validation doesn't automatically extend to this mixed-term
design.

### Pairs plots: the plateaus are genuine alternative parameter combinations

Every fitted `theta` (full 6-component native parameterization), `beta`,
and `sigma` was saved for all 200 starts × 3 methods, not just the
objective value. Pairs plots (`gap=0`) of all 8 values jointly, one per
method (`misc/GH643_multistart_pairs_lme4old.png`,
`_lme4current.png`, `_glmmTMB.png`; script `misc/GH643_multistart_pairs.R`):

- **Both lme4 versions** show visibly discrete clustering across theta
  *and* beta *and* sigma jointly — e.g. current lme4's `theta1` (group1's
  first diagonal entry) splits into a near-zero cluster and a ~0.6-0.8
  cluster, with `beta` and `sigma` shifting together depending on which
  `theta1` cluster a given start landed in. Confirms the ECDF plateaus
  are genuine alternative *parameter combinations* co-varying together,
  not just different objective values reached with similar parameters.
  Current lme4 shows visibly more cluster variety than old lme4, matching
  its longer ECDF tail.

- **glmmTMB** looks qualitatively different: `theta1`/`theta2` (group1's
  log-SD parameters, the truly-singular block) sit at strongly negative
  values (SD ≈ exp(-7) to exp(-1), i.e. essentially zero) across nearly
  all 200 starts — correctly and consistently recovering the singular
  structure regardless of starting point. `beta` and `sigma` are almost
  perfectly constant (beta ≈ 1.925-1.93, sigma ≈ 2.005-2.02) throughout,
  completely unaffected by exactly where group1's (functionally
  irrelevant, near-zero-variance) parameters land. `theta3` (group1's
  correlation transform) fans out to large values for some fits — the
  usual artifact of a correlation parameter becoming ill-defined as its
  corresponding variance shrinks to the boundary, and harmless here since
  a corner-case correlation on an ~zero-variance term doesn't affect the
  fit.

So the glmmTMB/lme4 contrast isn't just "tighter in NLL" — glmmTMB is
actually finding the *correct* (near-singular) structure consistently,
while both lme4 versions genuinely land on distinct, non-singular
alternative solutions with fixed effects and dispersion shifting together
with them.

### Fourth arm: phi as a first-class optimized parameter, not a nested fixed point

The working fix (§7-§9) profiles phi via a fixed-point iteration nested
*inside* PIRLS, at each trial value of theta/beta the outer optimizer
proposes. An alternative, tried earlier in this investigation as a
prototype (`misc/GH643_pirls_joint_phi.R`, §7) before the C++ fix existed:
put `log(phi)` directly in the parameter vector the *outer* nonlinear
optimizer sees, alongside theta and beta, and let it get optimized jointly
rather than profiled out at each outer step. The hypothesis: fixed-point
profiling adds an extra layer of (possibly imperfectly converged, possibly
locally-nonsmooth) nested optimization that the outer optimizer has to see
through, and removing that nesting might make the overall likelihood
surface better-behaved and more reliably optimized — matching, or getting
closer to, glmmTMB's TMB-style joint optimization over all parameters
(including dispersion) at once.

Reran that prototype on the *same* target dataset and starting-point
hypercube as the three arms above, switching its optimizer from
`optim(method="Nelder-Mead")` to `minqa::bobyqa` (lme4's own default
optimizer, for a fairer comparison), with explicit finite bounds
(`log(phi)` bounded to `[log(1e-4), log(1e4)]`). Its devfun computes
`-2*logLik` directly for a *given* phi (`aic_term(phi) + ‖u‖² + log|L|`,
no internal phi re-estimation at all) — so unlike the real lme4 fits above,
there's no `deviance()`-vs-`-2*logLik()` ambiguity to correct for here; the
objective value returned by `bobyqa` *is* the comparable quantity by
construction. Script: `misc/GH643_multistart_jointphi.R`.

**Results: this comparison confirms the hypothesis.** All 200/200 starts
converged "clean" (no errors, warnings, or sentinel hits — contrast both
lme4 fixed-point variants above, which had 15-37 warnings and, for 2.0-6,
29 outright errors). The dominant plateau is at **-2logLik ≈ 1329
(158/200 = 79%)** — much closer to glmmTMB's dominant plateau (1326,
190/200) than either fixed-point lme4 variant's dominant plateau (2.0-6:
1393, 118/200; current fix: 1360, 131/200). Dispersion recovery
(`misc/GH643_multistart_dispersion_hist.png`) is essentially unbiased at
the mode — median phi = 4.014 against a true value of 4, matching
glmmTMB's 4.009 — versus both fixed-point variants' persistent downward
bias (2.0-6 median 3.51, current-fix median 3.32). The pairs plot
(`misc/GH643_multistart_pairs_jointphi.png`) shows a tight dominant cluster
(theta1 ≈ 0.5-0.6, theta4 ≈ 0.4, beta ≈ 1.9, sigma ≈ 2.0) with only a
modest scatter of alternative-mode strays — visibly tighter than either
fixed-point lme4 pairs plot, though not as uniformly tight as glmmTMB's.

Caveat on scope: this is a standalone ~180-line R prototype validating the
*idea*, not a drop-in replacement for the C++ fix. It reimplements PIRLS
from scratch around a `phi`-parameterized weight matrix and calls
`bobyqa` directly; it isn't wired into `mkGlmerDevfun`/`optimizeGlmer` or
lme4's actual control-flow/optimizer-selection machinery, doesn't handle
multiple families or non-Gamma links, and its convergence/bounds choices
were hand-tuned for this one benchmark. Moving phi into the real outer
parameter vector for `glmer()` would be a materially larger change (widens
the parameter block the C++ optimizer interface exposes, changes what
`theta` means to callers, touches `refitML`, profiling, and starting-value
logic) than the nested-fixed-point fix already implemented — worth
pursuing given how much better it performs here, but a separate,
larger-scope follow-up, not a same-session swap-in.

### Scope/caveats

Only one dataset has been tested this deeply. Worth checking whether the
degree of multimodality tracks how close the true structure is to the
singular boundary, or is present more generally for Gamma GLMMs with
multiple random-effect terms regardless of singularity. Also worth
checking whether glmmTMB's much tighter clustering (in both -2logLik and
dispersion) holds up on other datasets, or is specific to this one, and
whether the joint-phi parameterization's advantage over the nested
fixed-point fix persists there too.

## Practical takeaway (regardless of whether/when this gets fixed upstream)

When advising on Gamma GLMMs (or other estimated-dispersion families) in
`lme4`: warn that `glmer`'s default Laplace fit can substantially inflate
random-effect variance estimates, especially with more than one grouping
factor (where `nAGQ>1` isn't available as a cross-check, since it's
restricted to single-scalar-RE models). Cross-check with `glmmTMB` when the
precision of variance components matters. The bias is worst when the true
dispersion is small (shape large) — i.e. in the regime most real Gamma GLMM
applications actually live in.

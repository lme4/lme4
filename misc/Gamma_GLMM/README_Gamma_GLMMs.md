# Gamma GLMM random-effect variance bias (GH #643)

Investigation notes, 2026-07-28, comparing `lme4pureR`, `lme4` (dev HEAD,
v2.1-0), and `glmmTMB`. Originating issue:
<https://github.com/lme4/lme4/issues/643> ("Reliability of glmer to fit
Gamma-distributed model - Strange Ranef Variance"), most recently revived by
<https://github.com/TiagoAMarques/Report4BB>. Scratch scripts referenced below
live in the session scratchpad (not committed); `misc/Gamma_GLMM/GH643.R` in this repo
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

Two further diagnostics (reprex scripts `misc/Gamma_GLMM/reprex_pirls_fragility.R`
and `misc/Gamma_GLMM/reliability_comparison.R`, run against lme4 2.0-6 installed
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

**Reliability comparison, completed** (`misc/Gamma_GLMM/reliability_comparison_params.R`,
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

Follow-up to §9's finding above (`misc/Gamma_GLMM/multimodal_check.R`): is the
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
joint-phi arm (run afterward, alone). Scripts: `misc/Gamma_GLMM/multistart_lme4.R`,
`misc/Gamma_GLMM/multistart_glmmTMB.R`, `misc/Gamma_GLMM/multistart_jointphi.R`,
`misc/Gamma_GLMM/multistart_analysis.R`; plots: `misc/Gamma_GLMM/multistart_ecdf.png`,
`misc/Gamma_GLMM/multistart_dispersion_hist.png`.

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

Visually (`misc/Gamma_GLMM/multistart_ecdf.png`): lme4-old rises in one fairly
sharp jump near its own dominant plateau with a short tail. **lme4-current
also jumps early to its (better/lower) dominant plateau, but then has a
visibly longer, more gradual staircase tail out to ~1500** — i.e. the fix's
dominant mode is objectively a better fit than old lme4's dominant mode,
but the fix's optimization landscape has *more* reachable alternative local
optima than old lme4's, not fewer. That's a genuine, previously-unstated
cost of the fix worth being aware of, distinct from the boundary-collapse
question in §9-§10.

### Dispersion (phi = sigma²), true value 4

`misc/Gamma_GLMM/multistart_dispersion_hist.png`:

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
method (`misc/Gamma_GLMM/multistart_pairs_lme4old.png`,
`_lme4current.png`, `_glmmTMB.png`; script `misc/Gamma_GLMM/multistart_pairs.R`):

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
prototype (`misc/Gamma_GLMM/pirls_joint_phi.R`, §7) before the C++ fix existed:
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
construction. Script: `misc/Gamma_GLMM/multistart_jointphi.R`.

**Results: this comparison confirms the hypothesis.** All 200/200 starts
converged "clean" (no errors, warnings, or sentinel hits — contrast both
lme4 fixed-point variants above, which had 15-37 warnings and, for 2.0-6,
29 outright errors). The dominant plateau is at **-2logLik ≈ 1329
(158/200 = 79%)** — much closer to glmmTMB's dominant plateau (1326,
190/200) than either fixed-point lme4 variant's dominant plateau (2.0-6:
1393, 118/200; current fix: 1360, 131/200). Dispersion recovery
(`misc/Gamma_GLMM/multistart_dispersion_hist.png`) is essentially unbiased at
the mode — median phi = 4.014 against a true value of 4, matching
glmmTMB's 4.009 — versus both fixed-point variants' persistent downward
bias (2.0-6 median 3.51, current-fix median 3.32). The pairs plot
(`misc/Gamma_GLMM/multistart_pairs_jointphi.png`) shows a tight dominant cluster
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

## 12. Parameter-recovery survey: six methods, four real datasets

Broader successor to the multistart diagnostic (§11): rather than
starting-point sensitivity on one dataset, fit **six** methods to **B**
resimulated datasets from **four** real designs, each simulated from
"pretty" (rounded) parameters near a real reference fit, and compare
parameter recovery, reliability (clean/warning/error/singular), timing,
and -2*logLik. Toolkit and scripts: `misc/Gamma_GLMM/paramsurvey/`
(`toolkit.R` plus `01_prep_*.R`, `02_fit_glmmTMB.R`, `03_fit_jointphi.R`,
`04_fit_lme4current.R`, `05_fit_lme4old.R`, `06_fit_pirls_phi.R`,
`07_analysis.R`, `08_summary_plots.R`; see that directory's own README
for a per-script rundown). **Caveat**: these scripts still
hardcode this session's scratch-directory paths (`wd <-
"/tmp/claude-.../scratchpad/param_survey"` etc.) — update those before
rerunning outside this session.

### Six methods

1. **glmmTMB** — reference/gold standard, full joint TMB optimization.
2. **joint-phi** — R-level devfun (§7/§11's prototype, generalized), phi
   as a first-class outer `bobyqa` parameter alongside theta/beta.
3. **PIRLS/digamma** — R-level nested fixed-point (like the C++ fix, but
   phi re-estimated via the true Gamma conditional MLE each outer
   iteration, not the moment plug-in); resurrects §6's "ALSO FALSIFIED"
   idea to check whether that single-RE-scenario null result generalizes.
4. **PIRLS/moment** — same nested fixed-point architecture, moment
   estimator (`phi = dev/n`) — the *same algorithm* as the C++ fix, in R;
   comparing timing against method 5 isolates pure R-vs-C++ overhead.
5. **lme4 current** — the compiled C++ fix (`glmer()`, main env install).
6. **lme4 2.0-6** — unmodified, isolated install (§9's baseline).

### Four example datasets

- **epil2 (simple)**: `y ~ trt + (1|subject)`, Gamma(log), the reduced
  form from §1-§2 (58 subjects, 213 obs after dropping y=0). Non-singular
  without a prior.
- **epil2 (complex)**: the full model from `?glmmTMB::epil2` (`y ~
  Base*trt + Age + Visit + (Visit|subject)`) — 2-term correlated random
  slope. Needed a regularizing prior (`normal(0,3)` on ranef) for the
  *reference* fit only (to pin down non-degenerate true parameters);
  per-replicate fits use no prior and are frequently singular by design.
- **Report4BB**: the real dataset/model that started this whole
  investigation (§1), `crate ~ (1|location) + (1|fyear)`, n=103, 7/10
  levels — two independent single-term REs, no correlation. Reconstructed
  from a fresh clone of github.com/TiagoAMarques/Report4BB (data-prep
  code in `01_prep_report4bb.R` mirrors `testing_lme4/Testing_lme4.R`
  there). Frequently singular (few levels per factor) even at
  non-degenerate true parameters — expected, matches the well-known
  small-number-of-groups instability noted since project memory.
- **schizophrenia**: bundled directly in lme4 (`data(schizophrenia)`,
  `man/schizophrenia.Rd`'s own example model, `imps79 ~ TxDrug*Week +
  (1|id)`), 1603 obs, 437 subjects — much larger/higher-powered than the
  other three, single random intercept, non-singular without a prior.

### Finding: a second, independent bug in `sigma()`'s dispersion reporting

While comparing PIRLS/moment (R) against lme4-current (C++) — supposedly
the *same* algorithm — on identical simulated data, their reported phi
values diverged by up to ~2.5x even though parameter estimates (sd, beta)
looked similar. Diagnosed by evaluating the R-level devfun at the exact
`(theta,beta)` lme4-current converged to: `-2*logLik` matched to 5
significant figures (confirming the fit/curvature itself is correct) but
`sigma(fit)^2` did not match the internally-converged phi at all.

Root cause: `sigma(fit)^2` for GLMMs is computed via a generic,
family-agnostic formula, `sigmaML <- (wrss + ||u||^2)/n` (`R/utilities.R`),
where `wrss = resp$wrss()` is the sum of squared *working* residuals —
computed using the same `1/phi`-reweighted working weights the Gamma_GLMM
fix (§9) introduced. Once phi moves away from 1, this formula becomes
self-referentially confounded with the very quantity it's meant to
report; it was only ever correct because old lme4 held phi fixed at 1.

Fixed (commit `1b7bfcd8`, pushed to `Gamma_GLMM`): exposed the internally-
converged phi to R directly (`glm_phi` in `src/external.cpp`, mirroring
the existing `glm_theta` accessor; `phi()` method on the `glmResp`
RefClass in `R/AllClass.R`) and use it directly in `R/utilities.R`'s
`sigmaML` computation whenever `isGLMM && !hasNoScale(resp$family)`,
instead of the generic `pwrss/n` path.

**At the same time, generalized the whole phi-profiling fix beyond
Gamma**, per explicit request: it was gated on `family()=="Gamma"` only,
but the underlying math (`resDev()`, `Laplace()`/`aic()`) already
dispatches generically per-family via `glmFamily`, so nothing else needed
to change. Added `glmDist::hasFreeDispersion()` (`src/glmFamily.h`),
overridden per subclass — `true` for Gamma/gaussian/inverse.gaussian,
`false` for poisson/binomial/negative.binomial, mirroring `hasNoScale()`
in `R/utilities.R` — delegated through `glmFamily`/`glmResp`, and used to
gate `glmerLaplace()`'s fixed-point loop instead of the old hardcoded
string check. This is a real, if rare, expansion of scope: gaussian
fit via `glmer()` with a non-identity link (e.g. `gaussian(link="log")`)
now also gets phi-profiled, which it never did before. The base
`glmDist::hasFreeDispersion()` **throws an informative error** for
unrecognized/user-supplied custom families rather than silently guessing
— **TODO, flagged for later**: reconsider this (a warning + a documented
default assumption, rather than a hard error, was suggested).

Confirmed a real, independent finding once the reporting bug was fixed:
even with `sigma()` reporting the true converged phi correctly, the
**moment/digamma estimators themselves carry a genuine ~15-25% bias**
in a shape≈0.5-1 regime, previously masked by the reporting bug
coincidentally landing closer to the true value by chance for some seeds.
Confirmed via `test-glmFamily.R`'s existing shape-recovery tests (2
tolerances loosened, with comments explaining why — not a regression).
Not yet root-caused.

Full GLMM-relevant test suite (`LME4_TEST_LEVEL=2`): 443 pass, 0 fail
after refreshing 8 stale hardcoded references across `test-glmFamily.R`,
`test-glmer.R`, `test-methods.R` (all either mechanical staleness from
the above, or newly-in-scope non-canonical-link-gaussian fits).

### Finding: a third, independent bug — "+2" AIC bookkeeping leaking into logLik/AIC/BIC

Found while checking why PIRLS/moment (C++) and PIRLS/moment (R) — again
supposedly the same algorithm — showed -2\*logLik values a suspiciously
constant 2 units apart on Report4BB, more visible there than elsewhere
because the differences are otherwise small. `glmResp::Laplace()`
(`src/respModule.cpp`) computes the Laplace deviance as
`ldL2 + sqrL + aic()`, where `aic()` dispatches to the family's C++
`aic()` method (`src/glmFamily.cpp`), a direct port of base R's
`family$aic()` (e.g. `Gamma()$aic`, `gaussian()$aic`,
`inverse.gaussian()$aic`). Those base-R functions bake in an additive
"+2" bookkeeping constant meant to be cancelled by the "p - aic/2" trick
in `stats::logLik.glm()` (`p` includes the dispersion parameter when
estimated) — a convention `glmResp::Laplace()` never applied, so it leaked
straight into `logLik()`/`AIC()`/`BIC()` for any GLMM with a freely-
estimated dispersion parameter (Gamma, inverse Gaussian, gaussian fit via
`glmer()`). Confirmed present on `master` (predates this whole
investigation by over a decade — `git blame` traces it to 2011), unrelated
to the RE-variance-bias or `sigma()` bugs above.

Fixed (commit `db6d172c`, pushed to `Gamma_GLMM`): subtract 2 in
`Laplace()` when `hasFreeDispersion()` is true. Verified against an
independent recomputation of `-2*logLik` straight from `dgamma()`/`dnorm()`
(bypassing `aic()` entirely, using the fitted model's own converged
conditional modes) — locked in as a regression test in
`test-gamma_glmm_bias.R`. One existing hardcoded test reference
(`test-glmer.R`, a `gaussian(link="log")` fit) needed refreshing; a full
audit of the test suite found no other hardcoded logLik/AIC/deviance
references for a free-dispersion-family GLMM. NEWS entry added under
2.1-0 BUG FIXES.

### Status as of this writing

B=500 completed for all 4 datasets × 6 methods (24 cells, up from the
original B=10); results, combined `.rds` files, and four summary plots
in `misc/Gamma_GLMM/paramsurvey/` (`param_summary_distrib.png`,
`param_summary_stderr.png`, `param_summary_stderr_newonly.png` —
per-parameter recovery, full distribution and mean±2SE variants;
`time-negll_summary.png` — elapsed time and paired per-replicate
Δ(-2\*logLik) vs glmmTMB, two panels via patchwork) committed; see
`paramsurvey/README.md` for a full script-by-script rundown. The
`-2*logLik` values were post-processed after the "+2" bug fix above (dev-
build/2.0-6 methods only; the R-level methods never went through `aic()`
and were already correct).

Headline pattern across all four datasets: lme4-2.0-6 shows large,
dataset-dependent bias (inflation in most regimes, e.g. epil2-complex sd2
~14x true, Report4BB sd1/sd2 ~2x); all five "modern" methods (glmmTMB,
joint-phi, PIRLS/digamma, PIRLS/moment, lme4-current) track each other
closely on parameter recovery, with PIRLS/moment and lme4-current
essentially identical (as expected, same algorithm) and consistently
~15-25% low on phi relative to glmmTMB/joint-phi (the estimator bias
found above). Timing is highly structure-dependent, not sample-size-
dependent: schizophrenia (n=1603, single RE) is *cheap* for every method
(glmmTMB and joint-phi both ~2.5s), while epil2-complex (n=213, one
correlated 2-term RE) is the expensive outlier (PIRLS/digamma,
PIRLS/moment ~100s/fit) — the correlated-RE structure drives cost, not n.

**Not yet done**: root-cause the moment/digamma ~15-25% phi bias; decide
on the custom-family error-vs-warning question above.

## Practical takeaway (regardless of whether/when this gets fixed upstream)

When advising on Gamma GLMMs (or other estimated-dispersion families) in
`lme4`: warn that `glmer`'s default Laplace fit can substantially inflate
random-effect variance estimates, especially with more than one grouping
factor (where `nAGQ>1` isn't available as a cross-check, since it's
restricted to single-scalar-RE models). Cross-check with `glmmTMB` when the
precision of variance components matters. The bias is worst when the true
dispersion is small (shape large) — i.e. in the regime most real Gamma GLMM
applications actually live in.

## Backward-compatibility switch, and `nAGQ>1` (AGQ) scoping (2026-08-05)

### `glmerControl(disp_method=)`

Added `glmerControl(disp_method = c("moment", "old/buggy"))` and
`glmerControl(maxPhiIter = 100L)`. `"moment"` (default) is the existing
nested-fixed-point fix described above; `"old/buggy"` reproduces the
original, pre-fix behaviour exactly (phi held fixed at 1 throughout PIRLS,
matching old CRAN `glmer` output bit-for-bit on a spot-checked example) for
users who need it for back-compatibility. `maxPhiIter` caps the nested
loop (was previously hardcoded to 100 in `external.cpp`).

Implementation notes:
- One PIRLS implementation (`pwrssUpdate()`) is shared by both code paths;
  `disp_method="old/buggy"` just takes the same single-call branch
  fixed-dispersion families (binomial/Poisson) already used, rather than
  needing a second PIRLS chunk.
- `sigma()`/`sigmaML` (`R/utilities.R`, `mkMerMod`) had to be made
  conditional on `disp_method`: it reads `resp$phi()` directly when
  `dispProfile=TRUE` (profiled value), but falls back to `pwrss/n` when
  `dispProfile=FALSE` -- under `"old/buggy"`, `resp$phi()` is stuck at 1 and
  meaningless, but `pwrss/n` is *not* confounded in that case (working
  weights were never reweighted by 1/phi either), so it correctly
  reproduces old CRAN's reported dispersion.
- `disp_method`/`maxPhiIter` are stored in `devcomp$dims` so `refit()` can
  recover the original fit's setting rather than silently reverting to
  package defaults.
- Build gotcha hit during development: after editing `external.cpp` and
  running a normal (incremental) `R CMD INSTALL`, only the touched
  translation units got recompiled, reusing stale `.o` files for the rest
  -- this produced a real segfault (null vtable dispatch through
  `hasFreeDispersion()`) from an ABI mismatch between the fresh and stale
  object files, for *every* `glmer()` call, not just Gamma ones. Fixed by
  a full clean rebuild (`rm src/*.o src/*.so && R CMD INSTALL --preclean`).
  Worth remembering if a similarly inexplicable native crash shows up again
  after a C++-only edit.

### `nAGQ>1` (AGQ) fixed (2026-08-06)

Originally, `glmerAGQ` (`external.cpp`, the `nAGQ>1` code path) called
`pwrssUpdate()` directly, bypassing `glmerLaplace()` entirely, and never
touched `rp->setPhi()`/`disp_method`/`maxPhiIter`. What was actually
happening (verified empirically, not just "phi==1 always"): phi got
profiled once during the `nAGQ=0` init stage (which does go through
`glmerLaplace`), then stayed *frozen* at that stale value for the entire
`nAGQ>1` optimization -- it didn't track theta as the AGQ-based
Nelder-Mead search moved it. This turned out to be a two-layer problem:

1. **Mode-finding/weights**: `pwrssUpdate()` and the nested-phi-loop
   pattern from `glmerLaplace` needed to be reusable for the
   u-mode-finding step inside `glmerAGQ` (which already calls
   `pwrssUpdate()` once before its quadrature loop).

2. **The returned marginal-likelihood formula itself**: `glmerAGQ`'s
   final value was `devc0.sum() + ldL2 - 2*log(mult.prod())`, a
   McCullagh & Nelder deviance-based AGQ approximation -- `devc0.sum()`
   is the raw, phi-independent unit deviance, never passed through the
   family's actual density. This is fine for binomial/Poisson
   (deviance-based and density-based -2logLik coincide there) but wrong
   for Gamma: deviance is `2*phi*(l_sat - l_fit)`, missing the
   phi-dependent normalizing terms (`log(mu*phi)/phi`, `lgamma(1/phi)`)
   that a real `dgamma()`-based likelihood carries.

**The fix**: both layers turned out to reduce to reuse, not a from-scratch
re-derivation. `external.cpp` gained a shared `profilePhi()` helper
(factored out of `glmerLaplace()`'s existing nested fixed-point loop),
called by both `glmerLaplace()` and `glmerAGQ()` before mode-finding.
`devcCol()` (the per-grouping-level squared-mode-plus-deviance helper used
by the AGQ quadrature loop) now divides its deviance-residual sum by the
profiled `phi`, same split as everywhere else in this fix (the `u^2`
penalty stays unscaled). And the key insight for layer 2: the missing
phi-dependent normalizing constant is *constant with respect to u* (a
function of `y` and `phi` only, standard exponential-dispersion-family
decomposition), so it cancels automatically inside the AGQ ratio and only
needs to be added once, to the final returned value -- which is exactly
what `rp->Laplace(ldL2, 0., sqrL)` already computes (reusing the same
`aic()`-based machinery `glmerLaplace()` uses, self-consistent with the
profiled phi at the converged mode). So `glmerAGQ()`'s final line is now
`rp->Laplace(pp->ldL2(), 0., pp->sqrL(1.)) - 2*log(mult.prod())` instead of
`devc0.sum() + pp->ldL2() - 2*log(mult.prod())`.

This also fixes the same constant-offset bug for binomial/Poisson
`nAGQ>1` fits, which was never Gamma-specific -- `devc0.sum()` omitted the
same y-only density-normalizing constant for every family, just harmlessly
so, since it doesn't affect where the optimizer converges (a fact that
was already, if implicitly, acknowledged: `anova()` had a hardcoded guard
refusing to compare an `nAGQ>1` glmer fit against a `glm()` object,
`stop("...incommensurate with glm() objects")`, `R/lmer.R`, now removed
since it's no longer true).

**Validated**: on a single-scalar-RE Gamma fit, `nAGQ=1` (Laplace) and
`nAGQ=5` (AGQ) now agree closely on sigma/logLik/VarCorr, where they
previously diverged sharply (frozen/stale phi). `nAGQ=2` agrees even more
closely with `nAGQ=1` than `nAGQ=5` does (~10x tighter on logLik),
consistent with AGQ genuinely converging toward the marginal likelihood as
quadrature order increases, rather than the agreement being a coincidence
at one particular order. The nAGQ>1 warning is removed
(`updateGlmerDevfun()`, `R/modular.R`); the now-obsolete
`test-gamma_glmm_bias.R` warning-existence test was rewritten into a real
correctness check.

One side effect surfaced by full-suite validation (unrelated to the AGQ
fix itself, but only visible once `glmerAGQ()` started calling
`hasFreeDispersion()`): fitting a GLMM with a deliberately malformed/
custom `family` object now throws `glmer()`'s existing "not implemented
for this family" hard error immediately at fit time for `nAGQ>1`, not just
for free-dispersion families under the Laplace path -- resolved below
(hard error vs. warning-with-default-assumption question).

## `hasFreeDispersion()` resolved generically via `$dispersion`, not a hard error (2026-08-06)

The open hard-error-vs-warning question above is resolved, and not by
picking one of those two options: unrecognized-by-name families are now
resolved generically via the family's own `$dispersion` component
(`?stats::family`, R >= 4.3.0's "Value" section: `NA_real_` = free
dispersion, a fixed numeric value = fixed dispersion) instead of
unconditionally erroring. Only genuinely incomplete custom families
(no `$dispersion` component either) still hard-error, now with a message
that includes the offending family's name and points at `?stats::family`.
`src/glmFamily.h`/`.cpp`: the base (unrecognized-name) `glmDist` class
gains `d_hasDispersionField`/`d_dispersionField`, read from the family
list's `dispersion` element in the constructor; `hasFreeDispersion()`
uses `ISNA()` on that value instead of always throwing. The six
name-dispatched subclasses (binomial/Gamma/gaussian/inverse.gaussian/
negative.binomial/poisson) are untouched -- they never reach the base
class -- so `MASS::negative.binomial()` (whose `$dispersion` field is
`NULL`, an unrelated gap in that family constructor) is unaffected.

Nicely, this also fixes `test-predict.R`'s "junk"-family test (a
mangled-name `inverse.gaussian()` object, `ig$family <- "junk"`) with
**no test-file edit needed**: `ig$dispersion` (inherited, untouched by
the name mangling) is still `NA_real_`, so `glmer()` now succeeds
(dispersion is unambiguously free) and only `simulate()` fails, on name
dispatch, with its own original, more specific message -- exactly the
pre-AGQ-fix behaviour the test always expected.

## `refit()`/`bootMer()` crash on pre-branch saved fit objects (2026-08-06)

Found via `confint(fit_cbpp_2, method="boot")` failing with "*all*
bootstrap runs failed" (`fit_cbpp_2` is a pre-baked fit loaded from
`inst/testdata/lme-tst-fits.rda`) -- initially assumed unrelated/
pre-existing since it reproduced on a commit that looked "clean," but
that commit (`5099a385`) postdates the `disp_method`/`maxPhiIter` dims
addition (`a46c21e3`), so it wasn't actually clean of this branch's
changes. `bootMer()` swallows the real per-replicate error; a direct
`refit()` call surfaced it: `Error in dc$dims[["dispProfile"]] :
subscript out of bounds`.

`refit.merMod` (`R/lmer.R`, ~line 1520) tries to recover the original
fit's `disp_method`/`maxPhiIter` settings from `devcomp$dims`, falling
back to package defaults for fits from before those fields existed --
intentional, per the code comment already there. But `dc$dims[["..."]]`
(double-bracket) throws for a genuinely-absent name in a named atomic
vector rather than returning `NULL`, so the fallback (`%||% TRUE`) never
gets a chance to run; base R's `%||%` also only rescues `NULL`, not
`NA`, so a single-bracket swap alone wouldn't have been enough either.
Fixed with single-bracket lookup (`NA` instead of an error) plus an
explicit `is.na()` fallback check, for both `dispProfile` and
`maxPhiIter`. Pure R change; verified `refit()` and
`confint(method="boot")` both now succeed on `fit_cbpp_2`.

## GH #936: `gaussian(link="log")` sigma gap -- confirmed real and systematic (2026-08-06)

[GH #936](https://github.com/lme4/lme4/issues/936) ("matching up
gaussian(link = "log") results across packages") reports that on
`nlme::Rail` (n=18, 6 groups), `glmer` and `glmmTMB` agree on fixed
effects but disagree sharply on `sigma` (3.33 vs 4.0) and logLik (-73 vs
-64.5); `mgcv::gam` and `MASS::glmmPQL` agree with `glmmTMB`. This
branch's dispersion fix already closed most of the logLik gap on the
real data (-72.3 to -65.1, per an earlier comment on the issue), but left
`sigma` still visibly off, with the open question of whether that's a
real remaining bug or just small-sample noise on one n=18 dataset.

Resolved that question with a mini parameter-recovery survey
(`misc/Gamma_GLMM/paramsurvey_loggaussian/`, modeled on `paramsurvey/`
but much more limited: one dataset, one family, four methods -- glmmTMB,
glmer, mgcv, glmmPQL, no R-level joint-phi/PIRLS-phi devfun arms).
Reference parameters rounded from a real glmmTMB fit to the real Rail
data (beta=4.1, RE sd=0.4, sigma=4.0), B=100 datasets resimulated from
those parameters via `glmmTMB::simulate_new()` using Rail's real design,
all four methods refit to each replicate (all 100 converged cleanly for
every method, no singular fits, no warnings):

| method | beta (4.10) | sd (0.40) | **sigma (4.00)** | negll |
|---|---|---|---|---|
| glmmTMB | 4.099 | 0.321 | **3.76** | 125.28 |
| **glmer** | 4.100 | 0.320 | **3.08** | 126.37 |
| mgcv | 4.101 | 0.321 | **3.76** | 125.28 |
| glmmPQL | 4.101 | 0.319 | **3.76** | NA (not a real ML method) |

**Confirmed real and systematic, not small-sample noise.** `beta` and RE
`sd` are essentially unbiased and agree tightly across all four methods.
glmmTMB/mgcv/glmmPQL agree with each other almost to 3 decimal places on
`sigma` (a modest, expected ~6% downward ML bias with only 6 groups) --
but `glmer`'s median `sigma` is 3.08, roughly 18% further off than the
other three, about 23% below the true value. Since the logLik gap is
already fixed (this branch's core contribution) and the remaining `sigma`
gap is glmer-specific (not shared by glmmTMB, which also uses a Laplace
approximation), this points at something in how `glmer` estimates or
reports `sigma` specifically for `gaussian()` fit with a non-identity
link, not the phi-profiling/AGQ mechanism this branch already fixed.

### Mechanism, precisely characterized: a missing degrees-of-freedom correction (2026-08-06)

Prompted by the question of whether the gap is an `n` vs. `n-p` (number
of fixed-effect parameters) denominator issue -- the intercept-only Rail
model has `n=18`, `p=1`, so this predicts too small a correction to
matter (`RSS/(n-p) = RSS/17` vs. `RSS/n = RSS/18`, barely different).
Checked directly on the real Rail data: glmer's reported `phi`
(`sigma^2`) *is* essentially exactly `RSS/n` at the conditional mode
(10.829 vs. reported 10.828) -- confirming glmer's moment-based
phi-profiling loop already uses the plain `n`, with no `p`-correction of
any kind, so an `n` vs. `n-p` explanation is ruled out (and even if
applied, points in far too small a direction to explain an 18%+ gap).

But the *right* correction, empirically: not `n-p` (fixed-effect
parameter count) but `n-q`, where `q` is the number of random-effect
**levels/groups** (here, `q=6` Rail groups). glmmTMB's `phi` on the real
data (16.146) is almost exactly `RSS/(n-q) = RSS/12 = 16.33`, vs.
glmer's `RSS/n = RSS/18 = 10.83`. Checked this isn't a one-dataset
coincidence: across all B=100 replicates in the sweep above, the
per-replicate ratio `phi_glmmTMB / phi_glmer` is **1.487 (SD 0.0097)** --
essentially constant -- matching the predicted `n/(n-q) = 18/12 = 1.5`
to within ~1%. This is a strikingly tight match for a quantity derived
from a completely different estimation procedure (glmmTMB's joint
marginal MLE via Laplace-approximated integration over `u`, vs. glmer's
PIRLS-conditional-mode moment plug-in), and precise enough to plausibly
be *fixable*, not just diagnosable.

**Interpretation**: glmer's phi-profiling loop estimates phi as
`deviance/n` at the current PIRLS conditional mode of `u` -- i.e., it
treats the `q` fitted random-effect levels as if they were fixed/known,
the same way an ordinary GLM's `deviance/n` treats fixed effects. It
never discounts for the degrees of freedom effectively consumed by
*estimating* those `q` levels (with shrinkage) rather than knowing them
outright. This is conceptually the same correction that separates REML
from ML variance estimation in classical (L)MMs -- except this isn't a
REML-vs-ML distinction (glmmTMB here is doing ML, not REML): it's that
glmer's moment/plug-in shortcut for phi applies *no* correction at all
(behaves as though `q=0`), while a proper joint marginal likelihood
naturally bakes it in. This is very likely the same underlying mechanism
as the still-unexplained "moment/digamma estimators carry a genuine
~15-25% bias" finding for Gamma (README, "Finding: a second, independent
bug in `sigma()`'s dispersion reporting") -- not previously connected to
a concrete, quantifiable degrees-of-freedom formula until this check.

**Caveat, scope not yet checked**: `q` was unambiguous here (one
grouping factor, no crossing). For crossed or multiple random-effect
terms, what the "right" `q` is is not obvious -- naively summing levels
across terms is one guess, but crossed designs share information across
factors in a way a simple level-count may not capture correctly, and
glmmTMB handles crossed REs natively while `mgcv`/`glmmPQL` (this
sweep's other two commensurate-`sigma` methods) have much more limited
or no support for them -- weakening the cross-method validation this
finding relied on. **Not yet extended beyond the single-grouping-factor
case**.

### Confirmed via a fifth method, joint-phi: same mechanism, isolated more cleanly (2026-08-06)

The `n-q` diagnosis above was cross-implementation (glmer's C++/PIRLS
stack vs. glmmTMB's independent TMB/AD stack) -- consistent with a real
effect, but not a controlled comparison, since the two implementations
differ in far more than just how phi is handled. Added a fifth method to
`misc/Gamma_GLMM/paramsurvey_loggaussian/` to isolate the one variable
that matters: **joint-phi**, an R-level devfun (adapted from
`../paramsurvey/toolkit.R`'s Gamma version, §7/§11) that shares glmer's
exact PIRLS mode-finding code and starting values, but promotes `phi` to
a first-class parameter optimized jointly with theta/beta via `bobyqa`,
rather than profiling it via glmer's nested `deviance/n` moment plug-in.
Prediction: if the moment plug-in is really the culprit, joint-phi's
`sigma` should match glmmTMB/mgcv/glmmPQL, not glmer, despite sharing
glmer's own mode-finding machinery.

**Confirmed, precisely, on the same B=100 replicates as the four-method
sweep above** (all 100 converged cleanly, no singular fits, ~0.53s/fit):

| method | sigma (4.00) | ratio to glmmTMB |
|---|---|---|
| glmmTMB | 3.758 | 1 (reference) |
| **joint-phi** | 3.759 | **1.0013 (SD 0.0006)** |
| mgcv | 3.757 | 0.9997 |
| glmmPQL | 3.764 | 1.0016 |
| glmer | 3.084 | 0.821 |

joint-phi lands essentially exactly on glmmTMB/mgcv/glmmPQL. More
tellingly, the per-replicate ratio `phi_glmer / phi_jointphi` is
**0.6716 (SD 0.0044)**, matching the predicted `(n-q)/n = 12/18 =
0.6667` even more tightly (~0.7% off) than the glmer-vs-glmmTMB
comparison was (~1% off) -- because isolating the *one* difference
between glmer and joint-phi (moment-plug-in phi vs. jointly-optimized
phi, holding the mode-finding code itself fixed) removes the
"different software stack" confound that the original glmmTMB
comparison couldn't rule out. This is about as clean a confirmation of
the mechanism as this kind of simulation study can give.

Scripts: `misc/Gamma_GLMM/paramsurvey_loggaussian/toolkit.R`
(`make_joint_phi_devfun_gaussian`/`fit_jointphi_one`),
`03_fit_jointphi.R`. The crossed-RE caveat above still applies --
joint-phi confirms the mechanism for the single-grouping-factor case
specifically, not the general one.

### Crossed random effects: `q_eff = sum(q_k) - (K-1)`, confirmed empirically (2026-08-06)

The single-grouping-factor finding above left an open question: what
should `q` mean with multiple/crossed grouping factors, since naively
summing levels (`q1+q2`) was only a guess -- flagged explicitly as
untested in the TODO.

Tested directly with a new extension,
`misc/Gamma_GLMM/paramsurvey_loggaussian/08_prep_crossed.R` through
`12_analysis_crossed.R`: a 16x16-level fully-crossed factorial base
design (`Rail1 x Rail2`, 3 reps, n=768; RE sd=1 on the log scale for each
factor, residual sigma=0.01) simulated once via
`glmmTMB::simulate_new()`, then 8 fixed subsamples fit with
glmmTMB/glmer/mgcv (B=10 replicates each, all sharing the same
underlying simulated realizations since every case subsamples the same
master simulation):

- 3 matched pairs of "structured" (k x k x 3, full crossing, k=4/5/6) vs
  "random" (uniform row subsample at the same n) crossed designs, all
  fit with `~1+(1|Rail1)+(1|Rail2)`
- 2 "oneway" control cases (all levels of Rail1, only 1 level of Rail2,
  at two different Rail1 level counts) -- effectively a plain one-factor
  design, fit with `~1+(1|Rail1)`, where `q` is unambiguous

Because each case's design is fixed across replicates (only the response
is resimulated), the `glmer/glmmTMB` phi ratio (`phi = sigma^2`) turns
out to be **deterministic given the design, not just unbiased in
expectation** -- its SD across B=10 replicates is 4-5 orders of
magnitude smaller than the effect itself. That determinism meant B=10
was already enough to resolve the implied `q_eff = n*(1 -
phi_glmer/phi_glmmTMB)` to within ~0.03 of an integer on every design
tried, an early (RE sd=3) attempt at this design broke glmer's PIRLS
Cholesky step outright on every replicate -- see `08_prep_crossed.R`'s
header comment; RE sd=1 avoids it).

| case | n | q (= sum of levels) | q_eff | q1+q2-1 or q1 |
|---|---|---|---|---|
| structured (4x4x3) | 48 | 8 | 6.988 | 7 |
| random (n=48) | 48 | 30 | 28.994 | 29 |
| oneway (q1=16) | 48 | 16 | 15.991 | 16 |
| structured5 (5x5x3) | 75 | 10 | 8.984 | 9 |
| random75 | 75 | 32 | 30.988 | 31 |
| structured6 (6x6x3) | 108 | 12 | 10.976 | 11 |
| random108 | 108 | 32 | 30.980 | 31 |
| oneway10 (q1=10) | 30 | 10 | 9.995 | 10 |

**Both one-way controls match `q_eff = q1` almost exactly; all six
crossed designs match `q_eff = q1+q2-1`, not `q1+q2`** -- consistently
off by very close to one full level, regardless of design size, balance,
or crossing degree (structured vs. random subsample).

**Clean explanation**: classical ANOVA identifiability. With a shared
fixed intercept plus two crossed grouping factors, each factor's own set
of level-effects has exactly one redundant direction relative to the
shared intercept (the usual sum-to-zero constraint), so the count of
independent mean-structure parameters is `1 + (q1-1) + (q2-1) =
q1+q2-1`, not `1+q1+q2`. For a single grouping factor this collapses to
`1+(q1-1) = q1` -- exactly what the one-way controls show. Generalizes
naturally to `K` crossed intercept-only grouping factors: `q_eff =
sum(q_k) - (K-1)`.

**Not yet checked**: nested (rather than crossed) grouping factors,
random slopes (not just intercepts), and correlated multi-term random
effects -- all likely need a different or more general
(rank-of-design-matrix-based) treatment. This result is specifically for
additive, intercept-only, crossed grouping factors.

Scripts: `08_prep_crossed.R` (simulates the shared master dataset,
derives all 8 case subsamples), `09_fit_glmmTMB_crossed.R` /
`10_fit_glmer_crossed.R` / `11_fit_mgcv_crossed.R` (dispatch to crossed
vs. one-way fit wrappers by case name -- glmmPQL and joint-phi are
skipped here, see `toolkit.R`'s comments), `12_analysis_crossed.R`
(combines and summarizes).

### `q_eff` finding confirmed family-general: Gamma matches gaussian (2026-08-06)

The `q_eff = sum(q_k) - (K-1)` result above was found under
`gaussian(link="log")`; there's no obvious reason it should be
gaussian-specific (it's an argument about the mean-structure design
matrix's rank, independent of the response distribution), but that's
worth actually checking rather than assuming. `08_prep_crossed.R`
through `12_analysis_crossed.R` now take a `family` argument (CLI arg 2
to `08_prep_crossed.R`, "gaussian" or "Gamma"; propagated through
`toolkit.R`'s crossed/one-way fit wrappers and saved/loaded via a
`crossed_simdata_<family>.rds` / `crossed_results_..._<family>.rds`
naming convention so the two families' results don't clobber each
other) rather than hardcoding `gaussian(link="log")` throughout.

Two snags surfaced getting `Gamma(link="log")` working, both fixed in
the toolkit rather than by rewriting the pipeline:

- glmmTMB's current RTMB backend doesn't implement `Gamma` yet
  (`distribution not implemented yet for use with RTMB backend: Gamma`)
  -- fixed by calling `glmmTMB:::useRTMB(FALSE)` once at the top of
  `08_prep_crossed.R` and `toolkit.R`, forcing the legacy TMB backend for
  *both* families so the comparison stays apples-to-apples.
- glmmTMB's `betadisp` parameter means different things for different
  families: for gaussian, `exp(betadisp)` is sigma itself (matching
  `sigma(fit)`, already used by the single-factor Rail survey above);
  for Gamma, `exp(betadisp)` is the **shape** parameter (`1/sigma^2`,
  where `sigma(fit)` is the CV) -- confirmed empirically by fitting a
  known-shape simulated Gamma sample and checking `exp(betadisp)`
  against the true shape. Missing this the first time produced a
  wildly overdispersed (near-zero/huge-outlier) simulated response;
  `08_prep_crossed.R` now dispatches the true-parameter-to-`betadisp`
  mapping on `family_name`.

With those fixed, re-running the identical 8-case sweep (B=10) under
`Gamma(link="log")` gives:

| case | n | q | q_eff | q1+q2-1 or q1 |
|---|---|---|---|---|
| structured (4x4x3) | 48 | 8 | 6.986 | 7 |
| random (n=48) | 48 | 30 | 28.986 | 29 |
| oneway (q1=16) | 48 | 16 | 15.984 | 16 |
| structured5 (5x5x3) | 75 | 10 | 8.978 | 9 |
| random75 | 75 | 32 | 30.976 | 31 |
| structured6 (6x6x3) | 108 | 12 | 10.970 | 11 |
| random108 | 108 | 32 | 30.968 | 31 |
| oneway10 (q1=10) | 30 | 10 | 9.989 | 10 |

Essentially identical to the gaussian table (same designs, same
precision) -- the one-way controls match `q_eff = q1`, the crossed
designs match `q_eff = q1+q2-1`, confirming the mechanism is not
gaussian-specific. (Also notably: glmer had zero fit failures across
all 240 Gamma fits, cleaner than the 2/240 stray `Downdated VtV`
failures under gaussian on the same designs -- not investigated
further, but consistent with Gamma's canonical-ish log link being a
numerically easier case for PIRLS than gaussian's non-canonical one.)

### First implementation: `glmerControl(disp_dof_correction=TRUE)` (2026-08-06)

First pass at actually applying the `q_eff = sum(q_k) - (K-1)` finding
inside glmer, rather than just diagnosing it externally. New opt-in
`glmerControl()` argument, default `FALSE` -- deliberately not on by
default, since validation so far is simulation-only (one dataset
family, gaussian/Gamma, crossed scalar-intercept designs) with no real-
data check and no check of how it interacts with REML or `nAGQ>1`.

**Where `n`, `sum(q_k)`, and `K` live**: investigated before writing any
code (this was the actual ask that started this round). `n` was already
available in `profilePhi()` (`rp->weights().sum()`, `external.cpp`).
`sum(q_k)` (= total random-effects dimension) is also already computed
in C++ -- `merPredD`'s constructor sets `d_q = d_Zt.rows()`
(`predModule.cpp`), reachable from `external.cpp` via the existing
public `pp->Zt().rows()` with no header changes at all. `K` (the number
of independent random-effects *terms*, as opposed to the number of
covariance parameters) is **not** reliably recoverable from `Lambdat`/
`Lind` alone: for scalar-intercept terms you could count distinct
`Lind` values, but a single correlated random-slope term (`(1+x|g)`)
contributes 3 theta values (var, cov, var) for what is structurally one
term, so that approach overcounts `K` in general. The robust source of
truth is `reTrms$Gp` (R-level, already built by `mkReTrms()`):
`diff(Gp)` gives each term's `q_k`, `length(Gp)-1` gives `K` -- and this
was **not** previously passed into the C++ side at all (`merPredD$new()`
only ever received `Zt`/`theta`/`Lambdat`/`Lind`). One simplifying fact
found along the way: `nAGQ>1` already hard-restricts to a single scalar
random-effects term (`updateGlmerDevfun()`, `R/modular.R`, pre-existing
check), so the crossed-term correction only actually matters for the
Laplace (`nAGQ<=1`) path -- `glmerAGQ()` will only ever see `K=1`, where
the formula trivially gives `q_eff=q1`.

**Wiring** (mirrors how `disp_method`/`maxPhiIter` were already threaded
through as extra params for the original phi-profiling fix): `qEff` is
computed once in R, in `mkGlmerDevfun()` where `reTrms` is fully in
scope, via a new `computeQEff(reTrms, enable)` helper (`R/modular.R`);
stored on `rho$qEff`; threaded through `mkdevfun()`'s devfun closure and
`glmerPwrssUpdate()` exactly like `dispProfile`/`maxPhiIter`; passed as
a new trailing arg to `.Call(glmerLaplace, ...)` / `.Call(glmerAGQ,
...)` (both `CALLDEF` counts bumped); used in `profilePhi()` as `n -
qEff` instead of `n` when `R_finite(qEff)` (an `NA_real_` sentinel means
"no correction", the default). `glmerLaplaceHandle()` (the exported
handle for calling `glmerLaplace` directly) got the same new `qEff`
argument, defaulting to `NA_real_` for backward compatibility.

**Scalar-RE-only gate**: `computeQEff()` only computes a real correction
when `all(lengths(reTrms$cnms) == 1L)` -- every random-effects term is a
scalar random intercept, the only case the formula has been validated
for (see the crossed-RE section above). Otherwise it returns `NA_real_`
(no correction applied, matching pre-existing behaviour) and emits a
`warning()` explaining that sigma/dispersion estimates may be more
biased than necessary as a result. Deliberately conservative: rather
than guessing whether `q_k = nlevels_k * ncol_k` generalizes correctly
to random slopes or correlated blocks, it just doesn't apply anything
there yet -- see the TODO for the follow-up simulation study
(`ar1()`/random-slope designs) intended to test that generalization
before trying it in code. The working guess (untested) is that it
*will* generalize cleanly, because `q_k` here means the dimension of
the *spherical* random-effects vector `u` (always independent by
construction, `u ~ N(0, I_q)`, regardless of how `Lambdat` correlates
`Lambdat %*% u` afterward) rather than anything about the correlation
structure itself.

**Spot-checked** (not a full validation) on `sleepstudy` (Gamma,
`(1|Subject)`, n=180, q=18): `disp_dof_correction=TRUE` raises sigma by
a factor of 1.0547, close to the predicted `sqrt(n/(n-q)) = 1.0541`.
`disp_dof_correction=FALSE` (default) reproduces the pre-existing
`sigma` exactly, as expected. A random-slope model
(`(Days|Subject)`) with `disp_dof_correction=TRUE` correctly warns and
falls back to the uncorrected estimate rather than applying anything
unvalidated.

**Explicitly not done in this pass** (see TODO): a full test-suite run
to check for both real bugs and reference-value drift, regression tests
for the new control option, and a NEWS entry. This was intentionally
stopped at a working checkpoint rather than pushed through to a fully
validated, tested state in one sitting.

### `randomslope10`: `q_eff` generalizes to a correlated random-slope term (2026-08-07)

Directly tests the working hypothesis above (TODO item, and the
"spherical random effects" framing in the `disp_dof_correction`
scalar-RE-only gate): does `q_eff = sum(q_k) - (K-1)` still hold when a
single term (`K=1`) has more than one column per level, rather than
requiring `K` separate scalar-intercept terms?

New case `randomslope10` in `misc/Gamma_GLMM/paramsurvey_loggaussian/`:
reuses `oneway10`'s 30-row design (first 10 `Rail1` levels, 3 reps) and
adds a small-variance covariate `x = rnorm(n, 0, 0.1)`, fit with
`travel ~ 1 + (1 + x | Rail1)` -- a single correlated random-intercept-
and-slope term, true intercept/slope SD = 1, true correlation = 0.1.
`q` for this term is `nlevels(Rail1) * 2 = 20` (2 columns/level); with
`K=1`, the formula predicts `q_eff = q = 20` (no extra `-1`, since
that only applies *between* multiple terms sharing the fixed
intercept). `mgcv` is skipped here -- its `s(g, x, bs="re")` gives a
slope-only random effect, not the correlated intercept+slope structure
this case needs, and no other easy mgcv equivalent was found.

Getting the *true* covariance right for `glmmTMB::simulate_new()` took
care: glmmTMB's `theta` for an unstructured ("us") covariance term is
**not** the same convention as lme4's (`reTrms$theta`, a raw relative-
Cholesky-factor of the covariance itself) -- glmmTMB's first `n`
elements are log-SDs, the remaining `n(n-1)/2` are a *scaled* Cholesky
term of the *correlation* matrix, not the correlation directly (see
glmmTMB's
[covstruct vignette](https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html#mappings)).
For SD=1, SD=1, correlation=`rho`: `theta = c(0, 0, rho/sqrt(1-rho^2))`.
Caught before implementing (not after, this time) by asking rather than
assuming, given the earlier `betadisp` mapping mistake in the crossed-RE
work above.

**Result, B=10, both families**: all 10 glmmTMB and glmer fits clean for
both gaussian and Gamma (`q=20` confirmed identical across methods, as
usual). The `glmer/glmmTMB` phi ratio is -- again -- essentially
deterministic given the fixed design: gaussian mean ratio 0.33347 (sd
1.6e-5) against a predicted `(n-q)/n = (30-20)/30 = 0.3333`, giving
implied `q_eff = 19.996`; Gamma mean ratio 0.33353 (sd 3.8e-3), implied
`q_eff = 19.994`. Both essentially exactly `20`, confirming the
prediction with no extra correction needed beyond treating the term's
full raw (spherical) dimension as `q`. (The correlation estimates
themselves are noisy across the B=10 replicates, as expected with only
10 groups -- correlation parameters are notoriously hard to pin down
with few groups -- but that's orthogonal to the dispersion-scale `q_eff`
check, which doesn't depend on how well `corr` itself is recovered.)

This is a real, if single-case, confirmation that `q_eff` generalizes
beyond scalar-intercept terms when measured on the spherical `u` scale,
supporting (but not yet fully justifying -- nested designs and
multi-term correlated models are still untested) relaxing
`computeQEff()`'s current scalar-RE-only gate in a future pass.

### The naive formula breaks for nesting and shared-covariate crossed slopes (2026-08-07)

Two more cases, `nested5` and `crossedslopes5`, both on the same 75-row
5x5x3 `structured5` design (full crossing, first 5 levels of each
factor):

- `nested5`: `~1 + (1|Rail1/Rail2)` (nested, not crossed), true SD=1 for
  both `Rail1` and `Rail1:Rail2`. `q1=nlevels(Rail1)=5`,
  `q2=nlevels(Rail1:Rail2)=25`, `K=2`, naive prediction `q_eff=29`.
- `crossedslopes5`: `~1 + (1+x|Rail1) + (1+x|Rail2)`, same per-term
  parameters as `randomslope10` (SD=1, SD=1, corr=0.1) for *both*
  terms, `x` shared across both (same covariate, not two different
  ones). `q1=2*5=10`, `q2=2*5=10`, `K=2`, naive prediction `q_eff=19`.

**Both broke the naive formula, in both families:**

| case | naive prediction | observed q_eff (gaussian / Gamma) |
|---|---|---|
| nested5 | 29 | **24.99 / 24.98** |
| crossedslopes5 | 19 | **17.99 / 17.94** |

**Why, in both cases -- the real mechanism is the *rank* of the design
matrix you'd get treating every RE term as a saturated fixed effect,
not a simple one-redundancy-per-term count:**

- Nested: `Rail1`'s 5 dummy columns are each exactly reconstructible by
  summing the `Rail1:Rail2` columns that fall inside them
  (`Rail1_i = sum_j Rail1:Rail2_{i,j}`). `Rail1`'s whole block is
  redundant with the finer nested term, not just "one direction" -- it
  contributes **zero** extra rank beyond `Rail1:Rail2` (25 levels)
  minus the one usual redundancy with the intercept:
  `q_eff = 1 + (25-1) = 25`, matching observed ~24.99 almost exactly.
  The outer level of a nested design costs nothing beyond the finest
  level.
- Crossed slopes: beyond the usual one redundancy per term's intercept
  column (both sum to the shared fixed intercept), the two terms'
  *slope* columns are also redundant with **each other**:
  `sum(Rail1's x-slope columns) = x = sum(Rail2's x-slope columns)`,
  even with no explicit fixed `x` term to blame it on -- the two
  random-slope blocks are linearly dependent on each other directly.
  One extra redundancy beyond the naive count: `20 - 2 = 18`, matching
  observed ~17.94-17.99.

**This mattered beyond the paramsurvey**: `computeQEff()`'s previous
gate ("every RE term is scalar") would have let a nested model like
`(1|school/class)` -- both terms scalar -- through to the naive formula
and gotten a silently-wrong (too-large) `q_eff`. Real correctness gap
in already-shipped (if experimental/opt-in) code, not just a
paramsurvey curiosity.

### Structural safety gate replaces the "all-scalar" gate (2026-08-07)

Rather than a blanket "every term must be scalar" gate (which was both
too conservative -- it blocked the now-validated single-random-slope
case -- and not conservative enough -- it let nested models through to
a formula known to be wrong for them), `R/modular.R` gained a new
standalone `qEffUnsafeReason(reTrms)`: returns `NULL` if the naive
formula is safe to apply, or a reason string otherwise. Two checks,
directly targeting the two failure modes just found:

1. **Nesting** between any pair of RE terms, via
   `reformulas::isNested()` (a fast, already-existing check -- not
   written new for this) on every pair of terms' grouping factors.
2. **Shared slope covariates** across more than one multi-column term,
   via `reTrms$cnms` (already lists every slope covariate per term;
   just checking for a name duplicated across terms). A *single*
   multi-column term (`K=1`) is unaffected by this check -- no second
   term for its slope to be redundant with, and that case is validated
   separately (spot-checked on `sleepstudy`, `(Days|Subject)`,
   `disp_dof_correction=TRUE`: sigma ratio 1.1228 vs. predicted
   `sqrt(n/(n-q))=1.1180` for `n=180, q=36`).

`computeQEff()` now calls this once and only falls back (with a
`warning()` naming the specific reason) when it fires -- otherwise
applies the same `sum(q_k)-(K-1)` formula as before, now also correctly
covering the single-random-slope-term case that the old gate
incorrectly blocked.

**Not checked** (no paramsurvey case yet, but consistent with the same
rank argument and left enabled rather than blocked): a mix of scalar
and multi-column terms where the multi-column term's slope covariate
isn't shared with any other term, e.g. `(1+x|Rail1)+(1|Rail2)`.

**Deliberately not attempted in this pass**: a fully general numerical
rank computation (which would catch redundancy patterns beyond these
two known ones, without needing them enumerated by hand) -- see TODO.
The two structural checks here are cheap (`O(n)`/`O(q)` factor-level
bookkeeping, no linear algebra) but only catch what they're explicitly
checking for.

### Replaced by a real rank computation -- structural gate retired (2026-08-07)

The TODO item above ("what if some other redundancy pattern shows up
that neither structural check catches?") turned out to be answerable
directly, and cheaply. `computeQEff()` no longer uses the naive
`sum(q_k)-(K-1)` formula plus a hand-maintained list of known-bad
structures at all -- it computes `q_eff` as the actual rank of the
combined "saturated" design `[X, Z]` (fixed-effect columns plus every
random-effects term's raw dummy/product columns, exactly the classical
`RSS/(n-p)` logic generalized to whatever redundancy the random effects
happen to have), via `Matrix::rankMatrix(cbind(X, t(Zt)), method =
"qr")`. `qEffUnsafeReason()` and the whole structural-detection gate
above are gone -- nesting and shared-slope-covariates aren't special
cases to detect anymore, they just fall out of the rank computation
like any other redundancy would.

**Why `method = "qr"` specifically**: `rankMatrix()`'s default method
(`"tolNorm2"`) calls `svd(x, 0, 0)` internally, which has no sparse
path and silently densifies (with an explicit warning past a size
threshold: `"rankMatrix(<large sparse Matrix>, ...) coerces to dense
matrix"`). `method="qr"` calls `qr(x, ...)` directly, which *does*
dispatch to a real sparse QR for `Matrix` sparse classes, confirmed via
direct comparison at `q=1000` (two 500-level crossed factors, `n=5000`):
`method="qr"` took 0.72s and stayed sparse; the default method took
3.96s after densifying. 0.72s is already *less* than that model's own
fit time (0.83s) -- and the rank only needs computing once per fit (the
design's rank doesn't depend on theta/phi/beta, unlike everything else
`profilePhi()` iterates on), so in context it's a rounding error. At
the much smaller paramsurvey scale (`n<=108`) the rank computation cost
was consistently under 0.1% of the corresponding fit's own time.

**Verified against every case already in this document**: the rank
computation reproduces `structured5`'s `q_eff~9`, `nested5`'s `~25`,
and `crossedslopes5`'s `~18` exactly (computed once from the actual
design matrices, not re-run through the full B=10 paramsurvey sweep --
this is a direct structural check, not a new simulation study). Spot-
checked live on `glmer()` fits too (Gamma family): the nested and
crossed-slopes-sharing-`x` cases above, which the structural gate used
to decline with a warning, now get a real correction applied
automatically, with observed/predicted sigma ratios close (nested:
1.234 observed vs. 1.225 predicted; crossed-slopes: 1.149 vs. 1.147 --
single-dataset spot checks, not averaged, so not expected to match to
the same precision as the B=10 paramsurvey numbers elsewhere in this
doc). Full test suite re-run clean after this change
(`LME4_TEST_LEVEL=2`).

### OLRE/OOM-avoidance regression test added (2026-08-07)

Closes the last untested piece of `computeQEff()`'s safety net (see
"First implementation" above): `tests/testthat/test-sigma-dof.R` gained
two regression tests exercising the `tryCatch` fallback around the rank
computation -- the same code path that catches the "QR factorization ...
out of memory" failure observed at `q~50000` for very large random-effect
dimensions. Rather than reproducing an actual `q~50000` failure (slow,
expensive, and not portable across CI machines), `testthat::
local_mocked_bindings()` forces `rankMatrix()` to fail deterministically
and cheaply -- possible because `rankMatrix` is imported directly into
lme4's own namespace (`importFrom("Matrix", ..., rankMatrix, ...)`,
`NAMESPACE`), which is exactly what `local_mocked_bindings()` needs to
intercept a call inside `computeQEff()`.

Uses a real observation-level-random-effect (OLRE) Gamma structure (one
RE level per observation, `q=n=200`), built via `glFormula()` rather than
a full `glmer()` fit: an actual OLRE Gamma model turned out to be
numerically pathological enough (each level has exactly one observation,
so within-level variance is fully confounded with the residual) that
`glmer()` errors outright during PIRLS ("Downdated VtV is not positive
definite") -- unrelated to `computeQEff()` or the dispersion correction
at all, just a genuinely ill-posed fit, consistent with the standing note
that an OLRE on a free-dispersion model is probably a modeling mistake to
begin with. `glFormula()` only builds the design matrices, so it's
unaffected by that and gives a stable structure to test `computeQEff()`
against directly. A companion unmocked test confirms the real computation
succeeds normally on the same structure and returns `qEff == n` exactly
-- concretely illustrating *why* OLRE plus a free-dispersion family is a
degenerate combination for this correction specifically: the correction's
own denominator (`n - qEff`) would be exactly zero.

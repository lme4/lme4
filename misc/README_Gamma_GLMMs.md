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
not a separate unexplained phenomenon. Root cause not yet re-derived for
this side (the §7 mechanism was worked out for the disp<1/inflation
direction); worth revisiting given how large and consistent the effect is.
Locked in as a characterization test in
`tests/testthat/test-gamma_glmm_bias.R` (gated behind
`LME4_TEST_LEVEL > 1`), alongside the disp=0.05 inflation case, so both
directions are covered.

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

## Practical takeaway (regardless of whether/when this gets fixed upstream)

When advising on Gamma GLMMs (or other estimated-dispersion families) in
`lme4`: warn that `glmer`'s default Laplace fit can substantially inflate
random-effect variance estimates, especially with more than one grouping
factor (where `nAGQ>1` isn't available as a cross-check, since it's
restricted to single-scalar-RE models). Cross-check with `glmmTMB` when the
precision of variance components matters. The bias is worst when the true
dispersion is small (shape large) — i.e. in the regime most real Gamma GLMM
applications actually live in.

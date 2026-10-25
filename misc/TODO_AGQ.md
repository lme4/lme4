# AGQ for vector-valued random effects: open issues

Branch `agq-vector`. Done so far: `glmer(nAGQ > 1)` works for a single random-effects term of
any width (`(1+x|g)`, `diag(x|g)`, `cs(x|g)`, ...). Quadrature rules come from `AGQrule()`
(`R/GHrule.R`): a tensor product (optionally pruned) or a sparse grid. The C++ entry point is
`glmerAGQ()` (`src/external.cpp`). Scalar fits are bit-identical to master.

## Open issues

### Correctness / validation
- [ ] **Free-dispersion families (Gamma etc.): constant offset in the deviance.** At fixed
  (theta, beta), the glmer devfun converges as nAGQ grows, but to a value about 15 below
  -2 log L computed by brute-force quadrature with the same profiled phi (2-D `(x|g)` Gamma-log
  example). Master shows the same kind of offset for scalar `(1|g)` (about -6.7), so this is
  existing behaviour, not an AGQ bug. Likely causes: `Gamma()$aic` adds +2 and plugs in its own
  dispersion, dev/n, instead of the profiled phi, which uses dev/(n - qEff). This needs sorting
  out before AGQ log-likelihoods for these families can be compared with other software.
- [ ] Validate Gamma/inverse-Gaussian vector AGQ against an independent reference once the
  offset above is understood. So far only binomial and Poisson are checked against the
  brute-force reference and GLMMadaptive.
- [ ] Check binomial models with `cbind(successes, failures)` responses and prior weights.
- [ ] Check models with offsets, and `nAGQ0initStep = FALSE`.
- [ ] Check sparse grids on strongly non-Gaussian integrands (small clusters, binary data). With
  negative weights the per-level sum can become non-positive. That currently returns the 1e10
  sentinel deviance, which the optimizer could get stuck against; consider a warning, or
  falling back to the product rule.

### Scope extensions
- [ ] **Several terms with the same single grouping factor**, e.g. `(1|g) + (0+x|g)`, including
  `(x||g)`, which expands into separate terms. The per-level blocks still factor by level, but
  each level's u entries are not contiguous. This needs a level-to-indices map, built from
  `reTrms$Gp` / `cnms`, in place of the `(j-1)k + 1..k` assumption.
- [ ] **Nested terms, `(1|a/b)`** (stretch goal, separate `agq-nested` branch), following
  Pinheiro & Chao (2006, JCGS 15:58) and Rabe-Hesketh et al. (2005):
  - Per top-level cluster, the precision has an arrowhead structure. Use the Schur complement
    for u_a, and condition u_b on u_a.
  - Nested sum: over z_a, then a product over the b levels of the sums over z_b. Cost is
    n_a x n_b linear-predictor evaluations.
  - Separate function `glmerAGQnested()`.
- [ ] Crossed terms: out of scope (the integral doesn't factor).
- [ ] `nlmer(nAGQ > 1)`: still unimplemented (stub in `R/lmer.R`).

### Performance
- [ ] The rule is rebuilt (`mkAGQmat()`) on every `refit()`, so `bootMer()` on a large grid
  repeats the > 1000-node warning once per simulation. Consider storing the rule in the fit,
  or warning only once.
- [ ] Each node costs a full linear-predictor evaluation (`pp->linPred`). Updating only Z
  Lambda delta would be cheaper, and the per-level `R_j^{-1} z` products could be vectorized.
- [ ] Explore other pruning rules for large k (e.g. dropping nodes by total weight, as in
  MixedModels.jl), and adaptive choice of nAGQ.

### Interface / docs
- [ ] Decide whether to export `AGQrule()` and document it (currently internal; `GHrule()` and
  `GQdk()` are exported).
- [ ] The scaling uses the expected (Fisher/PIRLS) Hessian rather than the observed Hessian, as
  GLMMadaptive does. These are equivalent for canonical links. For non-canonical links this
  only affects efficiency at finite nAGQ, but it should be documented, and maybe made optional.
- [ ] `profile()`, `confint(method = "profile")`, `allFit()` and `drop1()`/`anova()` with vector
  AGQ fits are untested.
- [ ] Add a worked example (e.g. a random-slope binomial model) to `?glmer` or the glmer
  vignette, comparing nAGQ = 1, 5 and 9.

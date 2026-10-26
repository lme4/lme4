# AGQ for vector-valued random effects: open issues

Branch `agq-vector`. Done so far: `glmer(nAGQ > 1)` works for a single random-effects term of
any width (`(1+x|g)`, `diag(x|g)`, `cs(x|g)`, ...). Quadrature rules come from `AGQrule()`
(`R/GHrule.R`): a tensor product (optionally pruned) or a sparse grid. The C++ entry point is
`glmerAGQ()` (`src/external.cpp`). Scalar fits are bit-identical to master.

## Open issues

### Correctness / validation
- [ ] **Free-dispersion families (Gamma etc.): the deviance uses two different phi values.** At
  fixed (theta, beta), the glmer devfun converges as nAGQ grows, but to a value about 15 below
  -2 log L computed by brute-force quadrature with the same profiled phi (2-D `(x|g)` Gamma-log
  example). Master shows the same kind of offset for scalar `(1|g)` (about -6.7), so this is
  existing behaviour, not an AGQ bug.
  - **Cause:** since `disp_dof_correction=TRUE` became the default (`8e1e24e1`), PIRLS, the
    conditional modes, `L` and `sigma()` use the profiled phi = dev/(n - qEff)
    (`profilePhi()`, `src/external.cpp`). `gammaDist::aic()` (`src/glmFamily.cpp`) still plugs in
    its own dev/n. The `+2` is not the cause: `glmResp::Laplace()` already subtracts it
    (`db6d172c`).
  - **Evidence** (`test-gamma_glmm_bias.R` example, n=200, 20 groups, `(1|f)`): with
    `disp_dof_correction=FALSE`, the reported -2 log L, the recomputed Laplace value and brute force
    agree to about 0.01. With the default `TRUE`, the reported value is 772.95, against 776.38
    from brute force at the profiled phi and 773.51 at dev/n. The gap grows with qEff/n.
  - **Why earlier checks missed it:** dev/n is close to the likelihood-maximising phi, so the
    reported value stays close to glmmTMB's maximum (772.84). Laplace and AGQ share the same
    `aic()`, so their agreement doesn't test it. The regression test in `test-gamma_glmm_bias.R`
    recomputes with `phihat <- dev/n`, not `sigma(fit)^2`, so it locks in the plug-in. That test
    was written before `8e1e24e1`, when the two phi values were the same.
  - **Effect on estimates:** the devfun is the optimizer's objective. Minimizing a consistent
    objective (density at the profiled phi) in the example above moves theta 0.1541 -> 0.1575
    (about 2%); beta and phi barely change, and -2 log L goes up by about 3.4. So the main effect
    is on `logLik`/`AIC`/`BIC`/`anova()`/profiles, but a fix will slightly change theta and the
    test reference values. Not yet checked for small-cluster designs, where qEff/n is larger.
  - **Neither phi is the ML estimate.** The n - qEff denominator is a Bessel-type
    residual-df correction. It treats the conditional modes as fitted fixed effects (qEff = rank
    of [X Z]), so it is stronger than REML, which integrates out beta only and gives n - p. Even
    dev/n is not the Gamma ML estimate: the MLE of the shape nu solves
    log(nu) - digamma(nu) = D/(2n), and D/n ~ phi(1 + phi/6), so dev/n overestimates phi.
    (In the example, glmmTMB phi = 0.5625 and 0.5625 * (1 + 0.5625/6) ~ 0.615, against
    dev/n = 0.595.) For `gaussian(link="log")`, dev/n is the conditional ML. Either way, phi comes
    from an estimating equation, not from maximizing the Laplace likelihood over phi, so a
    consistent log-likelihood at the reported phi will sit below the maximum. At finite sample
    sizes we don't know which of dev/n, dev/(n - qEff) or ML gives less bias or lower MSE for phi
    or for the other parameters; that would need a simulation study.
  - **To decide before fixing:** either pass the profiled phi to `aic()`/`Laplace()` (giving a
    consistent but non-ML log-likelihood, about 3.4 worse than glmmTMB in the example), or use the
    dof correction only for reporting `sigma()` and keep it out of the likelihood. Either way,
    change the regression test to use `sigma(fit)^2`. This needs sorting out before AGQ
    log-likelihoods for these families can be compared with other software.
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

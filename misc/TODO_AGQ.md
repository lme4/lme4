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
  Diagnosis and a simulation study of the dispersion estimators are on branch
  `GLMM_dispersion_redux` (`misc/TODO_dispersion.md`).
- [ ] Validate Gamma/inverse-Gaussian vector AGQ against an independent reference once the
  offset above is understood. So far only binomial and Poisson are checked against the
  brute-force reference and GLMMadaptive.
- [x] Binomial `cbind(s, n - s)`, `s/n` + `weights = n` and disaggregated 0/1 data agree
  (the last differs by the constant `2 * sum(lchoose(n, s))`) and converge to the brute-force
  reference; Poisson with non-integer prior weights too (`misc/agq_checks/agq_checks.R`, tests).
- [x] Offsets (formula `offset()` and `offset =` agree, match the reference; `refit()` OK) and
  `nAGQ0initStep = FALSE` (same estimates as the default, product and sparse grids).
- [x] Sparse grids on strongly non-Gaussian integrands (binary, 2-4 obs per level, k = 2-3;
  `misc/agq_checks/`). The 1e10 sentinel (non-positive per-level sum) is rare and was never hit
  during optimization. The real failure is a deviance that is far too low: devfun errors of 5-50
  that often don't shrink with nAGQ, and 2 of 72 glmer fits ran away (deviance 28-39 below the
  accurate value). Those two had mean per-level log corrections of about 0.4 against at most 0.15 for
  all other fits, so `glmerAGQcheck()` (`R/lmer.R`) warns for sparse grids when the mean exceeds
  0.2 or any sum is non-positive (from the `AGQcorr` attribute that `glmerAGQ()` now returns).
  Also, sparse grids only have fewer nodes than the product rule for k >= 4 (for k = 2 they
  have more). Docs now call them experimental.
- [x] `refit()` (and so `bootMer()`, which collects the warnings in `attr(, "boot.all.msgs")`)
  runs `glmerAGQcheck()`. `refit()` also now keeps `AGQ_grid`/`AGQ_prune`; before, a refitted
  sparse or pruned fit was recorded as a plain product rule, so a second refit changed rule.
- [ ] Consider whether sparse grids are worth keeping at all, or a better rule (e.g. nested
  Kronrod-Patterson / pruned product).

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
- [x] `profile()`, `confint()` (profile, Wald, boot), `allFit()`, `drop1()` and `anova()` work
  with vector AGQ fits (test in test-agq-vector.R). `anova()` now warns when the models were
  fitted with different nAGQ.
- [x] Worked example in `?glmer` (`\donttest{}`): Contraception (mlmRev),
  `(1 + urban | district)`, nAGQ = 1, 5, 9; about 9 s.
- [ ] Not AGQ-specific: `allFit()` fails ("object 'd' not found") when the model's data are
  visible only inside a function (also in CRAN lme4 2.1.1). It refits with
  `update(object, control = ctrl)` and puts `ctrl` in an environment attached to a *copy* of the
  formula, so `update.merMod()` can't find both `ctrl` and the data in one environment. Likely
  fix: inline the values, `bquote(update(object, control = .(ctrl)))` (and `start = .(pars)`).

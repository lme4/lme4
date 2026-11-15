# Free-dispersion GLMMs (Gamma etc.): dispersion estimate and log-likelihood

Branch `GLMM_dispersion_redux`. `glmerControl(disp_dof_correction=)` accepts `"n"`, `"n-p"`
or `"n-q"`; the simulation study is in `misc/Gamma_GLMM/dispsim/` (see its README).

## Open issues

- [ ] **The deviance uses two different phi values.** At fixed (theta, beta), the glmer
  devfun converges as nAGQ grows, but not to -2 log L computed by brute-force quadrature with
  the same profiled phi: about -6.7 for scalar `(1|g)` on master, and about -15 for a 2-D
  `(x|g)` Gamma-log example with the vector-valued AGQ of branch `agq-vector` (PR #1015).
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
  - **Simulation study:** `misc/Gamma_GLMM/dispsim/` (see its README; `glmerControl(
    disp_dof_correction = "n" / "n-p" / "n-q")`). For Gamma `(1|g)`, phi in {0.2, 1, 2} x
    cluster size {4, 20}: dev/n is biased low at small phi with large q/n (-17%) and high at
    large phi with small q/n (+24%), which explains the conflicting paramsurvey results. n - q
    over-corrects in every cell (+4% to +57%). Post-hoc digamma keeps the overfitting bias
    (-19% in small clusters). Only glmmTMB ML/REML is unbiased for phi everywhere, but glmer
    dev/n has the best RMSE for the RE sd in small clusters. Suggests not keeping n - q as the
    default; an effective-df (hat-matrix trace) correction is untried.
  - **To decide before fixing:** either pass the profiled phi to `aic()`/`Laplace()` (giving a
    consistent but non-ML log-likelihood, about 3.4 worse than glmmTMB in the example), or use the
    dof correction only for reporting `sigma()` and keep it out of the likelihood. Either way,
    change the regression test to use `sigma(fit)^2`. This needs sorting out before AGQ
    log-likelihoods for these families can be compared with other software.

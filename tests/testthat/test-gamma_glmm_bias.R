## Regression tests for the Gamma GLMM random-effect variance bias fix
## (GH#643). See misc/README_Gamma_GLMMs.md (branch Gamma_GLMM) for the full
## writeup and root-cause analysis: glmer's PIRLS/Laplace working weights
## used to be computed disp-free (as if the Gamma dispersion phi were always
## 1), which biased random-effect variance estimates -- upward when the true
## dispersion was small (shape > 1), and downward/toward singular collapse
## when the true dispersion was large (shape < 1). Fixed effects and the
## dispersion parameter itself were essentially unbiased in both regimes.
##
## The fix (src/respModule.{h,cpp}, src/external.cpp: glmResp gets a phi
## multiplier applied to the working weights, profiled via a nested
## fixed-point loop in glmerLaplace() using the moment estimator
## phi = deviance/n) is unconditional for the Gamma family -- there is no
## opt-out. These tests lock in that the fix keeps working, i.e. that both
## the inflation and collapse regimes stay close to unbiased. If these
## assertions start failing, that means the fix regressed, not that the
## test needs loosening.
##
## Not lightweight (each fits O(30) GLMMs), so gated behind
## LME4_TEST_LEVEL > 1.

testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1

if (testLevel > 1) {

  gammaBiasSim <- function(disp, sigma_true = 0.3, mean_true = 1,
                            ngrp = 30, nrep = 20, B = 30, seed = 2026) {
    set.seed(seed)
    dd0 <- expand.grid(group = factor(1:ngrp), rep = 1:nrep)
    est <- numeric(B)
    nsingular <- 0
    for (b in 1:B) {
      b_re <- rnorm(ngrp, 0, sigma_true)
      mu <- exp(mean_true + b_re[dd0$group])
      dd0$y <- rgamma(nrow(dd0), shape = 1 / disp, scale = mu * disp)
      m <- suppressWarnings(
        glmer(y ~ (1 | group), data = dd0, family = Gamma(link = "log"))
      )
      est[b] <- sqrt(as.numeric(VarCorr(m)$group))
      if (isSingular(m)) nsingular <- nsingular + 1
    }
    list(est = est, mean_est = mean(est), pct_bias = 100 * (mean(est) - sigma_true) / sigma_true,
         nsingular = nsingular, B = B, sigma_true = sigma_true)
  }

  test_that("Gamma GLMM RE variance is not inflated when shape > 1 (disp < 1)", {
    ## shape = 1/disp = 20 -- matches misc/GH643_disp_scan.R's disp=0.05 case;
    ## pre-fix bias here was roughly +100-120%. Observed post-fix bias is
    ## around -2%; use a loose bound to avoid a flaky test.
    r <- gammaBiasSim(disp = 0.05)
    expect_lt(abs(r$pct_bias), 30)
  })

  test_that("Gamma GLMM RE variance does not collapse when shape < 1 (disp > 1)", {
    ## shape = 1/disp = 0.5 -- pre-fix bias here was roughly -60-70%, with a
    ## majority of fits singular. Observed post-fix bias is around +10%, with
    ## no singular fits; use loose bounds to avoid a flaky test.
    r <- gammaBiasSim(disp = 2)
    expect_lt(abs(r$pct_bias), 30)
    expect_lt(r$nsingular, r$B * 0.3)
  })

}

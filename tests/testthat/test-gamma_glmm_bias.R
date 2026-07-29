## Known-bias characterization tests for Gamma GLMM random-effect variance
## estimation (GH#643). See misc/README_Gamma_GLMMs.md (branch Gamma_GLMM)
## for the full writeup and root-cause analysis: glmer's PIRLS/Laplace
## working weights are computed disp-free (as if the Gamma dispersion phi
## were always 1), which biases random-effect variance estimates -- upward
## when the true dispersion is small (shape > 1), and downward/toward
## singular collapse when the true dispersion is large (shape < 1). Fixed
## effects and the dispersion parameter itself are essentially unbiased in
## both regimes.
##
## These tests intentionally LOCK IN the current (buggy) behavior, so that
## if/when the fix described in misc/README_Gamma_GLMMs.md #7 is
## implemented, these assertions will start failing -- that's the signal to
## update or remove this file, not a regression to chase.
##
## Not lightweight (each fits O(30-60) GLMMs), so gated behind
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

  test_that("Gamma GLMM RE variance is inflated when shape > 1 (disp < 1)", {
    ## shape = 1/disp = 20 -- matches misc/GH643_disp_scan.R's disp=0.05 case
    r <- gammaBiasSim(disp = 0.05)
    ## known bias is roughly +100-120%; use a looser bound to avoid a flaky test
    expect_gt(r$pct_bias, 40)
  })

  test_that("Gamma GLMM RE variance collapses/deflates when shape < 1 (disp > 1)", {
    ## shape = 1/disp = 0.5
    r <- gammaBiasSim(disp = 2)
    ## known bias is roughly -60-70%, with a majority of fits singular
    expect_lt(r$pct_bias, -30)
    expect_gt(r$nsingular, r$B * 0.3)
  })

}

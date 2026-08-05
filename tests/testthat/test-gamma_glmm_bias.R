## Regression tests for the Gamma GLMM random-effect variance bias fix
## (GH#643). See misc/Gamma_GLMM/README_Gamma_GLMMs.md (branch Gamma_GLMM) for the full
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
    ## shape = 1/disp = 20 -- matches misc/Gamma_GLMM/disp_scan.R's disp=0.05 case;
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

## Regression test for a second, independent bug (unrelated to the RE-variance
## bias above, but found while investigating it): glmResp::Laplace()
## (src/respModule.cpp) used to add the family's aic() value straight into the
## deviance criterion. Base R's family$aic() (Gamma, gaussian, inverse.gaussian)
## bakes in an additive "+2" bookkeeping constant meant to be cancelled by the
## "p - aic/2" trick in stats::logLik.glm(); lme4's Laplace() never applied that
## cancellation, so logLik()/AIC()/BIC() for any GLMM with a freely-estimated
## dispersion parameter was off by a constant (logLik too large by 1, i.e.
## deviance too small by 2). Fixed by subtracting 2 in Laplace() when
## hasFreeDispersion() is true.
##
## This locks the fix in by recomputing -2*logLik from scratch -- straight from
## dgamma()/dnorm(), never touching the family's aic() method -- using the
## fitted model's own converged conditional modes/means, and checking it
## matches logLik(fit) exactly (not just up to an additive constant).
test_that("Laplace deviance for free-dispersion GLMMs has no aic() bookkeeping leak", {
  ldL2 <- function(fit) {
    L <- getME(fit, "L")
    ld <- 2 * Matrix::determinant(L, logarithm = TRUE, sqrt = TRUE)$modulus
    as.numeric(ld)
  }

  set.seed(101)
  n <- 200
  dd <- data.frame(x = rnorm(n), f = factor(rep(1:20, each = 10)))
  dd$y <- rgamma(n, shape = 2, scale = exp(1 + 0.5*dd$x)/2)
  gfit <- glmer(y ~ x + (1|f), data = dd, family = Gamma(link = "log"))

  mu <- fitted(gfit)
  dev <- sum(Gamma()$dev.resids(dd$y, mu, 1))
  phihat <- dev / n
  recomputed <- -2 * sum(dgamma(dd$y, shape = 1/phihat, scale = mu*phihat, log = TRUE)) +
    sum(getME(gfit, "u")^2) + ldL2(gfit)
  expect_equal(recomputed, -2*as.numeric(logLik(gfit)), tolerance = 1e-6)

  set.seed(202)
  dd$y <- pmax(rnorm(n, exp(1 + 0.3*dd$x), 0.3*exp(1 + 0.3*dd$x)), 0.01)
  ggfit <- glmer(y ~ x + (1|f), data = dd, family = gaussian(link = "log"))

  mu <- fitted(ggfit)
  dev <- sum(gaussian()$dev.resids(dd$y, mu, 1))
  phihat <- dev / n
  recomputed <- -2 * sum(dnorm(dd$y, mu, sqrt(phihat), log = TRUE)) +
    sum(getME(ggfit, "u")^2) + ldL2(ggfit)
  expect_equal(recomputed, -2*as.numeric(logLik(ggfit)), tolerance = 1e-6)
})

## Regression test for the nAGQ>1 caveat noted in
## misc/Gamma_GLMM/README_Gamma_GLMMs.md: glmerAGQ() (the nAGQ>1 code path)
## never profiles phi at all (phi gets frozen at whatever the nAGQ=0 init
## step left it, and the AGQ marginal-likelihood formula itself is a
## deviance-based approximation that isn't even the right kind of quantity
## for a free-dispersion family's density). Until that's fixed, glmer()
## warns instead. This locks in that the warning fires exactly when it
## should: nAGQ>1 + free dispersion (Gamma) yes; nAGQ<=1 no; nAGQ>1 with a
## fixed-dispersion family (binomial) no.
test_that("nAGQ>1 warns for free-dispersion families but not fixed-dispersion ones or nAGQ<=1", {
  warns <- function(expr) {
    w <- character(0)
    withCallingHandlers(force(expr),
                         warning = function(cnd) {
                           w <<- c(w, conditionMessage(cnd))
                           invokeRestart("muffleWarning")
                         })
    w
  }

  set.seed(303)
  n <- 60
  dd <- data.frame(f = factor(rep(1:15, each = 4)))
  b_re <- rnorm(15, 0, 0.3)
  dd$y <- rgamma(n, shape = 5, scale = exp(1 + b_re[dd$f]) / 5)
  dd$s <- rbinom(n, size = 10, prob = 0.4)

  msg <- "nAGQ>1 handles GLMMs with estimated dispersion parameters"

  w1 <- warns(glmer(y ~ (1 | f), data = dd, family = Gamma(link = "log"), nAGQ = 5))
  expect_true(any(grepl(msg, w1, fixed = TRUE)))

  w2 <- warns(glmer(y ~ (1 | f), data = dd, family = Gamma(link = "log"), nAGQ = 1))
  expect_false(any(grepl(msg, w2, fixed = TRUE)))

  w3 <- warns(glmer(cbind(s, 10 - s) ~ (1 | f), data = dd, family = binomial, nAGQ = 5))
  expect_false(any(grepl(msg, w3, fixed = TRUE)))
})

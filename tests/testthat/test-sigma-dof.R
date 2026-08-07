## Regression tests for glmerControl(disp_dof_correction=TRUE): a
## degrees-of-freedom correction to glmer's moment-based dispersion
## estimator, phi = deviance/(n - q_eff) instead of the plain
## deviance/n, where q_eff is the rank of the combined fixed+random
## "saturated" design [X, Z] (computeQEff(), R/modular.R, via
## Matrix::rankMatrix(cbind(X, t(Zt)), method="qr")). See
## misc/Gamma_GLMM/README_Gamma_GLMMs.md (branch Gamma_GLMM) for the
## full derivation, including the two structures (nested terms, and
## multiple random-slope terms sharing a covariate) where an earlier
## naive closed-form guess (q_eff = sum(q_k) - (K-1)) turned out to be
## wrong -- both covered here as regression cases, alongside the
## originally-validated crossed-scalar-intercept and single-random-
## slope-term cases.
##
## One simulated dataset per random-effects structure (not a full
## parameter-recovery sweep -- see misc/Gamma_GLMM/paramsurvey_loggaussian/
## for that), checking that turning the correction on moves sigma by
## the theoretically expected factor sqrt(n/(n-q_eff)). This is a
## deterministic, closed-form relationship (predicted, not just
## observed) if disp_dof_correction only rescaled sigma post-hoc; in
## fact it also perturbs PIRLS's working weights during mode-finding
## (profilePhi() uses n-q_eff as the moment-loop denominator at every
## outer iteration, not just once at the end), so the observed ratio is
## only approximately, not exactly, sqrt(n/(n-q_eff)) -- verified
## empirically to be within 1% for every case below, so a 2% tolerance
## is safe without being flaky.

test_that("disp_dof_correction moves sigma by the expected factor", {
    checkDofRatio <- function(formula, data, family = gaussian(link = "log")) {
        n <- nrow(data)
        gm <- glFormula(formula, data = data, family = family)
        qEff <- unname(Matrix::rankMatrix(cbind(gm$X, t(gm$reTrms$Zt)), method = "qr")[1])
        fit0 <- suppressWarnings(glmer(formula, data = data, family = family))
        fit1 <- suppressWarnings(glmer(formula, data = data, family = family,
                                        control = glmerControl(disp_dof_correction = TRUE)))
        observed <- sigma(fit1) / sigma(fit0)
        predicted <- sqrt(n / (n - qEff))
        expect_equal(observed, predicted, tolerance = 0.02)
    }

    set.seed(1)

    ## single scalar random intercept: q_eff = nlevels(g) = 8
    d1 <- data.frame(g = factor(rep(1:8, each = 4)))
    d1$y <- simulate(~ 1 + (1 | g), newdata = d1,
                      newparams = list(beta = 2, theta = 1, sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]
    checkDofRatio(y ~ 1 + (1 | g), d1)

    ## single correlated random-slope term (K=1, 2 columns/level):
    ## q_eff = q = 2*8 = 16, no "-1" (only one term)
    d2 <- data.frame(g = factor(rep(1:8, each = 4)))
    d2$x <- rnorm(nrow(d2), 0, 0.1)
    d2$y <- simulate(~ 1 + (1 + x | g), newdata = d2,
                      newparams = list(beta = 2, theta = c(1, 0.1, sqrt(1 - 0.1^2)), sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]
    checkDofRatio(y ~ 1 + (1 + x | g), d2)

    ## crossed scalar intercepts (K=2): q_eff = q1+q2-(K-1) = 4+4-1 = 7
    d3 <- expand.grid(g1 = factor(1:4), g2 = factor(1:4), rep = 1:2)
    d3$y <- simulate(~ 1 + (1 | g1) + (1 | g2), newdata = d3,
                      newparams = list(beta = 2, theta = c(1, 1), sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]
    checkDofRatio(y ~ 1 + (1 | g1) + (1 | g2), d3)

    ## two crossed random-slope terms sharing covariate x: naive formula
    ## (q1+q2-(K-1) = 8+8-1 = 15) is WRONG here -- the two terms' slope
    ## blocks are redundant with *each other*, one extra redundancy,
    ## true q_eff = 14 (checkDofRatio recomputes this via rankMatrix
    ## directly, doesn't rely on the naive formula)
    d4 <- expand.grid(g1 = factor(1:4), g2 = factor(1:4), rep = 1:2)
    d4$x <- rnorm(nrow(d4), 0, 0.1)
    theta_rs <- c(1, 0.1, sqrt(1 - 0.1^2))
    d4$y <- simulate(~ 1 + (1 + x | g1) + (1 + x | g2), newdata = d4,
                      newparams = list(beta = 2, theta = rep(theta_rs, 2), sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]
    checkDofRatio(y ~ 1 + (1 + x | g1) + (1 + x | g2), d4)

    ## nested terms (1|g1/g2): naive formula (q1+q2-(K-1) = 4+16-1 = 19)
    ## is WRONG here -- the outer term (g1) contributes zero extra rank,
    ## entirely subsumed by the inner nested term (g1:g2); true
    ## q_eff = 16 (again recomputed via rankMatrix, not the naive formula)
    d5 <- expand.grid(g1 = factor(1:4), g2 = factor(1:4), rep = 1:2)
    d5$y <- simulate(~ 1 + (1 | g1 / g2), newdata = d5,
                      newparams = list(beta = 2, theta = c(1, 1), sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]
    checkDofRatio(y ~ 1 + (1 | g1 / g2), d5)
})

test_that("disp_dof_correction=FALSE (default) leaves sigma unchanged", {
    set.seed(1)
    d <- data.frame(g = factor(rep(1:8, each = 4)))
    d$y <- simulate(~ 1 + (1 | g), newdata = d,
                     newparams = list(beta = 2, theta = 1, sigma = 0.1),
                     family = gaussian(link = "log"))[[1]]
    fit_default <- glmer(y ~ 1 + (1 | g), data = d, family = gaussian(link = "log"))
    fit_explicit_off <- glmer(y ~ 1 + (1 | g), data = d, family = gaussian(link = "log"),
                               control = glmerControl(disp_dof_correction = FALSE))
    expect_equal(sigma(fit_default), sigma(fit_explicit_off))
})

## Regression tests for glmerControl(disp_dof_correction=): a
## degrees-of-freedom correction to glmer's moment-based dispersion
## estimator, phi = deviance/(n - q_eff) instead of the plain
## deviance/n, where q_eff is the rank of the combined fixed+random
## "saturated" design [X, Z] (computeQEff(), R/modular.R, via
## Matrix::rankMatrix(cbind(X, t(Zt)), method="qr")). TRUE by default
## (GH #936); FALSE reproduces the older, uncorrected behaviour. See
## misc/Gamma_GLMM/README_Gamma_GLMMs.md (branch Gamma_GLMM) for the
## full derivation, including the two structures (nested terms, and
## multiple random-slope terms sharing a covariate) where an earlier
## naive closed-form guess (q_eff = sum(q_k) - (K-1)) turned out to be
## wrong -- both covered here as regression cases, alongside the
## originally-validated crossed-scalar-intercept and single-random-
## slope-term cases. Both TRUE and FALSE are always passed explicitly
## below (rather than relying on either being the current default), so
## these tests stay meaningful regardless of which way the default is
## set.
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

checkDofRatio <- function(formula, data, family = gaussian(link = "log")) {
    n <- nrow(data)
    gm <- glFormula(formula, data = data, family = family)
    qEff <- unname(Matrix::rankMatrix(cbind(gm$X, t(gm$reTrms$Zt)), method = "qr")[1])
    fit0 <- suppressWarnings(glmer(formula, data = data, family = family,
                                    control = glmerControl(disp_dof_correction = FALSE)))
    fit1 <- suppressWarnings(glmer(formula, data = data, family = family,
                                    control = glmerControl(disp_dof_correction = TRUE)))
    observed <- sigma(fit1) / sigma(fit0)
    predicted <- sqrt(n / (n - qEff))
    expect_equal(observed, predicted, tolerance = 0.02)
}

## Five small simulated datasets, one per random-effects structure,
## shared by both test_that() blocks below. All include an "x"
## covariate (a random-slope predictor for random-slope/crossed-slopes,
## unrelated noise otherwise) for the fixed-effect-redundancy test.
dofTestData <- local({
    set.seed(1)

    ## single scalar random intercept: q_eff = nlevels(g) = 8
    d1 <- data.frame(g = factor(rep(1:8, each = 4)))
    d1$y <- simulate(~ 1 + (1 | g), newdata = d1,
                      newparams = list(beta = 2, theta = 1, sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]
    d1$x <- rnorm(nrow(d1))

    ## single correlated random-slope term (K=1, 2 columns/level):
    ## q_eff = q = 2*8 = 16, no "-1" (only one term)
    d2 <- data.frame(g = factor(rep(1:8, each = 4)))
    d2$x <- rnorm(nrow(d2), 0, 0.1)
    d2$y <- simulate(~ 1 + (1 + x | g), newdata = d2,
                      newparams = list(beta = 2, theta = c(1, 0.1, sqrt(1 - 0.1^2)), sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]

    ## crossed scalar intercepts (K=2): q_eff = q1+q2-(K-1) = 4+4-1 = 7
    d3 <- expand.grid(g1 = factor(1:4), g2 = factor(1:4), rep = 1:2)
    d3$y <- simulate(~ 1 + (1 | g1) + (1 | g2), newdata = d3,
                      newparams = list(beta = 2, theta = c(1, 1), sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]
    d3$x <- rnorm(nrow(d3))

    ## two crossed random-slope terms sharing covariate x: naive formula
    ## (q1+q2-(K-1) = 8+8-1 = 15) is WRONG here -- the two terms' slope
    ## blocks are redundant with *each other*, one extra redundancy,
    ## true q_eff = 14
    theta_rs <- c(1, 0.1, sqrt(1 - 0.1^2))
    d4 <- expand.grid(g1 = factor(1:4), g2 = factor(1:4), rep = 1:2)
    d4$x <- rnorm(nrow(d4), 0, 0.1)
    d4$y <- simulate(~ 1 + (1 + x | g1) + (1 + x | g2), newdata = d4,
                      newparams = list(beta = 2, theta = rep(theta_rs, 2), sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]

    ## nested terms (1|g1/g2): naive formula (q1+q2-(K-1) = 4+16-1 = 19)
    ## is WRONG here -- the outer term (g1) contributes zero extra rank,
    ## entirely subsumed by the inner nested term (g1:g2); true
    ## q_eff = 16
    d5 <- expand.grid(g1 = factor(1:4), g2 = factor(1:4), rep = 1:2)
    d5$y <- simulate(~ 1 + (1 | g1 / g2), newdata = d5,
                      newparams = list(beta = 2, theta = c(1, 1), sigma = 0.1),
                      family = gaussian(link = "log"))[[1]]
    d5$x <- rnorm(nrow(d5))

    list(
        `single-scalar`     = list(data = d1, formula = y ~ 1 + (1 | g)),
        `random-slope`      = list(data = d2, formula = y ~ 1 + (1 + x | g)),
        `crossed-intercept` = list(data = d3, formula = y ~ 1 + (1 | g1) + (1 | g2)),
        `crossed-slopes`    = list(data = d4, formula = y ~ 1 + (1 + x | g1) + (1 + x | g2)),
        `nested`            = list(data = d5, formula = y ~ 1 + (1 | g1 / g2))
    )
})

test_that("disp_dof_correction moves sigma by the expected factor", {
    for (nm in names(dofTestData)) {
        entry <- dofTestData[[nm]]
        checkDofRatio(entry$formula, entry$data)
    }
})

test_that("disp_dof_correction accounts for fixed/random redundancy (rank of combined [X,Z])", {
    ## Add poly(x,3) to the FIXED effects on top of each RE structure
    ## above (not part of the true generating model -- this is purely a
    ## mechanism check, not testing recovery of a real fixed effect).
    ## Where x is *not* already part of the RE structure
    ## (single-scalar, crossed-intercept, nested), this costs the naive
    ## +3 degrees of freedom, as expected. But where x is already a
    ## random-slope covariate (random-slope, crossed-slopes), it only
    ## costs +2: that random slope's columns, summed across all levels
    ## of its grouping factor, reconstruct x exactly, so poly(x,3)'s
    ## linear component is redundant with that sum even though no two
    ## *random* terms are involved at all -- the same "shared
    ## covariate" mechanism found between crossed random-slope terms
    ## above, but this time between a fixed effect and a random slope.
    ## A formula that computed q_eff as p_fixed + f(RE structure)
    ## rather than the rank of the combined [X,Z] would have missed
    ## this and predicted +3 for those two cases too -- checkDofRatio()
    ## recomputes q_eff directly via rankMatrix() on the actual [X,Z]
    ## for the augmented formula, so it automatically reflects whichever
    ## delta is really there.
    for (nm in names(dofTestData)) {
        entry <- dofTestData[[nm]]
        formula_poly <- update(entry$formula, . ~ poly(x, 3) + .)
        checkDofRatio(formula_poly, entry$data)
    }
})

test_that("disp_dof_correction defaults to TRUE", {
    entry <- dofTestData[["single-scalar"]]
    fit_default <- glmer(entry$formula, data = entry$data, family = gaussian(link = "log"))
    fit_explicit_on <- glmer(entry$formula, data = entry$data, family = gaussian(link = "log"),
                              control = glmerControl(disp_dof_correction = TRUE))
    fit_explicit_off <- glmer(entry$formula, data = entry$data, family = gaussian(link = "log"),
                               control = glmerControl(disp_dof_correction = FALSE))
    expect_equal(sigma(fit_default), sigma(fit_explicit_on))
    expect_false(isTRUE(all.equal(sigma(fit_default), sigma(fit_explicit_off))))
})

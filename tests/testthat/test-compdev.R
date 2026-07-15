## test equivalence of glmer fits with compDev = TRUE or FALSE

local({
    dd <- data.frame(g = factor(rep(1:3, each = 3)))
    y <- simulate(~ 1 + (1|g), family = poisson, seed = 101, newdata = dd,
                  newparams = list(beta = 0, theta = 2))[[1]]
    g1 <- glmer(y ~ 1 + (1|g), data = dd, family = poisson,
                control = glmerControl(nAGQ0initStep = FALSE))
    g2 <- update(g1, control = glmerControl(compDev = FALSE, nAGQ0initStep = FALSE))

    test_that("glmer VarCorr agrees between compDev=TRUE and compDev=FALSE", {
        expect_identical(VarCorr(g1), VarCorr(g2))
    })

    test_that("glmer fixef agrees between compDev=TRUE and compDev=FALSE", {
        expect_identical(fixef(g1), fixef(g2))
    })

    test_that("compDev is propagated to devfun environment by getME", {
        dfun1 <- getME(g1, "devfun")
        dfun2 <- getME(g2, "devfun")
        expect_true( environment(dfun1)$compDev)
        expect_false(environment(dfun2)$compDev)
    })

    test_that("devfuns with compDev=TRUE and compDev=FALSE agree at a test point", {
        dfun1 <- getME(g1, "devfun")
        dfun2 <- getME(g2, "devfun")
        pars <- c(1.16204380240176, -0.183842216759865)
        expect_identical(dfun1(pars), dfun2(pars))
    })
})

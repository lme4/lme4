## GH #998: refresh the final PIRLS factorization before evaluating Laplace.

test_that("the complete-block factorial example uses the final working weights", {
    d <- expand.grid(c=c(1, -1), b=c(1, -1), a=c(1, -1), block=1:180)
    d$block <- factor(d$block)
    set.seed(271019)
    h <- rnorm(180, sd=sqrt(10))
    d$y <- as.integer(runif(nrow(d)) < plogis(h[d$block]))

    m <- glmer(y ~ (a + b + c)^2 + (1 | block), d, family=binomial)
    ac <- coef(summary(m))["a:c", ]
    ## Before the refresh, SE = 0.07746 and p = 0.08352, without a warning.
    expect_equal(unname(ac["Std. Error"]), 0.087106, tolerance=1e-3)
    expect_equal(unname(ac["Pr(>|z|)"]), 0.119015, tolerance=1e-3)

    dev <- update(m, devFunOnly=TRUE)
    value <- dev(c(getME(m, "theta"), fixef(m)))
    rho <- environment(dev)
    p <- rho$resp$mu
    ## For independent block intercepts, the determinant is a sum of
    ## log(1 + theta^2 * sum(p * (1-p))) over blocks, at the accepted means.
    ldL2 <- sum(log1p(getME(m, "theta")^2 * rowsum(p * (1-p), d$block)))
    expect_equal(rho$pp$ldL2(), ldL2, tolerance=1e-10)
    expect_equal(value, rho$resp$Laplace(ldL2, 0, rho$pp$sqrL(1)),
                 tolerance=1e-10)
})

test_that("Hessian standard errors do not collapse in the reported examples", {
    m <- glmer(cbind(incidence, size-incidence) ~ period + (1 | herd),
               cbpp, family=binomial, control=glmerControl(tolPwrss=1e-8))
    ## Use the Hessian directly: an RX fallback must not make this test pass.
    se <- sqrt(2 * diag(solve(m@optinfo$derivs$Hessian))[-1])
    expect_equal(unname(se[1]), 0.23247, tolerance=1e-3)

    ## The reporter's 300-subject example fails at the default tolerance
    ## for seed 1, and at a much tighter tolerance for seed 29.
    cases <- data.frame(seed=c(1, 29), sd=c(0.5, 1),
                        tol=c(1e-7, 1e-13), se=c(0.16295, 0.16678))
    for (i in seq_len(nrow(cases))) {
        d <- expand.grid(visit=1:7, subject=factor(1:300))
        d$treatment <- as.numeric(d$subject) %% 2L
        set.seed(cases$seed[i])
        u <- rnorm(300)
        eta <- -1 - 0.7 * d$treatment - 0.25 * d$visit +
            cases$sd[i] * u[d$subject]
        d$y <- rbinom(nrow(d), 1, plogis(eta))
        m <- glmer(y ~ treatment + visit + (1 | subject), d, family=binomial,
                   control=glmerControl(tolPwrss=cases$tol[i]))
        se <- sqrt(2 * diag(solve(m@optinfo$derivs$Hessian))[-1])
        expect_equal(unname(se[1]), cases$se[i], tolerance=1e-3)
    }
})

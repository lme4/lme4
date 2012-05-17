library("testthat")
library("lme4")

set.seed(101)
d <- expand.grid(block=LETTERS[1:26], rep=1:100, KEEP.OUT.ATTRS = FALSE)
d$x <- runif(nrow(d))  ## sd=1
reff_f <- rnorm(length(levels(d$block)),sd=1)
## set intercept large enough to create a constant response
d$eta0 <- 4+3*d$x  ## fixed effects only
d$eta <- d$eta0+reff_f[d$block]
dBc <- d
cc <- binomial(link="cloglog")
dBc$mu <- cc$linkinv(d$eta)
dBc$y <- rbinom(nrow(d),dBc$mu,size=1)

m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = binomial, data = cbpp)
context("Errors and warnings from glmer")
test_that("glmer", {
    expect_error(glmer(y ~ 1 + (1|block), data=dBc, family=binomial(link="cloglog")))
    expect_warning(m2 <-lmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             family = binomial, data = cbpp),
                   "calling lmer with family\\(\\) is deprecated.*")
    expect_equal(m1,m2)
    expect_warning(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                         family = binomial, data = cbpp, REML=TRUE),
                   "extra argument.*REML.*disregarded")
    m3 <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    m4 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    m5 <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy, family=gaussian)
    expect_equal(fixef(m3),fixef(m5))
    expect_equal(m3,m4)
    isTRUE(all.equal(m3,m5))

    ## would like m3==m5 != m4 ??
    VarCorr(m4)
    VarCorr(m5)  ## wrong??? is this the report or the
    getME(m4,"theta")
    getME(m5,"theta")
})


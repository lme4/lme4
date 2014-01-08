library("testthat")
library("lme4")

source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
##-> gSim(), a general simulation function ...

set.seed(101)
dBc <- gSim(family=binomial(link="cloglog"), nbinom = 1) # {0,1} Binomial

## m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##             family = binomial, data = cbpp)
context("Errors and warnings from glmer")
test_that("glmer", {
    expect_error(glmer(y ~ 1 + (1|block), data=dBc, family=binomial(link="cloglog")),
                 "Response is constant")
    expect_warning(lmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                              family = binomial, data = cbpp),
                   "calling lmer with .*family.* is deprecated.*")
    ## expect_equal(m1,m2)
    expect_warning(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                         family = binomial, data = cbpp, REML=TRUE),
                   "extra argument.*REML")
    expect_warning(glmer(Reaction ~ Days + (Days|Subject), sleepstudy),"calling glmer.*family=gaussian.*deprecated")
    expect_warning(glmer(Reaction ~ Days + (Days|Subject), sleepstudy, family=gaussian),
                          "calling glmer.*family=gaussian.*deprecated")
    m3 <- suppressWarnings(glmer(Reaction ~ Days + (Days|Subject), sleepstudy))
    m4 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    m5 <- suppressWarnings(glmer(Reaction ~ Days + (Days|Subject), sleepstudy, family=gaussian))
    expect_equal(fixef(m3),fixef(m5))
    m3@call[[1]] <- m5@call[[1]] <- quote(lmer)  ## hack call
    expect_equal(m3,m4)
    expect_equal(m3,m5)
})


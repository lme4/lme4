library("testthat")
library("lme4")
L <- load(system.file("testdata/lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))
fm1 <- fit_sleepstudy_1
fm0 <- fit_sleepstudy_0

context("summary")
test_that("summary", {
    ## test for multiple-correlation-warning bug
    cc <- capture.output(summary(fit_agridat_archbold))
    expect_true(length(g <- grep("not shown by default",cc))==0 || g<=1)
})
context("anova")
test_that("lmer", {
    expect_that(anova(fm0,fm1),                        is_a("anova"))
    expect_warning(do.call(anova,list(fm0,fm1)),"assigning generic names")
})

context("bootMer")
test_that("bootMer", {
    ## testing bug-fix for ordering of sd/cor components in sd/cor matrix with >2 rows
    m1 <- lmer(strength~1+(cask|batch),Pastes)
    bb <- suppressWarnings(confint(m1,method="boot",nsim=3,quiet=TRUE))
    corvals <- bb[grep("^cor_",rownames(bb)),]
    expect_true(all(abs(corvals)<=1))
})

context("confint")
test_that("confint", {
    load(system.file("testdata","gotway_hessianfly.rda",package="lme4"))
    ## gotway_hessianfly_fit <- glmer(cbind(y, n-y) ~ gen + (1|block),
    ##              data=gotway.hessianfly, family=binomial,
    ##              control=glmerControl(check.nlev.gtreq.5="ignore"))
    ## gotway_hessianfly_prof <- profile(gotway_hessianfly_fit,which=1)
    ## save(list=ls(pattern="gotway"),file="gotway_hessianfly.rda")
    ## expect_equal(confint(gotway_hessianfly_prof)[1,1],0)
    ## FIXME: should add tests for {-1,1} bounds on correlations as well
})

context("refit")
test_that("refit", {
    s1 <- simulate(fm1)
    expect_is(refit(fm1,s1),"merMod")
    s2 <- simulate(fm1,2)
    expect_error(refit(fm1,s2),"refit not implemented for lists")
})


library("testthat")
library("lme4")

context("anova")
test_that("lmer", {
    fm1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
    fm0 <- update(fm1,.~.-Days)
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
    expect_equal(confint(gotway_hessianfly_prof)[1,1],0)
    ## FIXME: should add tests for {-1,1} bounds on correlations as well
})




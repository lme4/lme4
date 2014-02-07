library("testthat")
library("lme4")
L <- load(system.file("testdata/lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))

if (getRversion()>"3.0.0") {
    ## saved fits are not safe with old R versions

fm0 <- fit_sleepstudy_0
fm1 <- fit_sleepstudy_1
gm1 <- fit_cbpp_1
gm2 <- fit_cbpp_2

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

    dat <- data.frame(y=1:5,u=c(rep("A",2),rep("B",3)),
                      t=c(rep("A",3),rep("B",2)))

    datfun <- function(x) dat
    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa <-  dat
    expect_is(stats::anova(lmer(y ~ u + (1 | t),
    dat = aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa,
                                REML=FALSE),
                           lmer(y ~ 1 + (1 | t),
     dat = aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa,
                                REML=FALSE)),"anova")
    expect_equal(rownames(stats::anova(lmer(y ~ u + (1 | t),
     dat = aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa,
                                            REML=FALSE),
                                       lmer(y ~ 1 + (1 | t),
     dat = aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa,
                                            REML=FALSE),
                                       model.names=c("a","b"))),
                 c("b","a"))
    expect_error(rownames(stats::anova(lmer(y ~ u + (1 | t),
                                            dat = dat, REML=FALSE),
                                       lmer(y ~ 1 + (1 | t),
                                            dat = dat, REML=FALSE),
                                       model.names=c("a","b","c"))),
                 "different lengths")
    z <- 1
    stats::anova(lmer(y ~ u + (1 | t),
                      data = datfun(z), REML=FALSE),
                 lmer(y ~ 1 + (1 | t),
                      data = datfun(z), REML=FALSE))

    datList <- list(data,data,data)
    
})

context("bootMer")
test_that("bootMer", {
    ## testing bug-fix for ordering of sd/cor components in sd/cor matrix with >2 rows
    m1 <- lmer(strength~1+(cask|batch),Pastes)
    bb <- suppressWarnings(confint(m1,method="boot",nsim=3,quiet=TRUE))
    corvals <- bb[grep("^cor_",rownames(bb)),]
    expect_true(all(abs(corvals)<=1))
    ## test bootMer with GLMM, multiple RE
    gm2 <- fit_cbpp_2
    expect_is(suppressWarnings(confint(gm2,method="boot",nsim=2,quiet=TRUE)),
                               "matrix")
    expect_equal(nrow(suppressWarnings(confint(gm2,
                 method="boot",nsim=2,quiet=TRUE,parm=4:6))),3)
    ## bootMer with NA values
    PastesNA <- Pastes
    PastesNA$Sample[1:3] <- NA
    m2 <- update(m1,data=PastesNA)
    confint(m2,method="boot",nsim=10,quiet=TRUE,seed=101)
    fm1 <- lmer(Reaction~Days+(Days|Subject),sleepstudy)
    sleepstudyNA <- sleepstudy
    sleepstudyNA$Days[1:3] = NA
    fm2 <- update(fm1,data=sleepstudyNA)
    confint(fm2, method="boot", nsim=10, seed=101)

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

context("refit")
test_that("refit", {
    s1 <- simulate(fm1)
    expect_is(refit(fm1,s1),"merMod")
    s2 <- simulate(fm1,2)
    expect_error(refit(fm1,s2),"refit not implemented for lists")
})

context("predict")
test_that("predict", {
    d1 <- with(cbpp, expand.grid(period=unique(period), herd=unique(herd)))
    d2 <- data.frame(period="1", herd=unique(cbpp$herd))
    d3 <- expand.grid(period=as.character(1:3),
                      herd=unique(cbpp$herd))
    p0 <- predict(gm1)
    p1 <- predict(gm1,d1)
    p2 <- predict(gm1,d2)
    p3 <- predict(gm1,d3)
    expect_equal(p0[1],p1[1],tolerance=4e-5)
    expect_equal(p0[1],p2[1],tolerance=4e-5)
    expect_equal(p0[1],p3[1],tolerance=4e-5)
    expect_message(predict(gm1,ReForm=NA),"is now preferred to")
})

context("simulate")
test_that("simulate", {
    expect_is(simulate(gm2),"data.frame")
    expect_message(simulate(gm2,ReForm=NA),"is now preferred to")
    expect_message(simulate(gm2,REForm=NA),"is now preferred to")
    p1 <- simulate(gm2,re.form=NULL,seed=101)
    p2 <- simulate(gm2,re.form=~0,seed=101)
    p3 <- simulate(gm2,re.form=NA,seed=101)
    p4 <- simulate(gm2,ReForm=NULL,seed=101)
    p5 <- simulate(gm2,ReForm=~0,seed=101)
    p6 <- simulate(gm2,ReForm=NA,seed=101)
    p7 <- simulate(gm2,REForm=NULL,seed=101)
    p8 <- simulate(gm2,REForm=~0,seed=101)
    p9 <- simulate(gm2,REForm=NA,seed=101)
    p10 <- simulate(gm2,use.u=FALSE,seed=101)
    p11 <- simulate(gm2,use.u=TRUE,seed=101)
    expect_equal(p2,p3)
    expect_equal(p2,p5)
    expect_equal(p2,p6)
    expect_equal(p2,p8)
    expect_equal(p2,p9)
    expect_equal(p2,p10)
    expect_equal(p1,p4)
    expect_equal(p1,p7)
    expect_equal(p1,p11)
    expect_error(simulate(gm2,use.u=TRUE,re.form=NA),"should specify only one")
})
context("misc")
test_that("misc", {
    expect_equal(df.residual(fm1),176)
    expect_is(fortify(fm1),"data.frame")
    expect_is(fortify(gm1),"data.frame")
    expect_is(as.data.frame(VarCorr(fm1)),"data.frame")
})
}

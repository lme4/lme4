library("testthat")
library("lme4")

context("fitting lmer models")
test_that("lmer", {
    expect_warning(lmer(z~ 1|f, method="abc"),"Use the REML argument")
    expect_warning(lmer(z~ 1|f, method="Laplace"),"Use the REML argument")
    expect_warning(lmer(z~ 1|f, sparseX=TRUE),"has no effect at present")
    expect_that(fm1 <- lmer(Yield ~ 1|Batch, Dyestuff), is_a("lmerMod"))
    expect_that(fm1@resp,                               is_a("lmerResp"))
    expect_that(fm1@pp,                                 is_a("merPredD"))
    expect_that(fe1 <- fixef(fm1),                      is_equivalent_to(1527.5))
    expect_that(VarCorr(fm1)[[1]][1,1],                 equals(1764.07265427677))
    expect_that(isREML(fm1),                            equals(TRUE))
    expect_that(REMLfun <- as.function(fm1),            is_a("function"))
    expect_that(REMLfun(1),                             equals(319.792389042002))
    expect_that(REMLfun(0),                             equals(326.023232155879))
    expect_that(family(fm1),                            equals(gaussian()))
    expect_that(isREML(fm1ML <- refitML(fm1)),          equals(FALSE))
    ## expect_that(is.na(deviance(fm1)),                   equals(TRUE))
    expect_that(deviance(fm1ML),                        equals(327.327059881135))
    expect_that(sigma(fm1),                             equals(49.5100503990048))
    expect_that(sigma(fm1ML),                           equals(49.5100999308089))
    expect_that(extractAIC(fm1),                        equals(c(3, 333.327059881135)))
    expect_that(extractAIC(fm1ML),                      equals(c(3, 333.327059881135)))
    expect_that(vcov(fm1)[1,1],                         equals(375.720278729861))
    expect_that(vcov(fm1ML)[1,1],			equals(313.097218742665,
							       #was 313.097224695739
							       tol = 1e-7))
    expect_that(fm2 <- refit(fm1, Dyestuff2$Yield),     is_a("lmerMod"))
    expect_that(fixef(fm2),                             is_equivalent_to(5.6656))
    expect_that(VarCorr(fm2)[[1]][1,1],                 is_equivalent_to(0))
    expect_that(getME(fm2, "theta"),                    is_equivalent_to(0))
    expect_that(X  <- getME(fm1, "X"),                  is_equivalent_to(array(1, c(1, 30))))
    expect_that(Zt <- getME(fm1, "Zt"),                 is_a("dgCMatrix"))
    expect_that(dim(Zt),                                equals(c(6L, 30L)))
    expect_that(length(Zt@x),                           equals(30L))
    expect_that(Zt@x,                                   equals(rep.int(1, 30L)))
    expect_that(theta <- getME(fm1, "theta"),           is_equivalent_to(0.848330078125))
    expect_that(Lambdat <- getME(fm1, "Lambdat"),       is_a("dgCMatrix"))
    expect_that(as(Lambdat, "matrix"),                  is_equivalent_to(diag(theta, 6L, 6L)))
    expect_that(fm3 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy), is_a("lmerMod"))
    expect_that(getME(fm3,"n_rtrms"),                   equals(2L))
    expect_that(getME(fm3,"n_rfacs"),                   equals(1L))
    expect_error(fm4 <- lmer(Reaction ~ Days + (1|Subject),
                            subset(sleepstudy,Subject==levels(Subject)[1])), "must have > 1")
    expect_warning(fm4 <- lFormula(Reaction ~ Days + (1|Subject),
                             subset(sleepstudy,Subject==levels(Subject)[1]),
                             control=lmerControl(check.nlev.gtr.1="warning")), "must have > 1")
    expect_warning(fm4 <- lmer(Reaction ~ Days + (1|Subject),
                            subset(sleepstudy,Subject %in% levels(Subject)[1:4]),
                               control=lmerControl(check.nlev.gtreq.5="warning")),
                   "< 5 sampled levels")
    sstudy9 <- subset(sleepstudy, Days == 1 | Days == 9)
    expect_warning(m1 <- lmer(Reaction ~ 1 + Days + (1 + Days | Subject),
                            data = sleepstudy, subset = (Days == 1 | Days == 9)),
                 "number of observations.*rank.*unidentifiable")
    expect_warning(lFormula(Reaction ~ 1 + Days + (1 + Days | Subject),
                           data = sleepstudy, subset = (Days == 1 | Days == 9)),
                 "number of observations.*rank.*unidentifiable")
    ## test arguments: promote warning to error so that any errors will stop the test
    options(warn=2)
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, REML=TRUE), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, start=NULL), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, verbose=0L), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, subset=TRUE), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, weights=rep(1,nrow(Dyestuff))), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, na.action="na.exclude"), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, offset=rep(0,nrow(Dyestuff))), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, contrasts=NULL), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, devFunOnly=FALSE), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, control=lmerControl(optimizer="Nelder_Mead")), is_a("lmerMod"))
    expect_that(lmer(Yield ~ 1|Batch, Dyestuff, control=lmerControl()), is_a("lmerMod"))
    ## disable test ... should be no warning
    expect_is(lmer(Reaction ~ 1 + Days + (1 + Days | Subject),
                   data = sleepstudy, subset = (Days == 1 | Days == 9),
                   control=lmerControl(check.nobs.vs.rankZ="ignore")),
              "merMod")
    expect_error(lmer(Reaction ~ 1 + Days + (1|obs),
                      data = transform(sleepstudy,obs=seq(nrow(sleepstudy))),
                      "number of levels of each grouping factor"))
    expect_is(lmer(Reaction ~ 1 + Days + (1|obs),
                   data = transform(sleepstudy,obs=seq(nrow(sleepstudy))),
                   control=lmerControl(check.nobs.vs.nlev="ignore",
                   check.nobs.vs.rankZ="ignore")),
              "merMod")

    ## disable warning via options
    options(lmerControl=list(check.nlev.gtreq.5="ignore",check.nobs.vs.rankZ="ignore"))
    expect_is(fm4 <- lmer(Reaction ~ Days + (1|Subject),
                            subset(sleepstudy,Subject %in% levels(Subject)[1:4])), "merMod")
    expect_is(lmer(Reaction ~ 1 + Days + (1 + Days | Subject),
                   data = sleepstudy, subset = (Days == 1 | Days == 9)),
              "merMod")
    options(lmerControl=NULL)
    options(warn=0)
    expect_warning(lmer(Yield ~ 1|Batch, Dyestuff, junkArg=TRUE),"extra argument.*disregarded")
    expect_warning(lmer(Yield ~ 1|Batch, Dyestuff, control=list()),
                    "passing control as list is deprecated")
    expect_warning(lmer(Yield ~ 1|Batch, Dyestuff, control=glmerControl()),
                   "passing control as list is deprecated")
})


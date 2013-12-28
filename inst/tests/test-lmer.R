
stopifnot(require("testthat"), require("lme4"))

context("fitting lmer models")

test_that("lmer", {
    expect_warning(lmer(z~ 1|f, method="abc"),"Use the REML argument")
    expect_warning(lmer(z~ 1|f, method="Laplace"),"Use the REML argument")
    expect_warning(lmer(z~ 1|f, sparseX=TRUE),"has no effect at present")
    expect_is(fm1 <- lmer(Yield ~ 1|Batch, Dyestuff), "lmerMod")
    expect_is(fm1_noCD <- update(fm1,control=lmerControl(calc.derivs=FALSE)),
              "lmerMod")
    expect_equal(VarCorr(fm1),VarCorr(fm1_noCD))
    ## backward compatibility version
    expect_is(fm1.old <- update(fm1,control=lmerControl(use.last.params=TRUE)),
                                "lmerMod")
    expect_is(fm1@resp,				"lmerResp")
    expect_is(fm1@pp, 				"merPredD")
    expect_that(fe1 <- fixef(fm1),                      is_equivalent_to(1527.5))
    expect_that(VarCorr(fm1)[[1]][1,1],
                equals(1764.03751951707))
    ## back-compatibility ...
    expect_that(VarCorr(fm1.old)[[1]][1,1],
                equals(1764.07265427677))

    
    expect_that(isREML(fm1),                            equals(TRUE))
    expect_is(REMLfun <- as.function(fm1),	"function")
    expect_that(REMLfun(1),                             equals(319.792389042002))
    expect_that(REMLfun(0),                             equals(326.023232155879))
    expect_that(family(fm1),                            equals(gaussian()))
    expect_that(isREML(fm1ML <- refitML(fm1)),          equals(FALSE))
    expect_that(deviance(fm1)[["REML"]],                equals(319.654276842342))
    expect_that(deviance(fm1ML),                        equals(327.327059881135))
    expect_that(sigma(fm1),                             equals(49.5101272946856))
    expect_that(sigma(fm1.old),                         equals(49.5100503990048))
    expect_that(sigma(fm1ML),                           equals(49.5100999308089))
    expect_that(extractAIC(fm1),                        equals(c(3, 333.327059881135)))
    expect_that(extractAIC(fm1ML),                      equals(c(3, 333.327059881135)))
    expect_that(vcov(fm1)[1,1],                         equals(375.714676744045))
    expect_that(vcov(fm1.old)[1,1],                     equals(375.72027872986))
    expect_that(vcov(fm1ML)[1,1],			equals(313.09721874266, tol=1e-7))
					#		   was 313.0972246957
    expect_is(fm2 <- refit(fm1, Dyestuff2$Yield), "lmerMod")
    expect_that(fixef(fm2),                             is_equivalent_to(5.6656))
    expect_that(VarCorr(fm2)[[1]][1,1],                 is_equivalent_to(0))
    expect_that(getME(fm2, "theta"),                    is_equivalent_to(0))
    expect_that(X  <- getME(fm1, "X"),                  is_equivalent_to(array(1, c(1, 30))))
    expect_is(Zt <- getME(fm1, "Zt"),		"dgCMatrix")
    expect_that(dim(Zt),                                equals(c(6L, 30L)))
    expect_that(Zt@x,                                   equals(rep.int(1, 30L)))
    expect_that(theta <- getME(fm1, "theta"),           is_equivalent_to(0.8483203125))
    expect_that(theta.old <- getME(fm1.old, "theta"),       is_equivalent_to(0.848330078125))
    expect_is(Lambdat <- getME(fm1, "Lambdat"), "dgCMatrix")
    expect_that(as(Lambdat, "matrix"),                  is_equivalent_to(diag(theta, 6L, 6L)))
    expect_is(fm3 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy),
              					"lmerMod")
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
    expect_error(lmer(Reaction ~ 1 + Days + (1 + Days | Subject),
                        data = sleepstudy, subset = (Days == 1 | Days == 9)),
                   "no. random effects")
    expect_error(lFormula(Reaction ~ 1 + Days + (1 + Days | Subject),
                           data = sleepstudy, subset = (Days == 1 | Days == 9)),
                 "no. random effects")
    ## Promote warning to error so that warnings or errors will stop the test:
    options(warn=2)
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, REML=TRUE), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, start=NULL), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, verbose=0L), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, subset=TRUE), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, weights=rep(1,nrow(Dyestuff))), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, na.action="na.exclude"), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, offset=rep(0,nrow(Dyestuff))), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, contrasts=NULL), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, devFunOnly=FALSE), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, control=lmerControl(optimizer="Nelder_Mead")), "lmerMod")
    expect_is(lmer(Yield ~ 1|Batch, Dyestuff, control=lmerControl()), "lmerMod")
    expect_error(lmer(Yield ~ 1|Batch, Dyestuff, control=lmerControl(optimizer="optimx")),"must be loaded")
    expect_error(lmer(Yield ~ 1|Batch, Dyestuff, control=lmerControl(optimizer="junk")), "couldn't find optimizer function")
    ## disable test ... should be no warning
    expect_is(lmer(Reaction ~ 1 + Days + (1 + Days | Subject),
                   data = sleepstudy, subset = (Days == 1 | Days == 9),
                   control=lmerControl(check.nobs.vs.rankZ="ignore",
                   check.nobs.vs.nRE="ignore",
                   check.conv.hess="ignore")),
              "merMod")
    expect_error(lmer(Reaction ~ 1 + Days + (1|obs),
                      data = transform(sleepstudy,obs=seq(nrow(sleepstudy))),
                      "number of levels of each grouping factor"))
    expect_is(lmer(Reaction ~ 1 + Days + (1|obs),
                   data = transform(sleepstudy,obs=seq(nrow(sleepstudy))),
                   control=lmerControl(check.nobs.vs.nlev="ignore",
                   check.nobs.vs.nRE="ignore",
                   check.nobs.vs.rankZ="ignore")),
              "merMod")

    ## check for errors with illegal control options
    flags <- grep("^check\\.(?!conv)",names(formals(lmerControl)),
                  perl=TRUE,value=TRUE)
    for (i in seq_along(flags)) {
        ## cat(flags[i],"\n")
        expect_error(lFormula(Reaction~1+Days+(1|Subject),
                              data=sleepstudy,
                              control=do.call(lmerControl,
                              ## DELIBERATE TYPO
                              setNames(list(c("warnign")),flags[i]))),
                     "invalid control level")
    }

    ## disable warning via options
    options(lmerControl=list(check.nobs.vs.rankZ="ignore",check.nobs.vs.nRE="ignore"))
    expect_is(fm4 <- lmer(Reaction ~ Days + (1|Subject),
                            subset(sleepstudy,Subject %in% levels(Subject)[1:4])), "merMod")
    expect_is(lmer(Reaction ~ 1 + Days + (1 + Days | Subject),
                   data = sleepstudy, subset = (Days == 1 | Days == 9),
                   control=lmerControl(check.conv.hess="ignore")),
              "merMod")
    options(lmerControl=NULL)
    options(warn=0)
    expect_warning(lmer(Yield ~ 1|Batch, Dyestuff, junkArg=TRUE),"extra argument.*disregarded")
    expect_warning(lmer(Yield ~ 1|Batch, Dyestuff, control=list()),
                    "passing control as list is deprecated")
    expect_warning(lmer(Yield ~ 1|Batch, Dyestuff, control=glmerControl()),
                   "passing control as list is deprecated")

    ## test deparsing of very long terms inside mkReTrms
    set.seed(101)
    longNames <- sapply(letters[1:25],
                        function(x) paste(rep(x,8),collapse=""))
    tstdat <- data.frame(Y=rnorm(10),
                         F=factor(1:10),
                         matrix(runif(250),ncol=25,
                                dimnames=list(NULL,
                                longNames)))
    expect_is(lFormula(Y~1+(aaaaaaaa+bbbbbbbb+cccccccc+dddddddd+
                        eeeeeeee+ffffffff+gggggggg+hhhhhhhh+
                        iiiiiiii+jjjjjjjj+kkkkkkkk+llllllll|F),
                   data=tstdat,
                   control=lmerControl(check.nobs.vs.nlev="ignore",
                   check.nobs.vs.nRE="ignore",
                   check.nobs.vs.rankZ="ignore")),"list")

    ## do.call(new,...) bug
    new <- "foo"
    expect_is(refit(fm1),"merMod")
    rm("new")
})


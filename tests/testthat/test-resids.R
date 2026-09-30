#context("residuals")
test_that("lmer", {
    C1 <- lmerControl(optimizer="nloptwrap",
                      optCtrl=list(xtol_abs=1e-6, ftol_abs=1e-6))
    fm1 <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy,
                control=C1)
    fm2 <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy,
                control=lmerControl(calc.derivs=FALSE,
                                    optimizer="nloptwrap",
                                    optCtrl=list(xtol_abs=1e-6, ftol_abs=1e-6)))
    expect_equal(resid(fm1), resid(fm2))
    expect_equal(range(resid(fm1)),              c(-101.17996, 132.54664), tolerance=1e-6)
    expect_equal(range(resid(fm1, scaled=TRUE)), c(-3.9536067, 5.1792598), tolerance=1e-6)
    expect_equal(resid(fm1,"response"),resid(fm1))
    expect_equal(resid(fm1,"response"),resid(fm1,type="working"))
    expect_equal(resid(fm1,"deviance"),resid(fm1,type="pearson"))
    expect_equal(resid(fm1),resid(fm1,type="pearson"))  ## because no weights given
    expect_error(residuals(fm1,"partial"),
                 "partial residuals are not implemented yet")
    sleepstudyNA <- sleepstudy
    na_ind <- c(10,50)
    sleepstudyNA[na_ind,"Days"] <- NA
    fm1NA <- update(fm1,data=sleepstudyNA)
    fm1NA_exclude <- update(fm1,data=sleepstudyNA,na.action="na.exclude")
    expect_equal(length(resid(fm1)),length(resid(fm1NA_exclude)))
    expect_true(all(is.na(resid(fm1NA_exclude)[na_ind])))
    expect_true(!any(is.na(resid(fm1NA_exclude)[-na_ind])))
})

test_that("glmer", {
    gm1 <- glmer(incidence/size ~ period + (1|herd), cbpp,
                 family=binomial, weights=size)
    gm2 <- update(gm1,control=glmerControl(calc.derivs=FALSE))
    gm1.old <- update(gm1,control=glmerControl(calc.derivs=FALSE,
                          use.last.params=TRUE))
    expect_equal(resid(gm1),resid(gm2))
    ## Reference values use the factorization at the accepted PIRLS mean (GH #998).
    expect_equal(range(resid(gm1.old)), c(-3.19768734334282, 2.35698728407442), tolerance=1e-6)
    expect_equal(range(resid(gm1)), c(-3.19768734334282, 2.35698728407442), tolerance=1e-6)
    expect_equal(range(resid(gm1.old, "response")), c(-0.194642013418991, 0.318485961516536), tolerance=1e-6)
    expect_equal(range(resid(gm1,"response")),c(-0.194642013418991, 0.318485961516536))
    expect_equal(range(resid(gm1.old, "pearson")),  c(-2.38178751224499, 2.87959762644642),tolerance=1e-5)
    expect_equal(range(resid(gm1,"pearson")), c(-2.38178751224499, 2.87959762644642))
    expect_equal(range(resid(gm1.old, "working")),  c(-1.24168384328726, 5.4136064619901),tolerance=1e-5)
    expect_equal(range(resid(gm1, "working")),    c(-1.24168384328726, 5.4136064619901))
    expect_equal(resid(gm1),resid(gm1,scaled=TRUE))  ## since sigma==1

    expect_error(resid(gm1,"partial"),
                 "partial residuals are not implemented yet")

    cbppNA <- cbpp
    na_ind <- c(10,50)
    cbppNA[na_ind,"period"] <- NA
    gm1NA <- update(gm1,data=cbppNA)
    gm1NA_exclude <- update(gm1,data=cbppNA,na.action="na.exclude")
    expect_equal(length(resid(gm1)),length(resid(gm1NA_exclude)))
    expect_true(all(is.na(resid(gm1NA_exclude)[na_ind])))
    expect_true(!any(is.na(resid(gm1NA_exclude)[-na_ind])))
})

test_that("floating-point issues -> NaN dev resids", {
    dat <- data.frame(x=1:5, n=100, id=1:5)
    res <- glmer(cbind(x,n-x) ~ 1 + (1 | id), data=dat, family=binomial,
                 control = glmerControl(check.conv.singular = "ignore"))
    r <- residuals(res, type="deviance")
    expect_equal(unname(r[3]), 0)
})

test_that("weighted residuals", {
    skip_if(getRversion() < "4.6.0") ## methods have changed to match r-devel
    ss <- sleepstudy
    ## make sure napredict() is exercised
    ss$Reaction[1] <- NA_real_
    fm1 <- lmer(Reaction ~ 1 + (1|Subject), ss, na.action = na.exclude)
    expect_equal(head(weighted.residuals(fm1), 3),
                 structure(c(NA, -86.4020327538219, -94.3061327538219), names = c(NA, "2", "3")))

    gm1 <- glmer(round(Reaction) ~ 1 + (1|Subject), ss, na.action = na.exclude,
                 family = poisson)
    expect_equal(head(weighted.residuals(gm1), 3),
                 structure(c(NA, -4.92718365355386, -5.35397459386427), names = c(NA, "2", "3")))


})

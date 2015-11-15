library("testthat")
library("lme4")
context("residuals")
test_that("lmer", {
    fm1 <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy)
    ## fm2 <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy,
    ## control=lmerControl(calc.derivs=FALSE))
    ## all.equal(resid(fm1),resid(fm2))
    expect_equal(range(resid(fm1)), c(-101.1789,132.5466), tolerance=1e-6)
    expect_equal(range(resid(fm1, scaled=TRUE)), c(-3.953567,5.179260), tolerance=1e-6)
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
    ## y, wtres, mu change ??
    ## FIX ME:: why does turning on derivative calculation make these tests fail???
    expect_equal(range(resid(gm1.old)), c(-3.197512,2.356677), tolerance=1e-6)
    expect_equal(range(resid(gm1)), c(-3.1975034,2.35668826), tolerance=1e-6)
    expect_equal(range(resid(gm1.old, "response")), c(-0.1946736,0.3184579), tolerance=1e-6)
    expect_equal(range(resid(gm1,"response")),c(-0.194674747774946, 0.318458889275477))
    expect_equal(range(resid(gm1.old, "pearson")),  c(-2.381643,2.879069),tolerance=1e-5)
    expect_equal(range(resid(gm1,"pearson")), c(-2.38163599828335, 2.87908806084918))
    expect_equal(range(resid(gm1.old, "working")),  c(-1.241733,5.410587),tolerance=1e-5)
    expect_equal(range(resid(gm1, "working")),    c(-1.24173431447365, 5.41064465283686))
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

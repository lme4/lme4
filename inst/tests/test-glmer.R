library("testthat")
library("lme4")

context("fitting glmer models")
test_that("glmer", {
    expect_warning(glmer(z~ 1|f, family=binomial, method="abc"),"Use the nAGQ argument")
    expect_warning(glmer(z~ 1|f, family=binomial, method="Laplace"),"Use the nAGQ argument")
    expect_warning(glmer(z~ 1|f, sparseX=TRUE),"has no effect at present")
    expect_that(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             data = cbpp, family = binomial), is_a("glmerMod"))
    expect_that(gm1@resp,                               is_a("glmResp"))
    expect_that(gm1@pp,                                 is_a("merPredD"))
    expect_equal(ge1 <- unname(fixef(gm1)),             c(-1.39854982537216, -0.992335519118859,
                                                          -1.12867532780426, -1.58030423764517),
                                                         tol=5e-4)
    expect_equal(c(VarCorr(gm1)[[1]]),                  0.41245527438386, tol=6e-4)
### expect_that(family(gm1),                            equals(binomial()))
### ?? binomial() has an 'initialize' component ... and the order is different
    expect_equal(deviance(gm1),                         184.052674598026, tol=1e-5)
    expect_equal(sigma(gm1),                            1)
    expect_equal(extractAIC(gm1),                       c(5, 194.052674598026), tol=1e-5)
                
    expect_equal(theta <- unname(getME(gm1, "theta")),  0.642226809144453, tol=6e-4)
###expect_that(X  <- getME(gm1, "X"),                  is_equivalent_to(array(1, c(1, 30))))
    expect_that(Zt <- getME(gm1, "Zt"),                 is_a("dgCMatrix"))
    expect_equal(dim(Zt),                               c(15L, 56L))
    expect_equal(length(Zt@x),                          56L)
    expect_equal(Zt@x,                                  rep.int(1, 56L))
    expect_that(Lambdat <- getME(gm1, "Lambdat"),       is_a("dgCMatrix"))
    expect_equivalent(as(Lambdat, "matrix"),            diag(theta, 15L, 15L))
    expect_error(glFormula(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             data = subset(cbpp, herd==levels(herd)[1]), family = binomial),
                 "must have > 1")
    expect_warning(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                         data = subset(cbpp, herd %in% levels(herd)[1:4]),
                         family = binomial,
                         control=glmerControl(check.nlev.gtreq.5="warning")),
                   "< 5 sampled levels")
    expect_warning(fm1. <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy),
                   regexp="calling .* with family=gaussian .* as a shortcut")
    options(warn=2)
    cbppX <- transform(cbpp,prop=incidence/size)
    expect_is(glmer(prop ~ period + (1 | herd),
                             data = cbppX, family = binomial, weights=size), "glmerMod")
    expect_is(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, start=NULL),
                "glmerMod")
    expect_is(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, verbose=0L),
                "glmerMod")
    expect_is(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, subset=TRUE),
                "glmerMod")
    expect_is(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, na.action="na.exclude"),
                "glmerMod")
    expect_is(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, offset=rep(0,nrow(cbppX))),
                "glmerMod")
    expect_is(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, contrasts=NULL),
                "glmerMod")
    expect_is(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, devFunOnly=FALSE),
                "glmerMod")
    expect_is(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size,
                    control=glmerControl(optimizer="Nelder_Mead")),
                "glmerMod")
    expect_is(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, control=glmerControl()),
                 "glmerMod")
    options(warn=0)
    expect_warning(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, junkArg=TRUE),
                   "extra argument.*disregarded")
    expect_warning(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             data = cbpp, family = binomial,
                         control=list()),
                   "instead of passing a list of class")
    expect_warning(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                         data = cbpp, family = binomial,
                         control=lmerControl()),
                   "instead of passing a list of class")

})


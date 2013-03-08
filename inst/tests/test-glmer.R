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
    expect_that(ge1 <- fixef(gm1),                      is_equivalent_to(c(-1.39854982537216, -0.992335519118859,
                                                                           -1.12867532780426, -1.58030423764517)))
    expect_equal(c(VarCorr(gm1)[[1]]),                  0.41245527438386)
### expect_that(family(gm1),                            equals(binomial()))
### ?? binomial() has an 'initialize' component ... and the order is different
    expect_that(deviance(gm1),                          equals(184.052674598026))
    expect_that(sigma(gm1),                             equals(1))
    expect_that(extractAIC(gm1),                        equals(c(5, 194.052674598026)))
                
    expect_that(theta <- getME(gm1, "theta"),                    is_equivalent_to(0.642226809144453))
###expect_that(X  <- getME(gm1, "X"),                  is_equivalent_to(array(1, c(1, 30))))
    expect_that(Zt <- getME(gm1, "Zt"),                 is_a("dgCMatrix"))
    expect_that(dim(Zt),                                equals(c(15L, 56L)))
    expect_that(length(Zt@x),                           equals(56L))
    expect_that(Zt@x,                                   equals(rep.int(1, 56L)))
    expect_that(Lambdat <- getME(gm1, "Lambdat"),       is_a("dgCMatrix"))
    expect_that(as(Lambdat, "matrix"),                  is_equivalent_to(diag(theta, 15L, 15L)))
    expect_error(gm2 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             data = subset(cbpp, herd==levels(herd)[1]), family = binomial),
                 "must have > 1")
    expect_warning(gm2 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             data = subset(cbpp, herd %in% levels(herd)[1:4]), family = binomial),
                   "< 5 sampled levels")
    expect_warning(fm1. <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy),
                   regexp="calling .* with family=gaussian .* as a shortcut")
    options(warn=2)
    cbppX <- transform(cbpp,prop=incidence/size)
    expect_that(gm1 <- glmer(prop ~ period + (1 | herd),
                             data = cbppX, family = binomial, weights=size), is_a("glmerMod"))
    expect_that(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, start=NULL),
                is_a("glmerMod"))
    expect_that(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, verbose=0L),
                is_a("glmerMod"))
    expect_that(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, subset=TRUE),
                is_a("glmerMod"))
    expect_that(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, na.action="na.exclude"),
                is_a("glmerMod"))
    expect_that(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, offset=rep(0,nrow(cbppX))),
                is_a("glmerMod"))
    expect_that(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, contrasts=NULL),
                is_a("glmerMod"))
    expect_that(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, devFunOnly=FALSE),
                is_a("glmerMod"))
    expect_that(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, optimizer="Nelder_Mead"),
                is_a("glmerMod"))
    expect_that(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, control=list()),
                is_a("glmerMod"))
    options(warn=0)
    expect_warning(glmer(prop ~ period + (1 | herd),
                      data = cbppX, family = binomial, weights=size, junkArg=TRUE),
                   "extra argument.*disregarded")
})


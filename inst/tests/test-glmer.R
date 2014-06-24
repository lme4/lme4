library("testthat")
library("lme4")

testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL")))
                 as.numeric(s) else 1

gives_error_or_warning <- function (regexp = NULL, all = FALSE, ...) 
{
    function(expr) {
        res <- try(evaluate_promise(expr),silent=TRUE)
        no_error <- !inherits(res, "try-error")
        if (no_error) {
            warnings <- res$warnings

            if (!is.null(regexp) && length(warnings) > 0) {
                return(matches(regexp, all = FALSE, ...)(warnings))
            } else {
                return(expectation(length(warnings) > 0, "no warnings or errors given", 
                            paste0(length(warnings), " warnings created")))
            }
        }
        if (!is.null(regexp)) {
            return(matches(regexp, ...)(res))
        }
        else {
            expectation(TRUE, "no error thrown", "threw an error")
        }
    }
}
    ## expect_that(stop("foo"),gives_error_or_warning("foo"))
    ## expect_that(warning("foo"),gives_error_or_warning("foo"))
    ## expect_that(TRUE,gives_error_or_warning("foo"))
    ## expect_that(stop("bar"),gives_error_or_warning("foo"))
    ## expect_that(warning("bar"),gives_error_or_warning("foo"))

context("fitting glmer models")
test_that("glmer", {
    set.seed(101)
    d <- data.frame(z=rbinom(200,size=1,prob=0.5),
                    f=factor(sample(1:10,200,replace=TRUE)))
    expect_warning(glmer(z~ 1|f, d, family=binomial, method="abc"),"Use the nAGQ argument")
    expect_warning(glmer(z~ 1|f, d, family=binomial, method="Laplace"),"Use the nAGQ argument")
    expect_warning(glmer(z~ 1|f, d, sparseX=TRUE),"has no effect at present")
    expect_that(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             data = cbpp, family = binomial), is_a("glmerMod"))
    expect_that(gm1@resp,                               is_a("glmResp"))
    expect_that(gm1@pp,                                 is_a("merPredD"))
    expect_equal(ge1 <- unname(fixef(gm1)),             c(-1.39854982537216, -0.992335519118859,
                                                          -1.12867532780426, -1.58030423764517),
                                                         tolerance=5e-4)
    expect_equal(c(VarCorr(gm1)[[1]]),                  0.41245527438386, tolerance=6e-4)
### expect_that(family(gm1),                            equals(binomial()))
### ?? binomial() has an 'initialize' component ... and the order is different
    expect_equal(deviance(gm1),                         184.052674598026, tolerance=1e-5)
    expect_equal(sigma(gm1),                            1)
    expect_equal(extractAIC(gm1),                       c(5, 194.052674598026), tolerance=1e-5)

    expect_equal(theta <- unname(getME(gm1, "theta")),  0.642226809144453, tolerance=6e-4)
    expect_that(X  <- getME(gm1, "X"),                  is_equivalent_to(
        model.matrix(model.frame(~ period, data=cbpp), cbpp)))
    expect_that(Zt <- getME(gm1, "Zt"),                 is_a("dgCMatrix"))
    expect_equal(dim(Zt),                               c(15L, 56L))
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
    options(glmerControl=list(junk=1,check.conv.grad="ignore"))
    expect_warning(glmer(z~ 1|f, d, family=binomial),
                   "some options")
    options(glmerControl=NULL)
    cbppX <- transform(cbpp,prop=incidence/size)
    expect_is(glmer(prop ~ period + (1 | herd),
		      data = cbppX, family = binomial, weights=size),
	      "glmerMod")
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
if(FALSE) { ## Hadley broke this
    expect_warning(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             data = cbpp, family = binomial,
                         control=list()),
                   "instead of passing a list of class")
    expect_warning(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                         data = cbpp, family = binomial,
                         control=lmerControl()),
                   "instead of passing a list of class")
}
    ##
    load(system.file("testdata","radinger_dat.RData",package="lme4"))
    mod <- glmer(presabs~predictor+(1|species),family=binomial,
                 radinger_dat)
    expect_is(mod,"merMod")
    ## TODO: is this reliable across platforms or do we have to loosen?
    expect_equal(unname(fixef(mod)),c(0.5425528,6.4289962))
    set.seed(101)
    ## complete separation case
    d <- data.frame(y=rbinom(1000,size=1,p=0.5),
                    x=runif(1000),
                    f=factor(rep(1:20,each=50)),
                    x2=rep(0:1,c(999,1)))
    mod2 <- glmer(y~x+x2+(1|f),data=d,family=binomial,
                                   control=glmerControl(check.conv.hess="ignore",
                                                        check.conv.grad="ignore"))
    expect_equal(unname(fixef(mod2))[1:2],
                 c(-0.10036244,0.03548523), tolerance=1e-4)
    expect_true(unname(fixef(mod2)[3] < -10))
    mod3 <- update(mod2, family=binomial(link="probit"))
    # singular Hessian warning
    expect_equal(unname(fixef(mod3))[1:2], c(-0.062889, 0.022241), tolerance=1e-4)
    expect_true(fixef(mod3)[3] < -4)
    mod4 <- update(mod2, family=binomial(link="cauchit"),
                   control=glmerControl(check.conv.hess="ignore",
                                        check.conv.grad="ignore"))#--> singular Hessian warning

    ## on-the-fly creation of index variables
    if (FALSE) {
        ## FIXME: fails in testthat context -- 'd' is not found
        ##  in the parent environment of glmer() -- but works fine
        ## otherwise ...
        set.seed(101)
        d <- data.frame(y1=rpois(100,1),  x=rnorm(100), ID=1:100)
        fit1 <- glmer(y1 ~ x+(1|ID),data=d,family=poisson)
        fit2 <- update(fit1, .~ x+(1|rownames(d)))
        expect_equal(unname(unlist(VarCorr(fit1))),
                     unname(unlist(VarCorr(fit2))))
    }

    ##
    if(testLevel > 1) {
        load(system.file("testdata","mastitis.rda",package="lme4"))
        t1 <- system.time(g1 <-
                          glmer(NCM ~ birth + calvingYear + (1|sire) + (1|herd),
                                mastitis, poisson,
                                ## current (2014-04-24) default:
                                control=glmerControl(optimizer=c("bobyqa","Nelder_Mead"))))
        t2 <- system.time(g2 <- update(g1,
                         control=glmerControl(optimizer="bobyqa")))
        ## 20 (then 13.0) seconds N-M vs 8 (then 4.8) seconds bobyqa ...
        ## problem is fairly ill-conditioned so parameters
        ##  are relatively far apart even though likelihoods are OK
        expect_equal(logLik(g1),logLik(g2),tolerance=1e-7)
    }
    ## test bootstrap/refit with nAGQ>1
    gm1AGQ <- update(gm1,nAGQ=2)
    expect_equal(attr(bootMer(gm1AGQ,fixef),"bootFail"),0)

    ## do.call(new,...) bug
    new <- "foo"
    expect_that(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                             data = cbpp, family = binomial), is_a("glmerMod"))
    rm("new")

    ## test issue #47, from Wolfgang Viechtbauer
    ## create some data
    n <- 100
    ai <- rep(0:1, each = n/2)
    bi <- 1-ai
    ci <- c(rep(0,42), rep(1,8), rep(0,18), rep(1,32))
    di <- 1-ci
    event <- c(rbind(ai,ci))
    group <- rep(c(1,0), times=n)
    id    <- rep(1:n, each=2)
    gm3 <- glmer(event ~ group + (1 | id), family=binomial, nAGQ=21)
    sd3 <- sqrt(diag(vcov(gm3)))
    expect_equal(sd3, c(0.4254254, 0.424922), tolerance=1e-5)
    expect_warning(vcov(gm3,use.hessian=FALSE), "finite-difference Hessian")
    expect_equal(suppressWarnings(sqrt(diag(vcov(gm3,use.hessian=FALSE)))),
                 c(0.3840921, 0.3768747), tolerance=1e-7)
    expect_equal(sd3, unname(coef(summary(gm3))[,"Std. Error"]))
    ## test non-pos-def finite-difference Hessian ...
    if(getRversion() > "3.0.0") {
        ## saved fits are not safe with old R versions
        L <- load(system.file("testdata","polytomous_vcov_ex.RData",
                              package="lme4", mustWork=TRUE))
        expect_warning(vcov(polytomous_vcov_ex),"falling back to var-cov")
    }

    ## test convergence warnings
    L <- load(system.file("testdata","gopherdat2.RData",
                          package="lme4", mustWork=TRUE))
    g0 <- glmer(shells~prev + (1|Site)+offset(log(Area)),
                family=poisson, data=Gdat)
    ## fit year as factor: OK
    gc <- glmerControl(check.conv.grad="stop")
    expect_is(update(g0,.~.+factor(year), control=gc), "glmerMod")
    ## error/warning with year as numeric:
    ## don't have full knowledge of which platforms lead to which
    ##     results, and can't detect whether we're running on valgrind,
    ##     which changes the result on 32-bit linux ...
    ## SEGFAULT on MacOS? why?
    if (FALSE) {
       expect_that(update(g0,.~.+year),
         gives_error_or_warning("(failed to converge|pwrssUpdate did not converge)"))
    }
    ## ("(failed to converge|pwrssUpdate did not converge in)"))
    ## if (sessionInfo()$platform=="i686-pc-linux-gnu (32-bit)") {
    ##     expect_warning(update(g0, .~. +year), "failed to converge")
    ## } else {
    ##     ## MacOS x86_64-apple-darwin10.8.0 (64-bit)
    ##     ## MM's platform
    ##     ## "pwrssUpdate did not converge in (maxit) iterations"
    ##     expect_error(update(g0, .~. +year), "pwrssUpdate did not converge in")
    ## }
    ## OK if we scale & center it
    expect_is(update(g0,.~. + scale(year), control=gc), "glmerMod")
    ## not OK if we scale and don't center
    expect_warning(update(g0,.~. + scale(year,center=FALSE)),
                   "failed to converge with max|grad|")
    ## OK if center and don't scale
    expect_is(update(g0,.~. + scale(year,center=TRUE,scale=FALSE),
                     control=gc),
              "glmerMod")
})

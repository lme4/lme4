library("testthat")
library("lme4")
L <- load(system.file("testdata", "lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))

if (getRversion() > "3.0.0") {
    ## saved fits are not safe with old R versions

fm0 <- fit_sleepstudy_0
fm1 <- fit_sleepstudy_1
gm1 <- fit_cbpp_1
gm2 <- fit_cbpp_2
## More objects to use in all contexts :
set.seed(101)
dNA <- data.frame(
    xfac= sample(letters[1:10], 100, replace=TRUE),
    xcov= runif(100),
    resp= rnorm(100))
dNA[sample(1:100, 10), "xcov"] <- NA
##
rSimple <- function(rep = 2, m.u = 2, kinds = c('fun', 'boring', 'meh')) {
    stopifnot(is.numeric(rep), rep >= 1,
              is.numeric(m.u), m.u >= 1,
              is.character(kinds), (nk <- length(kinds)) >= 1)
    nobs <- rep * m.u * nk
    data.frame(kind= rep(kinds, each=rep*m.u),
               unit = gl(m.u, 1, nobs), y = round(50*rnorm(nobs)))
}
d12 <- rSimple()

data("Pixel", package="nlme")
nPix <- nrow(Pixel)
fmPix <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog), data = Pixel)


context("summary")
test_that("summary", {
    ## test for multiple-correlation-warning bug
    cc <- capture.output(summary(fit_agridat_archbold))
    expect_true(length(g <- grep("not shown by default",cc))==0 || g<=1)
})

context("anova")
test_that("lmer", {
    expect_that(anova(fm0,fm1), is_a("anova"))
    expect_warning(do.call(anova,list(fm0,fm1)), "assigning generic names")
    ##
    dat <- data.frame(y = 1:5,
                      u = c(rep("A",2), rep("B",3)),
                      t = c(rep("A",3), rep("B",2)))
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
    stats::anova(lmer(y ~ u + (1 | t), data = datfun(z), REML=FALSE),
                 lmer(y ~ 1 + (1 | t), data = datfun(z), REML=FALSE))
    ##
    ## from Roger Mundry via Roman Lustrik
    full <- lmer(resp ~ xcov + (1|xfac), data=dNA)
    null <- lmer(resp ~  1   + (1|xfac), data=dNA)
    expect_error(anova(null,full),
                 "models were not all fitted to the same size of dataset")
    ## Github issue #256  from Jonas Lindeløv -- issue is *not* specific for this dataset
    ## Two models with subset() within lmer()
    full3 <- lmer(y ~ kind + (1|unit), subset(d12, kind != 'boring'), REML=FALSE)
    null3 <- update(full3, .~. - kind)
    op <- options(warn = 2) # no warnings!
    ano3 <- anova(full3, null3)## issue #256: had warning in data != data[[1]] : ...
    o3 <- capture.output(ano3) # now prints with only one 'Data:'
    expect_equal(1, grep("^Data:", o3))
    ##
    ## no problem with subset()ting outside lmer() call:
    d12sub <- subset(d12, kind != 'boring')
    expect_is(full3s <- lmer(y ~ kind + (1|unit), d12sub, REML=FALSE), "lmerMod")
    expect_is(null3s <- update(full3s, .~. - kind), "lmerMod")
    expect_is(ano3s <- anova(full3s, null3s), "anova")
    expect_equal(ano3, ano3s, check.attributes=FALSE)
    options(op)
})

context("bootMer confint()")
set.seed(47)
test_that("bootMer", {
    CI.boot <- function(fm, nsim=10, ...)
        suppressWarnings(confint(fm, method="boot", nsim=nsim, quiet=TRUE, ...))
    ## testing bug-fix for ordering of sd/cor components in sd/cor matrix with >2 rows
    m1 <- lmer(strength~1+(cask|batch),Pastes)
    ci <- CI.boot(m1)
    corvals <- ci[grep("^cor_",rownames(ci)),]
    expect_true(all(abs(corvals) <= 1))
    ## test bootMer with GLMM, multiple RE
    ci1 <- CI.boot(gm2, nsim=3)
    ci2 <- CI.boot(gm2, nsim=5, parm=4:6)
    expect_true(nrow(ci2) == 3)
    expect_equal(ci1[4:6,], ci2, tolerance = 0.4)# 0.361
    ## bootMer with NA values
    PastesNA <- Pastes
    PastesNA$Sample[1:3] <- NA
    m2 <- update(m1, data=PastesNA)
    ci3 <- CI.boot(m2, seed=101)
    expect_equal(ci, ci3, tol = 0.06)# 0.0425
    fm1. <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    sleepstudyNA <- sleepstudy
    sleepstudyNA$Days[1:3] = NA
    fm2 <- update(fm1., data=sleepstudyNA)
    expect_true(nrow(ci4 <- CI.boot(fm2, seed=101)) == 6) # could check more
    ##
    ## semipar bootstrapping
    fm01 <- lmer(Yield ~ 1|Batch, Dyestuff)
    set.seed(1)
    suppressPackageStartupMessages(require(boot))
    boo01_sp <- bootMer(fm01, fixef, nsim = 100, use.u = TRUE,
                        type = "semiparametric")
    expect_equal(sd(boo01_sp$t), 8.215586, tol=1e-4)
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
    data(Orthodont,package="nlme")
    fmOrth <- fm <- lmer(distance ~ I(age - 11) + (I(age - 11) | Subject),
                         data = Orthodont)
    expect_equal(simulate(fm,newdata=Orthodont,seed=101),
                 simulate(fm,seed=101))
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
    expect_equal(p0[1], p1[1])
    expect_equal(p0[1], p2[1])
    expect_equal(p0[1], p3[1])
    expect_warning(predict(gm1,ReForm=NA),"is deprecated")
    ## matrix-valued predictors: Github #201 from Fabian S.
    sleepstudy$X <- cbind(1, sleepstudy$Days)
    m <- lmer(Reaction ~ -1 + X  + (Days | Subject), sleepstudy)
    pm <- predict(m, newdata=sleepstudy)
    expect_is(pm, "numeric")
    expect_equal(quantile(pm, names=FALSE),
                 c(211.006525, 260.948978, 296.87331, 328.638297, 458.155583))
    ## test spurious warning with factor as response variable
    data("Orthodont",package="MEMSS")
    silly <- glmer(Sex ~ distance + (1|Subject),
                   data=Orthodont, family=binomial)
    sillypred <- data.frame(distance=c(20, 25))
    options(warn=2) # no warnings!
    ps <- predict(silly, sillypred, re.form=NA, type="response")
    expect_is(ps, "numeric")
    expect_equal(unname(ps), c(0.999989632, 0.999997201))
    ## a case with interactions (failed in one temporary version):
    expect_warning(fmPixS <<- update(fmPix, .~. + Side), "nearly unidentifiable")
    set.seed(1); ii <- sample(nrow(Pixel), 16)
    expect_equal(predict(fmPix,  newdata=Pixel[ii,]), fitted(fmPix )[ii])
    expect_equal(predict(fmPixS, newdata=Pixel[ii,]), fitted(fmPixS)[ii])
    options(warn=0)

    set.seed(7); n <- 100; y <- rnorm(n)
    dd <- data.frame(id = factor(sample(10, n, replace = TRUE)),
                     x1 = 1, y = y, x2 = rnorm(n, mean = sign(y)))
    m <- lmer(y ~ x1 + x2 + (1 | id), data=dd)
    ##-> "fixed-effect model matrix is rank deficient so dropping 1 column / coefficient"
    summary(m)
    ii <- sample(n, 16)
    expect_equal(predict(m, newdata = dd[ii,]), fitted(m)[ii])
    ## predict(*, new..) gave Error in X %*% fixef(object) - now also drops col.
})

context("simulate")
test_that("simulate", {
    expect_is(simulate(gm2),"data.frame")
    expect_warning(simulate(gm2,ReForm=NA),"is deprecated")
    expect_warning(simulate(gm2,REForm=NA),"is deprecated")
    p1 <- simulate(gm2,re.form=NULL,seed=101)
    p2 <- simulate(gm2,re.form=~0,seed=101)
    p3 <- simulate(gm2,re.form=NA,seed=101)
    p4 <- simulate(gm2,re.form=NULL,seed=101)
    expect_warning(p5 <- simulate(gm2,ReForm=~0,seed=101),"is deprecated")
    p6 <- simulate(gm2,re.form=NA,seed=101)
    expect_warning(p7 <- simulate(gm2,REForm=NULL,seed=101),"is deprecated")
    p8 <- simulate(gm2,re.form=~0,seed=101)
    p9 <- simulate(gm2,re.form=NA,seed=101)
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
    ## hack: test with three REs
    p1 <- lmer(diameter ~ (1|plate) + (1|plate) + (1|sample), Penicillin,
               control=lmerControl(check.conv.hess="ignore"))
    expect_is(sp1 <- simulate(p1), "data.frame")
    expect_true(all(dim(sp1) == c(nrow(Penicillin), 1)))
    ## Pixel example
    expect_true(all(dim(simulate(fmPixS)) == c(nPix,1)))
    expect_true(all(dim(simulate(fmPix )) == c(nPix,1)))
})

context("misc")
test_that("misc", {
    expect_equal(df.residual(fm1),176)
    if (require(ggplot2)) {
        expect_is(fortify(fm1),"data.frame")
        expect_is(fortify(gm1),"data.frame")
    }
    expect_is(as.data.frame(VarCorr(fm1)),"data.frame")
})
}# R >= 3.0.0

context("plot")
test_that("plot", {
    ## test getData() within plot function: reported by Dieter Menne
    doFit <- function(){
        data(Orthodont,package="nlme")
        data1 <- Orthodont
        fm1 <- lmer(distance ~ age + (age|Subject), data=data1)
    }
    data(Orthodont, package="nlme")
    fm0 <- lmer(distance ~ age + (age|Subject), data=Orthodont)
    expect_is(plot(fm0),"trellis")
    suppressWarnings(rm("Orthodont"))
    fm1 <- doFit()
    pp <- plot(fm1, resid(., scaled=TRUE) ~ fitted(.) | Sex, abline = 0)
    expect_is(pp,"trellis")
})

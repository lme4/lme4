library("testthat")
library("lme4")

## use old (<=3.5.2) sample() algorithm if necessary
if ("sample.kind" %in% names(formals(RNGkind))) {
    suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
}

L <- load(system.file("testdata", "lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))

## FIXME: should test for old R versions, skip reloading test data in that
## case?
fm0 <- fit_sleepstudy_0
fm1 <- fit_sleepstudy_1
fm2 <- fit_sleepstudy_2
gm1 <- fit_cbpp_1
gm2 <- fit_cbpp_2
gm3 <- fit_cbpp_3
## More objects to use in all contexts :
set.seed(101)
dNA <- data.frame(
    xfac= sample(letters[1:10], 100, replace=TRUE),
    xcov= runif(100),
    resp= rnorm(100))
dNA[sample(1:100, 10), "xcov"] <- NA

CI.boot <- function(fm, nsim=10, seed=101, ...)
    suppressWarnings(confint(fm, method="boot", nsim=nsim,
                             quiet=TRUE, seed=seed, ...))

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
    ## test for multiple-correlation-warning bug and other 'correlation = *' variants
    ## Have 2 summary() versions, each with 3 print(.) ==> 6 x capture.output(.)
    sf.aa <- summary(fit_agridat_archbold)
    msg1 <- "Correlation.* not shown by default"
    ## message => *not* part of capture.*(.)
    expect_message(c1 <- capture.output(sf.aa), msg1)
                                        # correlation = NULL - default
    cF <- capture.output(print(sf.aa, correlation=FALSE))
    ## TODO? ensure the above gives *no* message/warning/error
    expect_identical(c1, cF)
    expect_message(
    cT <- capture.output(print(sf.aa, correlation=TRUE))
    , "Correlation.* could have been required in summary()")
    expect_identical(cF, cT[seq_along(cF)])
    sfT.aa <- summary(fit_agridat_archbold, correlation=TRUE)
    expect_message(cT2 <- capture.output(sfT.aa), msg1)
    expect_identical(cF, cT2)
    cT3 <- capture.output(print(sfT.aa, correlation=TRUE))
    expect_identical(cT, cT3)
    cF2 <- capture.output(print(sfT.aa, correlation=FALSE))
    expect_identical(cF, cF2)
})

context("anova")
test_that("lmer", {
    expect_that(suppressMessages(anova(fm0,fm1)), is_a("anova"))
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
                                REML=FALSE)), "anova")
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

    ## anova() of glmer+glm models:
    gm1 <- glmer(y~(1|u), data=dat[1:4,], family=poisson)
    gm0 <- glm(y~1,data=dat[1:4,],family=poisson)
    gm2 <- glmer(y~(1|u),data=dat[1:4,],family=poisson,nAGQ=2)
    aa <- anova(gm1,gm0)
    expect_equal(aa[2,"Chisq"],0)
    expect_error(anova(gm2,gm0),"incommensurate")

    ## anova() of lmer+glm models:
    dat2 <- dat
    set.seed(101)
    dat2$y <- rnorm(5)
    fm1 <- lmer(y~(1|u),data=dat2,REML=FALSE)
    fm0 <- lm(y~1,data=dat2)
    aa2 <- anova(fm1,fm0)
    expect_equal(aa2[2,"Chisq"],0)

})

context("bootMer confint()")
set.seed(47)
test_that("bootMer", {
    ## testing bug-fix for ordering of sd/cor components in sd/cor matrix with >2 rows
    ## FIXME: This model makes no sense [and CI.boot() fails for "nloptwrap"!]
    m1 <- lmer(strength ~ 1 + (cask|batch), Pastes,
               control = lmerControl(optimizer="bobyqa"))
    ci <- CI.boot(m1)
    corvals <- ci[grep("^cor_",rownames(ci)),]
    expect_true(all(abs(corvals) <= 1))
    ## test bootMer with GLMM, multiple RE
    expect_message(ci1 <- CI.boot(gm2, nsim=5), "singular")
    ci2 <- CI.boot(gm2, nsim=5, parm=3:6)
    ci2B <- CI.boot(gm2, nsim=5, parm="beta_")
    ## previously tested with nsim=5 vs nsim=3
    expect_true(nrow(ci2) == 4)
    expect_equal(ci2,ci2B)
    expect_equal(ci1[3:6,], ci2) ## , tolerance = 0.4)# 0.361
    ## bootMer with NA values
    PastesNA <- Pastes
    PastesNA$cask[1:3] <- NA
    ## previously set 'Sample' (sic) -- no effect!
    m2 <- update(m1, data=PastesNA)
    ci3 <- CI.boot(m2)
    expect_equal(ci, ci3, tolerance=0.2)
    sleepstudyNA <- sleepstudy
    sleepstudyNA$Days[1:3] <- NA
    m4 <- update(fm2, data = sleepstudyNA,
                 control=lmerControl(check.conv.grad = .makeCC("warning",tol=4e-3)))
    expect_true(nrow(ci4 <- CI.boot(m4)) == 6) # could check more
    ##
    ## semipar bootstrapping
    m5 <- lmer(Yield ~ 1|Batch, Dyestuff)
    set.seed(1)
    suppressPackageStartupMessages(require(boot))
    boo01.sp <- bootMer(m5, fixef, nsim = 100, use.u = TRUE,
                        type = "semiparametric")
    expect_equal(sd(boo01.sp$t), 8.215586, tolerance = 1e-4)
    ## passing FUN to confint: Torbjørn Håkan Ergon
    testFun <- function(fit) fixef(fit)[1]
    expect_equal(c(unclass(
        suppressWarnings(confint(fm1, method="boot", FUN=testFun, nsim=10,
                                 quiet=TRUE)))),
                 c(243.7551,256.9104),tolerance=1e-3)

    ## passing re.form to bootMer
    FUN <- function(.){
        predict(., type="response")
    }
    fm3 <- lmer(strength ~ (1|batch/cask), Pastes)
    expect_is(bootMer(fm3, predict, nsim=3),"boot")
    expect_is(bootMer(fm3, predict, re.form=NULL, nsim=3),"boot")
    expect_is(bootMer(fm3, predict, re.form=~(1|batch)+(1|cask:batch), nsim=3),
              "boot")
    expect_is(b3 <- bootMer(fm3, predict, re.form=~(1|batch), nsim=3),
              "boot")

    FUN_name <- function(.) getME(.,"theta")
    FUN_noname <- function(.) unname(getME(.,"theta"))

    c_name <- suppressWarnings(
        confint(fm3, method="boot", FUN=FUN_name, nsim=3, seed=101))

    c_noname <- suppressWarnings(
        confint(fm3, method="boot", FUN=FUN_noname, nsim=3, seed=101))

    expect_equal(unname(c_name),unname(c_noname))

    ## example from @Mark
    ## bootstrapping etc. on GLMMs with scale
    ## SO 37466771:
    ## https://stackoverflow.com/questions/37466771/using-profile-and-boot-method-within-confint-option-with-glmer-model
    df2 <- data.frame(
    prop1 = c(0.46, 0.471, 0.458, 0.764, 0.742, 0.746,
              0.569, 0.45,    0.491,    0.467, 0.464,
              0.556, 0.584, 0.478, 0.456, 0.46, 0.493, 0.704, 0.693, 0.651),
    prop2 = c(0.458, 0.438, 0.453, 0.731, 0.738, 0.722, 0.613,
              0.498, 0.452, 0.451, 0.447,
              0.48, 0.576, 0.484, 0.473, 0.467, 0.467, 0.722, 0.707, 0.709),
    site = factor(rep(1:8,c(2,1,3,5,3,2,1,3))))

    gaussmodel <- glmer(prop2 ~ prop1 + (1|site),
                        data=df2, family=gaussian(link="logit"))

    set.seed(101)
    bci <- suppressWarnings(confint(gaussmodel,method="boot",nsim=10,quiet=TRUE))
    expect_equal(bci,
                 structure(c(16.0861072699207, 0.0367496156026639,
                             -4.21025090053564,
                             3.1483596407467, 31.3318861915707,
                             0.0564761126030204, -1.00593089841924,
                             4.7064432268919), .Dim = c(4L, 2L),
                           .Dimnames = list(c(".sig01",
                                              ".sigma", "(Intercept)",
                                              "prop1"), c("2.5 %", "97.5 %"))),
                 tolerance=1e-5)


})

context("confint_other")
test_that("confint", {
    load(system.file("testdata", "gotway_hessianfly.rda", package = "lme4"))
    ## generated via:
    ## gotway_hessianfly_fit <- glmer(cbind(y, n-y) ~ gen + (1|block),
    ##              data=gotway.hessianfly, family=binomial,
    ##              control=glmerControl(check.nlev.gtreq.5="ignore"))
    ## gotway_hessianfly_prof <- profile(gotway_hessianfly_fit,which=1)
    ## save(list=ls(pattern="gotway"),file="gotway_hessianfly.rda")
    expect_equal(confint(gotway_hessianfly_prof)[1,1],0)
    ## FIXME: should add tests for {-1,1} bounds on correlations as well
    expect_equal(c(confint(fm1,method="Wald",parm="beta_")),
                 c(232.301892,8.891041,270.508318,12.043531),
                 tolerance=1e-5)
    ## Wald gives NA for theta values
    expect_true(all(is.na(confint(fm1,method="Wald",parm="theta_"))))

    ## check names
    ci1.p <- suppressWarnings(confint(fm1,quiet=TRUE))
    ci1.w <- confint(fm1,method="Wald")
    ci1.b <- CI.boot(fm1, nsim=2)
    expect_equal(dimnames(ci1.p),
                 list(c(".sig01", ".sigma", "(Intercept)", "Days"),
                      c("2.5 %", "97.5 %")))
    expect_equal(dimnames(ci1.p),dimnames(ci1.w))
    expect_equal(dimnames(ci1.p),dimnames(ci1.b))
    ci1.p.n <- suppressWarnings(confint(fm1,quiet=TRUE,oldNames=FALSE))
    ci1.w.n <- confint(fm1,method="Wald", oldNames=FALSE)
    ci1.b.n <- CI.boot(fm1, nsim=2, oldNames=FALSE)
    expect_equal(dimnames(ci1.p.n),
                         list(c("sd_(Intercept)|Subject", "sigma", "(Intercept)", "Days"),
                              c("2.5 %", "97.5 %")))
    expect_equal(dimnames(ci1.p.n),dimnames(ci1.w.n))
    expect_equal(dimnames(ci1.p.n),dimnames(ci1.b.n))

    ## test case of slightly wonky (spline fit fails) but monotonic profiles:
    ##
    simfun <- function(J,n_j,g00,g10,g01,g11,sig2_0,sig01,sig2_1){
        N <- sum(rep(n_j,J))
        x <- rnorm(N)
        z <- rnorm(J)
        mu <- c(0,0)
        sig <- matrix(c(sig2_0,sig01,sig01,sig2_1),ncol=2)
        u   <- MASS::mvrnorm(J,mu=mu,Sigma=sig)
        b_0j <- g00 + g01*z + u[,1]
        b_1j <- g10 + g11*z + u[,2]
        y <- rep(b_0j,each=n_j)+rep(b_1j,each=n_j)*x + rnorm(N,0,sqrt(0.5))
        sim_data <- data.frame(Y=y,X=x,Z=rep(z,each=n_j),
                               group=rep(1:J,each=n_j))
    }
    set.seed(102)
    dat <- simfun(10,5,1,.3,.3,.3,(1/18),0,(1/18))
    fit <- lmer(Y~X+Z+X:Z+(X||group),data=dat)

    if (Sys.info()["sysname"] != "SunOS" &&
        .Platform$OS.type != "windows") {
        ## doesn't produce warnings on Solaris, or win-builder ...
        expect_warning(pp <- profile(fit,"theta_"),
                       "non-monotonic profile")
        expect_warning(cc <- confint(pp),"falling back to linear interpolation")
        ## very small/unstable problem, needs large tolerance
        expect_equal(unname(cc[2,]), c(0, 0.509), tolerance=0.09) # "bobyqa" had 0.54276
    }

    badprof <- readRDS(system.file("testdata","badprof.rds",
                                   package="lme4"))
    expect_warning(cc <- confint(badprof), "falling back to linear")
    expect_equal(cc,
        array(c(0, -1, 2.50856219044636, 48.8305727797906, NA, NA,
                33.1204478717389, 1, 7.33374326592662, 68.7254711217912, -6.90462047196017, NA),
              dim = c(6L, 2L),
              dimnames = list(c(".sig01", ".sig02", ".sig03", ".sigma", "(Intercept)", "cYear"),
                              c("2.5 %", "97.5 %"))),
                 tolerance=1e-3)
})


context("refit")
test_that("refit", {
    s1 <- simulate(fm1)
    expect_is(refit(fm1,s1), "merMod")
    s2 <- simulate(fm1,2)
    expect_error(refit(fm1,s2), "refit not implemented .* lists")
    data(Orthodont,package = "nlme")
    fmOrth <- fm <- lmer(distance ~ I(age - 11) + (I(age - 11) | Subject),
                         data = Orthodont)
    expect_equal(s1 <- simulate(fm,newdata = Orthodont,seed = 101),
                 s2 <- simulate(fm,seed = 101))


    ## works *without* offset ...
    m5 <- glmer(round(Reaction) ~ Days + (1|Subject),
                data = sleepstudy, family=poisson,
                offset=rep(0,nrow(sleepstudy)))
    m5R <- refit(m5)
    ## lots of fussy details make expect_equal() on the whole object difficult
    expect_equal(coef(m5),coef(m5R),tolerance=3e-6)
    expect_equal(VarCorr(m5),VarCorr(m5R),tolerance=1e-6)
    expect_equal(logLik(m5),logLik(m5R))
})

context("predict")
test_that("predict", {
    d1 <- with(cbpp, expand.grid(period = unique(period), herd = unique(herd)))
    d2 <- data.frame(period = "1", herd = unique(cbpp$herd))
    d3 <- expand.grid(period = as.character(1:3),
                      herd = unique(cbpp$herd))
    p0 <- predict(gm1)
    p1 <- predict(gm1,d1)
    p2 <- predict(gm1,d2)
    p3 <- predict(gm1,d3)
    expect_equal(p0[1], p1[1])
    expect_equal(p0[1], p2[1])
    expect_equal(p0[1], p3[1])
    expect_warning(predict(gm1, ReForm=NA), "is deprecated")
    ## matrix-valued predictors: Github #201 from Fabian S.
    sleepstudy$X <- cbind(1, sleepstudy$Days)
    m <- lmer(Reaction ~ -1 + X  + (Days | Subject), sleepstudy)
    pm <- predict(m, newdata=sleepstudy)
    expect_is(pm, "numeric")
    expect_equal(quantile(pm, names = FALSE),
                 c(211.0108, 260.9496, 296.873, 328.6378, 458.1584), tol=1e-5)
    op <- options(warn = 2) # there should be no warnings!
    if (require("MEMSS",quietly=TRUE)) {
        ## test spurious warning with factor as response variable
        data("Orthodont", package = "MEMSS") # (differently "coded" from the 'default' "nlme" one)
        silly <- glmer(Sex ~ distance + (1|Subject),
                       data = Orthodont, family = binomial)
        sillypred <- data.frame(distance = c(20, 25))
        ps <- predict(silly, sillypred, re.form=NA, type = "response")
        expect_is(ps, "numeric")
        expect_equal(unname(ps), c(0.999989632, 0.999997201), tolerance=1e-6)
        detach("package:MEMSS")
    }
    ## a case with interactions (failed in one temporary version):
    expect_warning(fmPixS <<- update(fmPix, .~. + Side),
                   "nearly unidentifiable|unable to evaluate scaled gradient|failed to converge")
	## (1|2|3); 2 and 3 seen (as Error??) on CRAN's Windows 32bit
    options(op)

    set.seed(1); ii <- sample(nrow(Pixel), 16)
    expect_equal(predict(fmPix,  newdata = Pixel[ii,]), fitted(fmPix )[ii])
    expect_equal(predict(fmPixS, newdata = Pixel[ii,]), fitted(fmPixS)[ii])

    set.seed(7); n <- 100; y <- rnorm(n)
    dd <- data.frame(id = factor(sample(10, n, replace = TRUE)),
                     x1 = 1, y = y, x2 = rnorm(n, mean = sign(y)))
    expect_message(m <- lmer(y ~ x1 + x2 + (1 | id), data = dd),
                   "fixed-effect model matrix is rank deficient")
    expect_is(summary(m),"summary.merMod")
    ii <- sample(n, 16)
    expect_equal(predict(m, newdata = dd[ii,]), fitted(m)[ii])
    ## predict(*, new..) gave Error in X %*% fixef(object) - now also drops col.

    ## predict(*, new..) with NA in data {and non-simple model}, issue #246:
    m1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    sleepst.NA <- sleepstudy ; sleepst.NA$Days[2] <- NA
    m2 <- update(fm1, data = sleepst.NA)
    ## maybe tricky for evaluation; fm1 was defined elsewhere, so data
    expect_equal(length(predict(m2, sleepst.NA[1:4,])),4)

    ## Wrong 'b' constructed in mkNewReTrms() -- issue #257
    data(Orthodont,package="nlme")
    Orthodont <- within(Orthodont, nsex <- as.numeric(Sex == "Male"))
    m3 <- lmer(distance ~ age + (age|Subject) + (0 + Sex |Subject),
               data=Orthodont,
               control=lmerControl(check.conv.hess="ignore",
                                   check.conv.grad="ignore"))
    m4 <- lmer(distance ~ age + (age|Subject) + (0 + nsex|Subject), data=Orthodont)
    expect_equal(p3 <- predict(m3, Orthodont), fitted(m3), tolerance=1e-14)
    expect_equal(p4 <- predict(m4, Orthodont), fitted(m4), tolerance=1e-14)

    ## related to GH #275 (*passes*),
    ss <- sleepstudy
    set.seed(1)
    ss$noiseChar <- ifelse(runif(nrow(sleepstudy)) > 0.8, "Yes", "No")
    ss$noiseFactor <- factor(ss$noiseChar)
    fm4 <- lmer(Reaction ~ Days + noiseChar + (Days | Subject), ss)
    expect_equal(predict(fm4, newdata = model.frame(fm4)[2:3, ])[2],
                 predict(fm4, newdata = model.frame(fm4)[3, ]))
    fm3 <- lmer(Reaction ~ Days + noiseFactor + (Days | Subject), ss)
    expect_equal(predict(fm3, newdata = model.frame(fm3)[2:3, ])[2],
                 predict(fm3, newdata = model.frame(fm3)[3, ]))

    ## complex-basis functions in RANDOM effect
    fm5 <- lmer(Reaction~Days+(poly(Days,2)|Subject),sleepstudy)
    expect_equal(predict(fm5,sleepstudy[1,]),fitted(fm5)[1])
    ## complex-basis functions in FIXED effect
    fm6 <- lmer(Reaction~poly(Days,2)+(1|Subject),sleepstudy)
    expect_equal(predict(fm6,sleepstudy[1,]),fitted(fm6)[1])

    ## GH #414: no warning about dropping contrasts on random effects
    op <- options(warn = 2) # there should be no warnings!
    set.seed(1)
    dat <- data.frame(
        fac = rep(c("a", "b"), 100),
        grp = rep(1:25, each = 4))
    dat$y <- 0
    contr <- 0.5 * contr.sum(2)
    rownames(contr) <- c("a", "b")
    colnames(contr) <- "a"
    contrasts(dat$fac) <- contr
    m1_contr <- lmer(y~fac+(fac|grp),dat)
    pp <- predict(m1_contr,newdata=dat)
    options(op)

})

context("simulate")
test_that("simulate", {
    expect_is(simulate(gm2), "data.frame")
    expect_warning(simulate(gm2, ReForm = NA), "is deprecated")
    expect_warning(simulate(gm2, REForm = NA), "is deprecated")
    p1 <- simulate(gm2, re.form = NULL, seed = 101)
    p2 <- simulate(gm2, re.form = ~0, seed = 101)
    p3 <- simulate(gm2, re.form = NA, seed = 101)
    p4 <- simulate(gm2, re.form = NULL, seed = 101)
    expect_warning(p5 <- simulate(gm2, ReForm = ~0, seed = 101), "is deprecated")
    p6 <- simulate(gm2, re.form = NA, seed = 101)
    expect_warning(p7 <- simulate(gm2, REForm = NULL, seed = 101), "is deprecated")
    p8 <- simulate(gm2, re.form = ~0, seed = 101)
    p9 <- simulate(gm2, re.form = NA, seed = 101)
    p10 <- simulate(gm2,use.u = FALSE, seed = 101)
    p11 <- simulate(gm2,use.u = TRUE, seed = 101)
    ## minimal check of content:
    expect_identical(colSums(p1[,1]), c(incidence =  95, 747))
    expect_identical(colSums(p2[,1]), c(incidence = 109, 733))
    ## equivalences:
    ## group ~0:
    expect_equal(p2,p3)
    expect_equal(p2,p5)
    expect_equal(p2,p6)
    expect_equal(p2,p8)
    expect_equal(p2,p9)
    expect_equal(p2,p10)
    ## group 1:
    expect_equal(p1,p4)
    expect_equal(p1,p7)
    expect_equal(p1,p11)
    expect_error(simulate(gm2,use.u = TRUE, re.form = NA), "should specify only one")
    ##
    ## hack: test with three REs
    p1 <- lmer(diameter ~ (1|plate) + (1|plate) + (1|sample), Penicillin,
               control = lmerControl(check.conv.hess = "ignore",
                                     check.conv.grad = "ignore"))
    expect_is(sp1 <- simulate(p1, seed=123), "data.frame")
    expect_identical(dim(sp1), c(nrow(Penicillin), 1L))
    expect_equal(fivenum(sp1[,1]),
                 c(20.864, 22.587, 23.616, 24.756, 28.599), tolerance=0.01)
## Pixel example

    expect_identical(dim(simulate(fmPixS)), c(nPix, 1L))
    expect_identical(dim(simulate(fmPix )), c(nPix, 1L))

    ## simulation with newdata smaller/larger different from original
    fm <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin)
    expect_is(simulate(fm,newdata=Penicillin[1:10,],allow.new.levels=TRUE),"data.frame")
    expect_is(simulate(fm,newdata=do.call(rbind,replicate(4,Penicillin,simplify=FALSE))),"data.frame")

    ## negative binomial sims
    set.seed(101)
    dd <- data.frame(f=factor(rep(1:10,each=20)),
                     x=runif(200),
                     y=rnbinom(200,size=2,mu=2))
    g1 <- glmer.nb(y ~ x + (1|f), data=dd)
    th.g1 <- getME(g1, "glmer.nb.theta")
    ## changed to setting seed internally
    ts1 <- table(s1 <- simulate(g1,seed=101)[,1])
    ## ts1B <- table(s1 <- simulate(g1,seed=101)[,1])
    expect_equal(fixef(g1),
                 c("(Intercept)" = 0.630067, x = -0.0167248),
                 tolerance = 1e-4)
    ## ?? Travis is getting hung up here/ignoring tolerance spec??
    expect_equal(th.g1, 2.013, tolerance = 1e-4)
    expect_equal(th.g1, g1@call$family[["theta"]])# <- important for pkg{effects} eval(<call>)
    expect_identical(sum(s1), 413)
    expect_identical(as.vector(ts1[as.character(0:5)]),
                     ##                     c(51L, 54L, 36L, 21L, 14L, 9L))
                     c(49L,56L,32L,25L,11L,9L))

    ## de novo NB simulation ...
    s2 <- simulate(~x + (1|f),seed=101,
             family=MASS::negative.binomial(theta=th.g1),
             newparams=getME(g1,c("theta","beta")),
             newdata=dd)[,1]
    expect_equal(s1,s2)

    ## Simulate with newdata with *new* RE levels:
    d <- sleepstudy[-1] # droping the response ("Reaction")
    ## d$Subject <- factor(rep(1:18, each=10))
    ## Add 18 new subjects:
    d <- rbind(d, d)
    d$Subject <- factor(rep(1:36, each=10))
    d$simulated <- simulate(fm1, seed=1, newdata = d,
                            re.form=NULL,
                            allow.new.levels = TRUE)[,1]
    expect_equal(mean(d$simulated), 299.9384608)

    ## Simulate with weights:
    newdata <- with(cbpp, expand.grid(period=unique(period),
                                      herd=unique(herd)))
    ss <- simulate(gm1, newdata=newdata[1:3,],
                               weights=20, seed=101)[[1]]
    expect_equal(ss,
                 matrix(c(4,2,0,16,18,20),nrow=3,
                        dimnames=list(NULL,c("incidence",""))))
    ss <- simulate(gm3, newdata=newdata[1:3,],
                               weights=20, seed=101)[[1]]
    expect_equal(ss,c(0.2,0.1,0.0))
    ss <- simulate(gm1, newdata=newdata[1,],
                               weights=20, seed=101)[[1]]
    expect_equal(unname(ss),matrix(c(4,16),nrow=1))

    ## simulate Gamma, from function and de novo
    set.seed(102)
    dd <- data.frame(x=rep(seq(-2,2,length=15),10),
                     f=factor(rep(1:10,each=15)))
    u <- rnorm(10)
    dd$y <- with(dd,
                 rgamma(nrow(dd),shape=2,
                        scale=exp(2+1*x+u[as.numeric(f)])/2))
    g1 <- glmer(y~x+(1|f),family=Gamma(link="log"),dd)
    s1 <- simulate(g1,seed=101)
    s2 <- suppressMessages(simulate(~x+(1|f), family=Gamma(link="log"),
                   seed=101,
                   newdata=dd,
                   newparams=getME(g1,c("theta","beta","sigma"))))
    expect_equal(s1,s2)
    dd$y2 <- s2[[1]]
    g2 <- glmer(y2~x+(1|f),family=Gamma(link="log"),dd)
    expect_equal(fixef(g2), tolerance = 4e-7, # 32-bit windows showed 1.34e-7
                 c("(Intercept)" = 2.81887136759369, x= 1.06543222163626))
})

context("misc")
test_that("misc", {
    expect_equal(df.residual(fm1),176)
    if (require(ggplot2)) {
        expect_is(fortify(fm1), "data.frame")
        expect_is(fortify(gm1), "data.frame")
    }
    expect_is(as.data.frame(VarCorr(fm1)), "data.frame")
})

context("plot")
test_that("plot", {
    ## test getData() within plot function: reported by Dieter Menne
    doFit <- function(){
        data(Orthodont,package = "nlme")
        data1 <- Orthodont
        lmer(distance ~ age + (age|Subject), data = data1)
    }
    data(Orthodont, package = "nlme")
    fm0 <- lmer(distance ~ age + (age|Subject), data = Orthodont)
    expect_is(plot(fm0), "trellis")
    suppressWarnings(rm("Orthodont"))
    fm <- doFit()
    pp <- plot(fm, resid(., scaled = TRUE) ~ fitted(.) | Sex, abline = 0)
    expect_is(pp, "trellis")

    ## test qqmath/getIDLabels()
    expect_is(q1 <- lattice::qqmath(fm,id=0.05),"trellis")

    cake2 <- transform(cake,replicate=as.numeric(replicate),
                    recipe=as.numeric(recipe))
    fm2 <- lmer(angle ~ recipe + temp        +
                    (1|recipe:replicate), cake2, REML= FALSE)
    expect_is(lattice::qqmath(fm2,id=0.05), "trellis")
    expect_is(lattice::qqmath(fm2,id=0.05, idLabels=~recipe), "trellis")
})

context("misc")
test_that("summary", {
    ## test that family() works when $family element is weird
    ## FIXME: is convergence warning here a false positive?
    gnb <- suppressWarnings(glmer(TICKS~1+(1|BROOD),
                    family=MASS::negative.binomial(theta=2),
                    data=grouseticks))
       expect_is(family(gnb),"family")
   })


context("profile")
test_that("profile", {
    ## FIXME: can we deal with convergence warning messages here ... ?
    ## fit profile on default sd/cor scale ...
    p1 <- suppressWarnings(profile(fm1,which="theta_"))
    ## and now on var/cov scale ...
    p2 <- suppressWarnings(profile(fm1,which="theta_",
                                   prof.scale="varcov"))
    ## because there are no correlations, squaring the sd results
    ## gives the same result as profiling on the variance scale
    ## in the first place
    expect_equal(confint(p1)^2,confint(p2),
              tolerance=1e-5)
    ## or via built-in varianceProf() function
    expect_equal(unname(confint(varianceProf(p1))),
                 unname(confint(p2)),
              tolerance=1e-5)
    p3 <- profile(fm2,which=c(1,3,4))
    p4 <- suppressWarnings(profile(fm2,which="theta_",prof.scale="varcov",
                                   signames=FALSE))
    ## compare only for sd/var components, not corr component
    ## FAILS on r-patched-solaris-x86 2018-03-30 ???
    ##    2/6 mismatches (average diff: 4.62)
    ##    [1] 207 - 216 == -9.23697
    ##    [4] 1422 - 1422 == -0.00301

    if (Sys.info()["sysname"] != "SunOS") {
        expect_equal(unname(confint(p3)^2),
                     unname(confint(p4)[c(1,3,4),]),
                     tolerance=1e-3)
    }

    ## check naming convention properly adjusted
    expect_equal(as.character(unique(p4$.par)),
                 c("var_(Intercept)|Subject", "cov_Days.(Intercept)|Subject",
                   "var_Days|Subject", "sigma"))
})

context("model.frame")
test_that("model.frame", {
    ## non-syntactic names
    d <- sleepstudy
    names(d)[1] <- "Reaction Time"
    ee <- function(m,nm) {
        expect_equal(names(model.frame(m, fixed.only=TRUE)),nm)
    }
    m <- lmer(Reaction ~ 1 + (1 | Subject), sleepstudy)
    ee(m,"Reaction")
    m2 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
    ee(m2,c("Reaction","Days"))
    m3 <- lmer(`Reaction Time` ~ Days + (1 | Subject), d)
    ee(m3, c("Reaction Time","Days"))
    m4 <- lmer(Reaction ~ log(1+Days) + (1 | Subject), sleepstudy)
    ee(m4, c("Reaction","log(1 + Days)"))
})


context("influence measures")

d <- as.data.frame(ChickWeight)
colnames(d) <- c("y", "x", "subj", "tx")
dNAs <- d
dNAs$y[c(1, 3, 5)] <- NA
fitNAs <- lmer(y ~ tx*x + (x | subj), data = dNAs,
               na.action=na.exclude)

test_that("influence/hatvalues works", {
    ifm1 <- influence(fm1, do.coef=FALSE)
    expect_equal(unname(head(ifm1$hat)),
                 c(0.107483311203734, 0.102096105816528,
                   0.0980557017761242, 0.0953620990825215,
                   0.0940152977357202, 0.0940152977357202),
                 tolerance=1e-6)
    expect_equal(nrow(dNAs),length(hatvalues(fitNAs)))
})

test_that("rstudent", {
    rfm1 <- rstudent(fm1)
    expect_equal(unname(head(rfm1)),
                 c(-1.45598270922089, -1.49664543508657, -2.11747425025103,
                   -0.0729690066951975, 0.772716397142335, 2.37859408861768),
                 tolerance=1e-6)
    expect_equal(nrow(dNAs),length(rstudent(fitNAs)))
})

test_that("cooks distance", {
    expect_equal(
        unname(head(cooks.distance(fm1))),
        c(0.127645976734753, 0.127346548123793, 0.243724627125036, 0.000280638917214881,
          0.0309804642689636, 0.293554225380831),
        tolerance=1e-6)
        expect_equal(nrow(dNAs),length(cooks.distance(fitNAs)))
})

## tweaked example so estimated var = 0
zerodat <- data.frame(x=seq(0,1,length.out=120),
                      f=rep(1:3,each=40))
zerodat$y1 <- simulate(~x+(1|f),
                      family=gaussian,
                      seed=102,
                      newparams=list(beta=c(1,1),
                                     theta=c(0.001),
                                     sigma=1),
                      newdata=zerodat)[[1]]
zerodat$y2 <- simulate(~x+(1|f),
                      family=poisson,
                      seed=102,
                      newparams=list(beta=c(1,1),
                                     theta=c(0.001)),
                      newdata=zerodat)[[1]]

test_that("rstudent matches for zero-var cases",
{
    lmer_zero <- lmer(y1~x+(1|f), data=zerodat)
    glmer_zero <- glmer(y2~x+(1|f),family=poisson, data=zerodat)
    lm_zero <- lm(y1~x, data=zerodat)
    glm_zero <- glm(y2~x,family=poisson, data=zerodat)
    expect_equal(suppressWarnings(rstudent(glmer_zero)),
                 rstudent(glm_zero),
                 tolerance=0.01)
    expect_equal(suppressWarnings(rstudent(lmer_zero)),
                 rstudent(lm_zero),tolerance=0.01)
})

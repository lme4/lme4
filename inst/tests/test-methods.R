library("testthat")
library("lme4")
L <- load(system.file("testdata", "lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))

## FIXME: should test for old R versions, skip reloading test data in that
## case?
fm0 <- fit_sleepstudy_0
fm1 <- fit_sleepstudy_1
fm2 <- fit_sleepstudy_2
gm1 <- fit_cbpp_1
gm2 <- fit_cbpp_2
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
    gm1 <- glmer(y~(1|u),data=dat[1:4,],family=poisson)
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
    m1 <- lmer(strength~1+(cask|batch),Pastes)
    ci <- CI.boot(m1)
    corvals <- ci[grep("^cor_",rownames(ci)),]
    expect_true(all(abs(corvals) <= 1))
    ## test bootMer with GLMM, multiple RE
    ci1 <- CI.boot(gm2, nsim=5)
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
    expect_equal(ci, ci3, tol=0.2)
    sleepstudyNA <- sleepstudy
    sleepstudyNA$Days[1:3] <- NA
    m4 <- update(fm2, data = sleepstudyNA)
    expect_true(nrow(ci4 <- CI.boot(m4)) == 6) # could check more
    ##
    ## semipar bootstrapping
    m5 <- lmer(Yield ~ 1|Batch, Dyestuff)
    set.seed(1)
    suppressPackageStartupMessages(require(boot))
    boo01.sp <- bootMer(m5, fixef, nsim = 100, use.u = TRUE,
                        type = "semiparametric")
    expect_equal(sd(boo01.sp$t), 8.215586, tol = 1e-4)
    ## passing FUN to confint: Torbjørn Håkan Ergon
    testFun <- function(fit) fixef(fit)[1]
    expect_equal(c(unclass(
        suppressWarnings(confint(fm1, method="boot", FUN=testFun, nsim=10,
                                 quiet=TRUE)))),
                 c(243.7551,256.9104),tol=1e-3)

    ## passing re.form to bootMer
    FUN <- function(.){
        predict(., type="response")
    }
    fm2 <- lmer(strength ~ (1|batch/cask), Pastes)
    expect_is(bootMer(fm2, predict, nsim=3),"boot")
    expect_is(bootMer(fm2, predict, re.form=NULL, nsim=3),"boot")
    expect_is(bootMer(fm2, predict, re.form=~(1|batch)+(1|cask:batch), nsim=3),
              "boot")
    expect_is(b3 <- bootMer(fm2, predict, re.form=~(1|batch), nsim=3),
              "boot")
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
                 tol=1e-5)
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
    if (Sys.info()["sysname"] != "SunOS") {
        ## doesn't produce warnings on Solaris ...
        expect_warning(pp <- profile(fit,"theta_",quiet=TRUE),
                       "non-monotonic profile")
        expect_warning(cc <- confint(pp),"falling back to linear interpolation")
        ## very small/unstable problem, needs large tolerance
        expect_equal(unname(cc[2,]),c(0,0.5427609),tolerance=1e-2)
    }

    badprof <- readRDS(system.file("testdata","badprof.rds",
                                   package="lme4"))
    expect_warning(cc <- confint(badprof), "falling back to linear")
    expect_equal(cc,
        structure(c(0, -1, 2.50856219044636, 48.8305727797906, NA, NA,
                    33.1204478717389, 1, 7.33374326592662, 68.7254711217912,
                    -6.90462047196017,
                    NA), .Dim = c(6L, 2L),
                  .Dimnames = list(c(".sig01", ".sig02",
                  ".sig03", ".sigma", "(Intercept)", "cYear"),
                  c("2.5 %", "97.5 %"))),
                 tol=1e-3)
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
                 c(211.006525, 260.948978, 296.87331, 328.638297, 458.155583))
    if (require("MEMSS",quietly=TRUE)) {
        ## test spurious warning with factor as response variable
        data("Orthodont", package = "MEMSS") # (differently "coded" from the 'default' "nlme" one)
        silly <- glmer(Sex ~ distance + (1|Subject),
                       data = Orthodont, family = binomial)
    }
    sillypred <- data.frame(distance = c(20, 25))
    op <- options(warn = 2) # no warnings!
    ps <- predict(silly, sillypred, re.form=NA, type = "response")
    expect_is(ps, "numeric")
    expect_equal(unname(ps), c(0.999989632, 0.999997201))
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
               control=lmerControl(check.conv.hess="ignore"))
    m4 <- lmer(distance ~ age + (age|Subject) + (0 + nsex|Subject), data=Orthodont)
    expect_equal(p3 <- predict(m3, Orthodont), fitted(m3), tol=1e-14)
    expect_equal(p4 <- predict(m4, Orthodont), fitted(m4), tol=1e-14)

    ## related to GH #275 (*passes*),
    ss <- sleepstudy
    set.seed(1)
    ss$noiseChar <- ifelse(runif(nrow(sleepstudy)) > 0.8, "Yes", "No")
    ss$noiseFactor <- factor(ss$noiseChar)
    fm2 <- lmer(Reaction ~ Days + noiseChar + (Days | Subject), ss)
    expect_equal(predict(fm2, newdata = model.frame(fm2)[2:3, ])[2],
                 predict(fm2, newdata = model.frame(fm2)[3, ]))
    fm3 <- lmer(Reaction ~ Days + noiseFactor + (Days | Subject), ss)
    expect_equal(predict(fm3, newdata = model.frame(fm3)[2:3, ])[2],
                 predict(fm3, newdata = model.frame(fm3)[3, ]))

    ## complex-basis functions in RANDOM effect: (currently)
    fm5 <- lmer(Reaction~Days+(poly(Days,2)|Subject),sleepstudy)
    expect_equal(predict(fm5,sleepstudy[1,]),fitted(fm5)[1])
    ## complex-basis functions in FIXED effect are fine
    fm6 <- lmer(Reaction~poly(Days,2)+(1|Subject),sleepstudy)
    expect_equal(predict(fm6,sleepstudy[1,]),fitted(fm6)[1])
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
		 c(20.9412, 22.5805, 23.5575, 24.6095, 27.6997), tol=1e-5)
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
    ts1 <- table(s1 <- simulate(g1)[,1])
    expect_equal(fixef(g1), c("(Intercept)" = 0.630067, x = -0.0167248), tol = 1e-5)
    expect_equal(th.g1, 2.013, tol = 1e-4)
    expect_equal(th.g1, g1@call$family[["theta"]])# <- important for pkg{effects} eval(<call>)
    expect_identical(sum(s1), 403)
    expect_identical(as.vector(ts1[as.character(0:5)]),
                     c(51L, 54L, 36L, 21L, 14L, 9L))

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
    gnb <- glmer(TICKS~1+(1|BROOD),
                    family=MASS::negative.binomial(theta=2),
                    data=grouseticks)
       expect_is(family(gnb),"family")
   })


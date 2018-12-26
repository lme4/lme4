stopifnot(require("testthat"), require("lme4"))

set.seed(101)
n <- 500
d <- data.frame(x=rnorm(n),
                f=factor(sample(1:10,n,replace=TRUE),
                         labels=LETTERS[1:10]),
                g=factor(sample(1:25,n,replace=TRUE),
                         labels=letters[1:25]))
d$y <- suppressMessages(simulate(~1+x+(1|f)+(x|g),family=binomial,
                                 newdata=d,
                                 newparams=list(beta=c(0,1),
                                                theta=c(1,1,2,1)))[[1]])
fm1 <- glmer(y~(1|f)+(x|g),family=binomial,data=d)

context("ranef")
test_that("warn extra args", {
    expect_warning(ranef(fm1,transf=exp),"additional arguments")
})
test_that("dotplot_ranef", {
    rr <- ranef(fm1,condVar=TRUE)
    expect_is(lattice::dotplot(rr,scales=list(x = list(relation = 'free')))$g,
              "trellis")
    expect_is(lattice::dotplot(rr,transf=exp,
                     scales=list(x = list(relation = 'free')))$g,
              "trellis")
    expect_is(as.data.frame(rr),"data.frame")
    rr0 <- ranef(fm1)
    expect_is(as.data.frame(rr0),"data.frame")
})

test_that("Dyestuff consistent with lme4.0", {
    lme4.0condVarDyestuff <- c(362.3583, 362.3583, 362.3583, 362.3583, 362.3583, 362.3583)
    fm <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE)
    lme4condVarDyestuff <- drop(attr(ranef(fm,condVar=TRUE)$Batch,"postVar"))
    expect_equal(lme4.0condVarDyestuff, lme4condVarDyestuff, tolerance = 1e-3)
})


test_that("sleepstudy consistent with lme4.0", {
        lme4.0condVarsleepstudy <- matrix(c(145.71273, -21.440414,
                                            -21.44041,   5.310927), 2, 2)
        fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
        lme4condVarsleepstudy <- attr(ranef(fm,condVar=TRUE)$Subject,"postVar")[,,1]
        expect_equal(lme4.0condVarsleepstudy, lme4condVarsleepstudy, tolerance = 2e-4)
})


test_that("cbpp consistent with lme4.0", {
    lme4.0condVarcbpp <- c(0.12128867, 0.13363275, 0.08839850, 0.17337928, 0.12277914, 0.14436663,
                           0.10658333, 0.10309812, 0.21289738, 0.13740279, 0.09555677, 0.19460241,
                           0.14808316, 0.12631006, 0.15816769)
    gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                data = cbpp, family = binomial)
    lme4condVarcbpp <- as.numeric(attr(ranef(gm,condVar=TRUE)$herd,"postVar"))
    expect_equal(lme4.0condVarcbpp, lme4condVarcbpp, tolerance = 1e-3)
})

context("multiple terms per factor")
test_that("multiple terms work", {
    fm <- lmer(Reaction ~ Days + (1|Subject)+ (0+Days | Subject), sleepstudy)
    rr <- ranef(fm, condVar=TRUE)
    expect_equal(as.data.frame(rr)[c(1,19),],
                 structure(
                     list(grpvar = c("Subject", "Subject"),
                          term = structure(1:2, .Label = c("(Intercept)", "Days"), class = "factor"),
                          grp = structure(c(9L, 9L),
                                          .Label = c("309", "310", "370", "349", "350", "334", "335", "371", "308", "369",
                                                     "351", "332", "372", "333", "352", "331", "330", "337"), class = "factor"),
                          condval = c(1.5116973008, 9.32373076098),
                          condsd  = c(12.238845590, 2.33546851406)),
                     row.names = c(1L, 19L), class = "data.frame"),
                 tolerance = 1e-5)
    cv <- attr(rr$Subject, "postVar")
    expect_equal(lapply(cv, drop),
                 list(`(Intercept)` = rep(149.79166, 18),
                      Days          = rep(5.4543894, 18)), tolerance = 1e-4)
})

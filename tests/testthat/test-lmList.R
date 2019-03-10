library(lme4)
library(testthat)

## use old (<=3.5.2) sample() algorithm if necessary
if ("sample.kind" %in% names(formals(RNGkind))) {
    suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
}

context("lmList")
test_that("basic lmList", {
    set.seed(17)
    fm1. <- lmList(Reaction ~ Days | Subject, sleepstudy, pool=FALSE)
    fm1  <- lmList(Reaction ~ Days | Subject, sleepstudy)
    cf.fm1 <- data.frame(
        `(Intercept)` =
        c(244.19267, 205.05495, 203.48423, 289.68509, 285.73897, 264.25161,
          275.01911, 240.16291, 263.03469, 290.10413, 215.11177, 225.8346,
          261.14701, 276.37207, 254.96815, 210.44909, 253.63604, 267.0448),
    Days =
        c(21.764702, 2.2617855, 6.1148988, 3.0080727, 5.2660188, 9.5667679,
          9.1420455, 12.253141, -2.8810339, 19.025974, 13.493933, 19.504017,
          6.4334976, 13.566549, 11.348109, 18.056151, 9.1884448, 11.298073))
    expect_equal(signif(coef(fm1), 8), cf.fm1,
		    tolerance = 1e-7, check.attributes=FALSE)
    expect_equal(coef(fm1.), coef(fm1))
    expect_true(inherits(formula(fm1), "formula")) ## <- had been wrong till 2015-04-09
    
    sm1. <- summary(fm1.)
    sm1 <- summary(fm1)
    expect_equal(sm1$RSE, 25.5918156267, tolerance = 1e-10)
    cf1 <- confint(fm1)


    ## Calling the plot.lmList4.confint() method :
    expect_true(inherits(pcf1 <- plot(cf1), "trellis"))

})

test_that("orthodont", {
    data(Orthodont, package="nlme")
    fm2 <- lmList(distance ~ age | Subject, Orthodont)
    fe2 <- fixef(fm2)
    expect_equal(fe2, c("(Intercept)" = 16.7611111111111,
                               age = 0.660185185185185))
    expect_true(inherits(pairs(fm2), "trellis"))
})

test_that("simulated", {
    set.seed(12)
    d <- data.frame(
        g = sample(c("A","B","C","D","E"), 250, replace=TRUE),
        y1 = runif(250, max=100),
        y2 = sample(c(0,1), 250, replace=TRUE)
    )

    fm3.1 <- lmList(y1 ~ 1 | g, data=d)
    expect_equal(coef(fm3.1),    
    structure(list(`(Intercept)` = c(45.8945525606396, 50.1127995110841, 
                 49.5320538515225, 52.4286874305165, 48.7716343882989)),
              .Names = "(Intercept)",
              row.names = c("A", 
                            "B", "C", "D", "E"), class = "data.frame",
              label = "Coefficients", effectNames = "(Intercept)",
              standardized = FALSE))

    cf31 <- confint(fm3.1)
    expect_true(inherits(plot(cf31), "trellis"))

    fm3.2 <- lmList(y2 ~ 1 | g, data=d, family=binomial)
    ##                                ^^^^^^^^ "glmList"
    cf32 <- suppressMessages(confint(fm3.2,quiet=TRUE))
    expect_identical(dim(cf32), c(5L,2:1))
    expect_true(inherits(plot(cf32), "trellis"))
    expect_equal(unname(getDataPart(signif(drop(cf32), 6))),
	    cbind(c(-0.400041, -0.311489, -1.07774, -0.841075, -0.273828),
	  c( 0.743188,  0.768538, 0.0723138, 0.274392,  0.890795)))


})

test_that("cbpp", {
    ## "glmList" (2) -- here,  herd == 8 has only one observation => not estimable
    expect_warning(fm4 <- lmList(cbind(incidence, size - incidence) ~ period | herd,
                                 family=binomial, data=cbpp),
                   "Fitting failed for ")

    cf4 <- coef(fm4) # with some 5 NA's
    ## match NA locations
    expect_equal(dim(cf4),c(15,4))
    expect_identical(which(is.na(cf4)),
                     sort(as.integer(c(8+15*(0:3), 47))))

    expect_warning(fm4B <- lme4::lmList(incidence ~ period | herd,
                                        data=cbpp),
                   "Fitting failed")
    
    

if(FALSE) {
        ## FIXME: this is actually an nlme bug ...
        ## https://bugs.r-project.org/bugzilla/show_bug.cgi?id=16542
        try(summary(fm4))
        ## Error in `[<-`(`*tmp*`, use, use, ii, value = lst[[ii]]) : 
        ##   subscript out of bounds
        library(nlme)
        data("cbpp",package="lme4")
        fm6 <- nlme::lmList(incidence ~ period | herd, data=cbpp)
        try(coef(fm6))  ## coef does *not* work here
        try(summary(fm6))

        ## this is a slightly odd example because the residual df from
        ##  these fits are in fact zero ...  so pooled.SD fails, as it should

    }
})

test_that("NA,weights,offsets", {
 
    ## from GH #320
    set.seed(101)
    x <- 1:8
    y <- c(2,2,5,4,3,1,2,1)
    g <- c(1,1,1,2,2,3,3,3)
    dat <- data.frame(x=x, y=y, g=g)
    m1 <- lmList(y ~ x | g, data=dat)
    expect_false(any(is.na(coef(m1))))
    w <- runif(nrow(sleepstudy))
    m2 <- lmList(Reaction ~ Days | Subject,
                 weights=w, sleepstudy)
    ss <- subset(sleepstudy,Subject==levels(Subject)[1])
    m2X <- lm(Reaction ~ Days, ss, weights=w[1:nrow(ss)])
    expect_equal(coef(m2X),as.matrix(coef(m2))[1,])
    m3 <- lmList(Reaction ~ Days | Subject, sleepstudy)
    m4 <- lmList(Reaction ~ Days | Subject,
                 offset=w, sleepstudy)
    m4X <- lm(Reaction ~ Days, ss, offset=w[1:nrow(ss)])
    expect_equal(coef(m4X),as.matrix(coef(m4))[1,])
    expect_false(identical(m2,m3))
    expect_false(identical(m4,m3))
    m5 <- lmList(Reaction ~ Days + offset(w) | Subject, sleepstudy)
    expect_equal(coef(m5),coef(m4))

    ## more from GH 320

    dat2 <- data.frame(dat,xx=c(NA,NA,NA,1:4,NA))
    m5 <- lmList(y ~ x | g, data=dat2)
    expect_equal(unlist(coef(m5)[1,]),
                 coef(lm(y~x,subset=(g==1))))
    expect_equal(unlist(coef(m5)[3,]),
                 coef(lm(y~x,subset=(g==3))))
})

test_that("pooled",
{
    ## GH #26
    fm_lme4 <- lme4:::lmList(Reaction ~ Days | Subject, sleepstudy)
    fm_nlme <- nlme:::lmList(Reaction ~ Days | Subject, sleepstudy)
    fm_nlme_nopool <- nlme:::lmList(Reaction ~ Days | Subject, sleepstudy, pool=FALSE)
    ci_lme4_pooled <- confint(fm_lme4,pool=TRUE) #get low and high CI estimates and pooled sd
    ci_nlme_pooled <- nlme:::intervals(fm_nlme,pool=TRUE)
    expect_equal(unname(ci_lme4_pooled[,,1]),unname(ci_nlme_pooled[,c(1,3),1]))
    ci_lme4_nopool1 <- confint(fm_lme4,pool=FALSE)
    ci_lme4_nopool2 <- confint(fm_lme4)
    expect_identical(ci_lme4_nopool1,ci_lme4_nopool2)
    ## BUG in nlme::intervals ... ? can't get CIs on unpooled fits
    ## nlme::intervals(fm_nlme,pool=FALSE)
    ## nlme::intervals(fm_nlme_nopool)
    expect_equal(ci_lme4_nopool1[1:3,,1],
                 structure(c(179.433862895996, 193.026448122379, 186.785722998616, 
                             308.951475285822, 217.083442786712, 220.182727910474),
                           .Dim = c(3L, 2L), .Dimnames = list(c("308", "309", "310"),
                                                              c("2.5 %", "97.5 %"))))
})

test_that("derived variables",
          {
              fm_lme4 <- lme4:::lmList(log(Reaction) ~ Days | Subject, sleepstudy)
              fm_nlme <- nlme:::lmList(log(Reaction) ~ Days | Subject, sleepstudy)
              expect_equal(c(coef(fm_lme4)),c(coef(fm_nlme)),tolerance=1e-5)
          })
          
test_that("subset", {
    data(MathAchieve, package="nlme")
    data(MathAchSchool, package="nlme")
    RB <- merge(MathAchieve, MathAchSchool[, c("School", "Sector")],
                by="School")
    names(RB) <- tolower(names(RB))
    RB$cses <- with(RB, ses - meanses)
    cat.list.nlme <- nlme::lmList(mathach ~ cses | school,
                  subset = sector=="Catholic", 
                  data=RB)
    cat.list.lme4 <- lme4::lmList(mathach ~ cses | school,
                                  subset = sector=="Catholic", data=RB)
    expect_equal(c(coef(cat.list.lme4)),
                 c(coef(cat.list.nlme)),tolerance=1e-5)
})


library("testthat")
library("lme4")
library("lattice")

## use old (<=3.5.2) sample() algorithm if necessary
if ("sample.kind" %in% names(formals(RNGkind))) {
    suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
}

do.plots <- TRUE

L <- load(system.file("testdata/lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))

if (getRversion() > "3.0.0") {
    ## saved fits are not safe with old R versions
    gm1 <- fit_cbpp_1
    fm1 <- fit_sleepstudy_1
    fm2 <- fit_sleepstudy_2
    fm3 <- fit_penicillin_1
    fm4 <- fit_cake_1
} else {
    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)
    fm1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
    fm2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
    fm3 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
    fm4 <- lmer(angle ~ temp + recipe + (1 | replicate), data=cake)
}

context("predict")
test_that("fitted values", {
    p0 <- predict(gm1)
    p0B <- predict(gm1,newdata=cbpp)
    expect_equal(p0, p0B, tolerance=2e-5) ## ? not sure why high tolerance necessary ? works OK on Linux/R 3.5.0

    ## fitted values, unconditional (level-0)
    p1 <- predict(gm1, re.form=NA)
    expect_true(length(unique(p1))==length(unique(cbpp$period)))

    ## fitted values, random-only
    p1R <- predict(gm1, random.only=TRUE)
    expect_equal(p1+p1R,p0)

    if (do.plots) matplot(cbind(p0,p1),col=1:2,type="b")

    ## neither fixed nor random -- all zero
    expect_equal(unique(predict(gm1,re.form=NA,random.only=TRUE)),0)

})

test_that("predict with newdata", {

    newdata <- with(cbpp,expand.grid(period=unique(period),herd=unique(herd)))
    ## new data, all RE
    p2 <- predict(gm1,newdata)
    ## new data, level-0
    p3 <- predict(gm1,newdata, re.form=NA)
    p3R <- predict(gm1,newdata, random.only=TRUE)
    expect_equal(p3+p3R,p2)

    if (do.plots) matplot(cbind(p2,p3),col=1:2,type="b")

})

test_that("predict on response scale", {
    p0 <- predict(gm1)
    p5 <- predict(gm1,type="response")
    expect_equal(p5, plogis(p0))
})

test_that("predict with newdata and RE", {

    newdata <- with(cbpp,expand.grid(period=unique(period),herd=unique(herd)))
    ## explicitly specify RE
    p2 <- predict(gm1,newdata)
    p4 <- predict(gm1,newdata, re.form=~(1|herd))
    expect_equal(p2, p4)

})

test_that("effects of new RE levels", {

    newdata <- with(cbpp,expand.grid(period=unique(period),herd=unique(herd)))
    newdata2 <- rbind(newdata,
                      data.frame(period=as.character(1:4),herd=rep("new",4)))
    expect_error(predict(gm1,newdata2), "new levels detected")

    p2 <- predict(gm1,newdata)
    p6 <- predict(gm1,newdata2,allow.new.levels=TRUE)
    expect_equal(p2, p6[1:length(p2)])  ## original values should match
    ## last 4 values should match unconditional values
    expect_true(all(tail(p6,4) ==
                        predict(gm1, newdata=data.frame(period=factor(1:4)), re.form=NA)))
})

test_that("multi-group model", {

    ## fitted values
    p0 <- predict(fm3)
    expect_equal(head(round(p0,4)),
                 c(`1` = 25.9638, `2` = 22.7663, `3` = 25.7147, `4` = 23.6799,
                   `5` = 23.7629, `6` = 20.773))
    ## fitted values, unconditional (level-0)
    p1 <- predict(fm3, re.form=NA)
    expect_equal(unique(p1),22.9722222222251)

    if (do.plots) matplot(cbind(p0,p1),col=1:2,type="b")
})

test_that("multi-group model with new data", {

    newdata <- with(Penicillin,expand.grid(plate=unique(plate),sample=unique(sample)))
    ## new data, all RE
    p2 <- predict(fm3,newdata)
    ## new data, level-0
    p3 <- predict(fm3,newdata, re.form=NA)
    ## explicitly specify RE
    p4 <- predict(fm3,newdata, re.form=~(1|plate)+(~1|sample))
    p4B <- predict(fm3,newdata, re.form=~(1|sample)+(~1|plate)) ## ****
    expect_equal(p2,p4,p4B)

    p5 <- predict(fm3,newdata, re.form=~(1|sample))
    p6 <- predict(fm3,newdata, re.form=~(1|plate))
    if (do.plots) matplot(cbind(p2,p3,p5,p6),type="b",lty=1,pch=16)
})

test_that("random-slopes model", {
    p0 <- predict(fm2)
    p1 <- predict(fm2, re.form=NA)
    ## linear model, so results should be identical patterns but smaller --
    ##   not including intermediate days
    newdata <- with(sleepstudy,expand.grid(Days=range(Days),Subject=unique(Subject)))
    newdata$p2 <- predict(fm2,newdata)
    newdata$p3 <- predict(fm2,newdata, re.form=NA)
    newdata$p4 <- predict(fm2,newdata, re.form=~(0+Days|Subject))
    newdata$p5 <- predict(fm2,newdata, re.form=~(1|Subject))

    ## reference values from an apparently-working run
    refval <- structure(list(Days = c(0, 9, 0, 9, 0, 9), Subject = structure(c(1L,
                                                         1L, 2L, 2L, 3L, 3L), .Label = c("308", "309", "310", "330", "331",
                                                                              "332", "333", "334", "335", "337", "349", "350", "351", "352",
                                                                              "369", "370", "371", "372"), class = "factor"), p2 = c(253.663652396798,
                                                                                                                              430.66001930835, 211.006415533628, 227.634788908917, 212.444742696829,
                                                                                                                              257.61053840953), p3 = c(251.405104848485, 345.610678484848,
                                                                                                                                                251.405104848485, 345.610678484848, 251.405104848485, 345.610678484848
                                                                                                                                                       ), p4 = c(251.405104848485, 428.401471760037, 251.405104848485,
                                                                                                                                                          268.033478223774, 251.405104848485, 296.570900561186), p5 = c(253.663652396798,
                                                                                                                                                                                                                 347.869226033161, 211.006415533628, 305.211989169991, 212.444742696829,
                                                                                                                                                                                                                 306.650316333193)), .Names = c("Days", "Subject", "p2", "p3",
                                                                                                                                                                                                                                     "p4", "p5"), out.attrs = structure(list(dim = structure(c(2L,
                                                                                                                                                                                                                                                                             18L), .Names = c("Days", "Subject")), dimnames = structure(list(
                                                                                                                                                                                                                                                                                                                   Days = c("Days=0", "Days=9"), Subject = c("Subject=308",
                                                                                                                                                                                                                                                                                                                                                 "Subject=309", "Subject=310", "Subject=330", "Subject=331",
                                                                                                                                                                                                                                                                                                                                                 "Subject=332", "Subject=333", "Subject=334", "Subject=335",
                                                                                                                                                                                                                                                                                                                                                 "Subject=337", "Subject=349", "Subject=350", "Subject=351",
                                                                                                                                                                                                                                                                                                                                                 "Subject=352", "Subject=369", "Subject=370", "Subject=371",
                                                                                                                                                                                                                                                                                                                                                 "Subject=372")), .Names = c("Days", "Subject"))), .Names = c("dim",
                                                                                                                                                                                                                                                                                                                                                                                                   "dimnames")), row.names = c(NA, 6L), class = "data.frame")

    expect_equal(head(newdata), refval, tol=5e-7)
})

test_that("predict and plot random slopes", {

    tmpf <- function(data,...) {
        data$Reaction <- predict(fm2,data,...)
        if (do.plots) xyplot(Reaction~Days,group=Subject,data=data,type="l")
        return(unname(head(round(data$Reaction,3))))
    }
    expect_equal(tmpf(sleepstudy),c(253.664, 273.33, 292.996, 312.662, 332.329, 351.995))
    expect_equal(tmpf(sleepstudy, re.form=NA), c(251.405, 261.872, 272.34, 282.807, 293.274, 303.742))
    expect_equal(tmpf(sleepstudy, re.form= ~(0+Days|Subject)),
                 c(251.405, 271.071, 290.738, 310.404, 330.07, 349.736))
    expect_equal(tmpf(sleepstudy, re.form= ~(1|Subject)),
                 c(253.664, 264.131, 274.598, 285.066, 295.533, 306))

})

test_that("fewer random effect levels than original", {
    ## from 'Colonel Triq'
    summary(fm4)

    ## replicate 1 only appears in rows 1:18.
    ## rownames(cake[cake$replicate==1,])

    predict(fm4, newdata=cake[-1:-17,], re.form=~ (1 | replicate))
    predict(fm4, newdata=cake[-1:-18,], re.form=NA)
    predict(fm4, newdata=cake[-1:-18,], re.form=~ (1 | replicate))
    predict(fm4, newdata=cake[-1:-18,], re.form=~ (1 | replicate), allow.new.levels=TRUE)
    ##

    p0 <- predict(fm1,newdata=data.frame(Days=6,Subject=c("308","309")))
    p1 <- predict(fm1,newdata=data.frame(Days=rep(6,4),
                      Subject=c("308","309")))
    expect_equal(rep(unname(p0),2),unname(p1))
    p2 <- predict(fm1,newdata=data.frame(Days=6,Subject="308"))
    nd <- data.frame(Days=6,
                     Subject=factor("308",levels=levels(sleepstudy$Subject)))
    p3 <- predict(fm1,newdata=nd)
    expect_equal(p2,p3)
    expect_equal(p2,p0[1])
})

test_that("only drop columns when using new data", {
    ## Stack Overflow 34221564:
    ##  should only drop columns from model matrix when using *new* data
    ## NB: Fit depends on optimizer somewhat: "nloptwrap" is really better than "bobyqa"
    library(splines)
    sleep <- sleepstudy  #get the sleep data
    set.seed(1234567)
    sleep$age <- as.factor(sample(1:3,length(sleep),rep=TRUE))
    form1 <- Reaction ~ Days + ns(Days, df=4) +
        age + Days:age + (Days | Subject)
    m4 <- lmer(form1, sleep) # fixed-effect model matrix is rank deficient so dropping 1 column / coefficient
    expect_lte(REMLcrit(m4), 1713.171) # FIXME !? why this regression??  had 1700.6431;  "bobyqa" gave  1713.171
    expect_equal(unname(head(predict(m4, re.form=NA))),
                 c(255.203, 259.688, 265.71, 282.583, 294.784, 304.933),
                 tolerance = 0.008)
})

test_that("only look for columns that exist in re.form", {
    ## GH 457
    set.seed(101)
    n <- 200
    dd <- data.frame(x=1:n,
                     f=factor(rep(1:10,n/10)),
                     g=factor(rep(1:20,each=n/20)),
                     h=factor(rep(1:5,n/5)),
                     y=rnorm(n))
    m1 <- lmer(y~1 + f + (1|h/f) + (poly(x,2)|g), data=dd, control=lmerControl(calc.derivs=FALSE))
    expect_equal(unname(predict(m1,re.form= ~1 | h/f,       newdata=dd[1,])), 0.14786, tolerance=1e-4)
    expect_equal(unname(predict(m1,re.form= ~poly(x,2) | g, newdata=dd[1,])),
                 0.1533, tolerance=.001)
    ## *last* RE not included (off-by-one error)
    m1B <- lmer(y~1 + f + (1|g) + (1|h), data=dd, control=lmerControl(calc.derivs=FALSE))
    expect_equal(unname(predict(m1B,re.form=~(1|g),newdata=data.frame(f="1",g="2"))),0.1512895,tolerance=1e-5)
    set.seed(101)
    n <- 100
    xx <- c("r1", "r2", "r3", "r4", "r5")
    xxx <- c("e1", "e2", "e3")
    p <- 0.3
    School <- factor(sample(xxx, n, replace=TRUE), levels=xxx, ordered=FALSE)
    Rank <-   factor(sample(xx,  n, replace=TRUE), levels=xx,  ordered=FALSE)
    df1 <- data.frame(
        ID = as.integer(runif(n, min = 1, max = n/7)),
        xx1 = runif(n, min = 0, max = 10),
        xx2 = runif(n, min = 0, max = 10),
        xx3 = runif(n, min = 0, max = 10),
        School,
        Rank,
        yx = as.factor(rbinom(n, size = 1, prob = p))
    )
    df1 <- df1[order(df1$ID, decreasing=FALSE),]
    mm2 <- glmer(yx ~ xx1 + xx2 + xx3 + Rank +  (1 | ID) + (1 | School / Rank),
                 data = df1,
                 family = "binomial",control = glmerControl(calc.derivs =FALSE))
    n11 <-  data.frame(School= factor("e1", levels = levels(df1$School),ordered=FALSE),
                       Rank  = factor("r1", levels = levels(df1$Rank),  ordered=FALSE),
                       xx1=8.58, xx2=8.75, xx3=7.92)

    expect_equal(unname(predict(mm2, n11, type="response",re.form= ~(1 | School / Rank))),
                 0.1174628,tolerance=1e-5)

    ## bad factor levels
    mm3 <- update(mm2, . ~ . - (1|ID))
    n12 = data.frame(School="e3",Rank="r2",xx1=8.58,xx2=8.75,xx3=7.92)
    expect_equal(unname(predict(mm3, n12, type="response")),0.1832894,tolerance=1e-5)

    ## GH #452
    ## FIXME: would like to find a smaller/faster example that would test the same warning (10+ seconds)
    set.seed(101)
    n <- 300
    df2 <- data.frame(
        xx1 = runif(n, min = 0, max = 10),
        xx2 = runif(n, min = 0, max = 10),
        xx3 = runif(n, min = 0, max = 10),
        School = factor(sample(xxx, n,replace=TRUE)),
        Rank = factor(sample(xx, n, replace=TRUE)),
        yx = as.factor(rbinom(n, size = 1, prob = p))
    )
    mm4 <- suppressWarnings(glmer(yx ~ xx1 + xx2 + xx3 + Rank +  (Rank|School),
                 data = df2,
                 family = "binomial",control = glmerControl(calc.derivs =FALSE)))

    ## set tolerance to 0.1 (!) to pass win-builder on R-devel/i386 (only:
    ## tolerance = 3e-5 is OK for other combinations of (R-release, R-devel) x (i386,x64)
    expect_equal(unname(predict(mm4, n11, type="response")), 0.2675081, tolerance=0.1)
})

test_that("simulation works with non-factor", {

    set.seed(12345)
    dd <- data.frame(a=gl(10,100), b = rnorm(1000))
    test2 <- suppressMessages(simulate(~1+(b|a), newdata=dd, family=poisson,
                                       newparams= list(beta = c("(Intercept)" = 1),
                                                       theta = c(1,1,1))))
    expect_is(test2,"data.frame")
})


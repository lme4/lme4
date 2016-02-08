library(lme4)
library(testthat)

set.seed(17)
fm1. <- lmList(Reaction ~ Days | Subject, sleepstudy, pool=FALSE)
fm1  <- lmList(Reaction ~ Days | Subject, sleepstudy)
## coef(fm1)
cf.fm1 <- data.frame(
    `(Intercept)` =
        c(244.19267, 205.05495, 203.48423, 289.68509, 285.73897, 264.25161,
          275.01911, 240.16291, 263.03469, 290.10413, 215.11177, 225.8346,
          261.14701, 276.37207, 254.96815, 210.44909, 253.63604, 267.0448),
    Days =
        c(21.764702, 2.2617855, 6.1148988, 3.0080727, 5.2660188, 9.5667679,
          9.1420455, 12.253141, -2.8810339, 19.025974, 13.493933, 19.504017,
          6.4334976, 13.566549, 11.348109, 18.056151, 9.1884448, 11.298073))
stopifnot(all.equal(signif(coef(fm1), 8), cf.fm1,
		    tolerance = 1e-7, check.attributes=FALSE),
	  all.equal(coef(fm1.), coef(fm1), tolerance = 1e-15),
	  inherits(formula(fm1), "formula") # <- had been wrong till 2015-04-09
	  )

sm1. <- summary(fm1.)
sm1 <- summary(fm1)
stopifnot(all.equal(sm1$RSE, 25.5918156267, tolerance = 1e-10))
cf1 <- confint(fm1)

if(!dev.interactive(orNone=TRUE)) pdf("lmList_plots.pdf")

## Calling the plot.lmList4.confint() method :
stopifnot(inherits(pcf1 <- plot(cf1), "trellis"))
pcf1 # nice lattice plot

data(Orthodont, package="nlme")
Orthodont <- as.data.frame(Orthodont) # no "groupedData"
fm2 <- lmList(distance ~ age | Subject, Orthodont)
coef(fm2)
(fe2 <- fixef(fm2))# did fail, now fine
stopifnot(all.equal(fe2, c("(Intercept)" = 16.7611111111111,
                           age = 0.660185185185185)))
stopifnot(inherits(print(pairs(fm2)), "trellis"))
##                       ----- needs a non-trivial X matrix


set.seed(12)
d <- data.frame(
  g = sample(c("A","B","C","D","E"), 250, replace=TRUE),
  y1 = runif(250, max=100),
  y2 = sample(c(0,1), 250, replace=TRUE)
)

fm3.1 <- lmList(y1 ~ 1 | g, data=d)
coef(fm3.1)
(cf31 <- confint(fm3.1))
stopifnot(inherits(print(plot(cf31)), "trellis"))

fm3.2 <- lmList(y2 ~ 1 | g, data=d, family=binomial)
(cf32 <- confint(fm3.2)) #                 ^^^^^^^^ "glmList"
stopifnot(identical(dim(cf32), c(5L,2:1)),
	  inherits(print(plot(cf32)), "trellis"),
	  all.equal(unname(getDataPart(signif(drop(cf32), 6))),
		    cbind(c(-0.400041, -0.311489, -1.07774, -0.841075, -0.273828),
			  c( 0.743188,  0.768538, 0.0723138, 0.274392,  0.890795))))


## "glmList" (2) -- here,  herd == 8 has only one observation => not estimable
expect_warning(fm4 <- lmList(cbind(incidence, size - incidence) ~ period | herd,
             family=binomial, data=cbpp),
             "Fitting failed for ")

fm4 # no pooled SD for glm
(cf4 <- coef(fm4)) # with some 5 NA's
## match NA locations
stopifnot(dim(cf4) == c(15,4),
          identical(which(is.na(cf4)),
                    sort(as.integer(c(8+15*(0:3), 47)))))

fm5 <- lmList(incidence ~ period | herd, data=cbpp)
fm6 <- nlme::lmList(incidence ~ period | herd, data=cbpp)

## for this example coef() *does* work ...
ctab <- t(sapply(split(cbpp,cbpp$herd),
                 function(x) {
    if (nrow(x)==1) {
        rep(NA,4)
    } else {
        g <- glm(cbind(incidence, size-incidence) ~ period, data=x,
                 family=binomial)
        cc <- coef(g)
        length(cc) <- 4  ## pad with NAs
        cc
    }}))
stopifnot(all.equal(c(ctab),c(as.matrix(coef(fm4)))))

if(FALSE) {## FIXME: but I (BMB) think this is actually an nlme bug ...
    summary(fm4)

    library(nlme)
    data("cbpp",package="lme4")
    fm6 <- nlme::lmList(incidence ~ period | herd, data=cbpp)
    ## Warning message:
    ## In lmList.formula(incidence ~ period | herd, data = cbpp) :
    ##   An lm fit failed, probably because a factor only had one level
    try(coef(fm6))  ## coef does *not* work here
    try(summary(fm6))
}

## this is a slightly odd example because the residual df from
##  these fits are in fact zero ...  so pooled.SD fails, as it should

## library(reshape2)
## library(ggplot2)
## ggplot(melt(as.matrix(coef(fm3))),
##       aes(value,Var2,colour=factor(Var1)))+
##    geom_point()+
##       geom_path(aes(group=factor(Var1)))

if(getRversion() < "3.2.0") {
    if(interactive()) break # gives an error
    else q() # <- undesirable when interactive !
}


## Try all "standard" (statistical) S3 methods:
.S3generics <- function(class) {
    s3m <- .S3methods(class=class)
    ii <- attr(s3m, "info")
    ii[!ii[, "isS4"], "generic"]
}

s3fn <- .S3generics(class= class(fm3.1)[1]) ## works for "old and new" class

noquote(s3fn <- s3fn[s3fn != "print"])# <-- it is show() not print() that works
## [1] coef    confint fitted fixef     formula logLik  pairs  plot
## [9] predict qqnorm  ranef  residuals sigma   summary update

## In lme4 1.1-7 (July 2014), only these worked:
##  coef(), confint(), formula(), logLik(), summary(), update()

## pairs() is excluded for fm3.1 which has only intercept:
## no errors otherwise:
(evs <- sapply(s3fn[s3fn != "pairs"], do.call, args = list(fm3.1)))

cls <- sapply(evs, function(.) class(.)[1])
clsOk <- cls[c("confint", "fixef", "formula", "logLik",
               "ranef", "sigma", "summary", "update")]
stopifnot(identical(unname(clsOk),
                    c("lmList4.confint", "numeric", "formula", "logLik",
                      "ranef.lmList", "numeric", "summary.lmList", "lmList4")))

## --- fm2 --- non-trivial X: can use pairs(), too:
(evs2 <- sapply(s3fn, do.call, args = list(fm2)))

## --- fm3.2 --- no failures for this "glmList" :
evs3.2 <- sapply(s3fn[s3fn != "pairs"], do.call, args = list(fm3.2))
evs3.2

## --- fm4 ---
evs4 <- sapply(s3fn, function(fn)
    tryCatch(do.call(fn, list(fm4)), error=function(e) e))
length(warnings())
unique(warnings()) ##  glm.fit: fitted probabilities numerically 0 or 1 occurred

str(sapply(evs4, class)) # more errors than above
isok4 <- !sapply(evs4, is, class2="error")

## includes a nice pairs():
evs4[isok4]
## Error msgs of those with errors, first 5, now 3, then 2 :
str(errs4 <- lapply(evs4[!isok4], conditionMessage))
## $ logLik : chr "log-likelihood not available with NULL fits"
## $ summary: chr "subscript out of bounds"
stopifnot(length(errs4) <= 2)

## from GH #320
x <- c(1,2,3,4,5,6,7,8)
y <- c(2,2,5,4,3,1,2,1)
g <- c(1,1,1,2,2,3,3,3)
dat <- data.frame(x=x, y=y, g=g)
m1 <- lmList(y ~ x | g, data=dat)
stopifnot(!any(is.na(coef(m1))))
m2 <- lmList(Reaction ~ Days | Subject,
             weights=runif(nrow(sleepstudy)), sleepstudy)
m3 <- lmList(Reaction ~ Days | Subject, sleepstudy)
m4 <- lmList(Reaction ~ Days | Subject,
             offset=runif(nrow(sleepstudy)), sleepstudy)
stopifnot(!identical(m2,m3))
stopifnot(!identical(m4,m3))

## more from GH 320

dat2 <- data.frame(dat,xx=c(NA,NA,NA,1:4,NA))
m5 <- lmList(y ~ x | g, data=dat2)
stopifnot(all.equal(unlist(coef(m5)[1,]),
          coef(lm(y~x,subset=(g==1)))))
stopifnot(all.equal(unlist(coef(m5)[3,]),
          coef(lm(y~x,subset=(g==3)))))

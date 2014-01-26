## test for models containing data-defined bases

## ?makepredictcall
## ?model.frame
## ????

data(sleepstudy,package="lme4")
library(splines)

## lm0 <- lm(Reaction~ns(Days,2),sleepstudy)
## attr(terms(lm0),"predvars")
## library(nlme)
## lme1 <- lme(Reaction~ns(Days,2),random=~1|Subject,sleepstudy)
## attr(terms(lme1),"predvars")  ## no!
## attr(lme1$terms,"predvars")   ## yes
## detach("package:nlme")

library(lme4)
library(testthat)
fm1 <- lmer(Reaction ~ ns(Days,2) + (1|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ poly(Days,2) + (1|Subject), sleepstudy)
fm3 <- lmer(Reaction ~ poly(Days,2,raw=TRUE) + (1|Subject), sleepstudy)

newdat0 <- data.frame(Days=unique(sleepstudy$Days))
newdat <- data.frame(Days=5:12)
tmpf <- function(fit) {
    with(sleepstudy,plot(Reaction~Days,xlim=c(0,12)))
    with(sleepstudy,points(Days,predict(fit),col=2))
    with(newdat0,lines(Days,predict(fit,ReForm=NA,newdata=newdat0),col=4))
    with(newdat,lines(Days,predict(fit,ReForm=NA,newdata=newdat),col=5))
}

stopifnot(all.equal(predict(fm2,newdat,ReForm=NA),
                    predict(fm3,newdat,ReForm=NA)))

## pictures
tmpf(fm1)
tmpf(fm2)
tmpf(fm3)

## test for GLMMs
set.seed(101)
d <- data.frame(y=rbinom(10,size=1,prob=0.5),
                x=1:10,
                f=factor(rep(1:5,each=2)))
gm1 <- glmer(y ~ poly(x,2) + (1|f), d, family=binomial)
gm2 <- glmer(y ~ poly(x,2,raw=TRUE) + (1|f), d, family=binomial)

newdat <- data.frame(x=c(1,4,6))
stopifnot(all.equal(predict(gm1,newdat,ReForm=NA),
                    predict(gm2,newdat,ReForm=NA),tolerance=3e-6))


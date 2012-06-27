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
    with(newdat0,lines(Days,predict(fit,REform=NA,newdata=newdat0),col=4))
    with(newdat,lines(Days,predict(fit,REform=NA,newdata=newdat),col=5))
}

## pictures
tmpf(fm1)
tmpf(fm2)
tmpf(fm3)
stopifnot(all.equal(predict(fm2,newdat,REform=NA),
                    predict(fm3,newdat,REform=NA)))

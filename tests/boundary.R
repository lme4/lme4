## Stephane Laurent:
dat <- read.csv("dat20101314.csv")
library(lme4)
fit <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat)
fit_b <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              optimizer="bobyqa")
fit_c <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=list(restart=FALSE))
getME(fit_c,"theta") ## some are zero

if (FALSE) {
    ## FIXME: fails on r-forge test on Windows x64 ... ???
    stopifnot(all.equal(getME(fit,"theta"),getME(fit_b,"theta"),tol=1e-6))
}


## Manuel Koller
source("koller-data.R")

ldata <- getData(13)
fm4 <- lmer(y ~ (1|Var2), ldata)
getME(fm4,"theta")
fm4b <- lmer(y ~ (1|Var2), ldata, control=list(restart=FALSE))
getME(fm4b,"theta")
fm4c <- lmer(y ~ (1|Var2), ldata, optimizer="bobyqa")
stopifnot(all.equal(getME(fm4,"theta"),getME(fm4c,"theta"),tol=1e-4))

## dd <- lmer(y ~ (1|Var2), ldata, devFunOnly=TRUE)
## tvec <- 10^seq(-7,0,by=0.1)
## dvec <- sapply(tvec,dd)
## d0 <- dd(0)
## plot(tvec,dvec,type="b")
## plot(tvec,abs(dvec-d0),log="xy",col=ifelse(dvec<d0,"black","red"))


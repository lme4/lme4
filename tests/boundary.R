## Stephane Laurent:
dat <- read.csv(system.file("testdata","dat20101314.csv",package="lme4"))
library(lme4)

## FIXME: in lmeControl branch, bobyqa fails: is this also the case
## in master branch?  (i.e. did I mess up something the way the controls are passed?)

## boundary fit is correct
fit <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat)
## bobyqa fails with or without restart_edge here: nelder_mead succeeds either way
## FIXME: can we tweak bobyqa controls to get it to work?
fit_b <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=lmerControl(optimizer="bobyqa",restart_edge=FALSE))
fit_c <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=lmerControl(restart_edge=FALSE))
stopifnot(all(getME(fit_c,"theta")[1:2]==0)) ## some are zero

## Manuel Koller
source(system.file("testdata","koller-data.R",package="lme4"))

ldata <- getData(13)
fm4 <- lmer(y ~ (1|Var2), ldata)
stopifnot(getME(fm4,"theta")==0)
fm4b <- lmer(y ~ (1|Var2), ldata, control=lmerControl(restart=FALSE))
stopifnot(getME(fm4b,"theta")==0)
fm4c <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="bobyqa"))
stopifnot(all.equal(getME(fm4,"theta"),getME(fm4c,"theta"),tol=1e-4))

## dd <- lmer(y ~ (1|Var2), ldata, devFunOnly=TRUE)
## tvec <- 10^seq(-7,0,by=0.1)
## dvec <- sapply(tvec,dd)
## d0 <- dd(0)
## plot(tvec,dvec,type="b")
## plot(tvec,abs(dvec-d0),log="xy",col=ifelse(dvec<d0,"black","red"))


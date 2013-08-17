## In both of these cases boundary fit (i.e. estimate of zero RE
## variance) is *incorrect*. (Nelder_Mead, restart_edge=FALSE) is the
## only case where we get stuck; either optimizer=bobyqa or
## restart_edge=TRUE (default) works

## Stephane Laurent:
dat <- read.csv(system.file("testdata","dat20101314.csv",package="lme4"))
library(lme4)

fit <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat)
fit_b <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=lmerControl(optimizer="bobyqa",restart_edge=FALSE))
fit_c <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=lmerControl(restart_edge=FALSE))
## tol=1e-5 seems OK in interactive use but not in R CMD check ... ??
stopifnot(all.equal(getME(fit,"theta"),getME(fit_b,"theta"),tol=2e-5))
stopifnot(all(getME(fit,"theta")>0))

## Manuel Koller

source(system.file("testdata","koller-data.R",package="lme4"))

ldata <- getData(13)
fm4 <- lmer(y ~ (1|Var2), ldata)
fm4b <- lmer(y ~ (1|Var2), ldata, control=lmerControl(restart_edge=FALSE))
stopifnot(getME(fm4b,"theta")==0)
fm4c <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="bobyqa"))
stopifnot(all.equal(getME(fm4,"theta"),getME(fm4c,"theta"),tol=1e-4))
stopifnot(all(getME(fm4,"theta")>0))

## dd <- lmer(y ~ (1|Var2), ldata, devFunOnly=TRUE)
## tvec <- 10^seq(-7,0,by=0.1)
## dvec <- sapply(tvec,dd)
## d0 <- dd(0)
## plot(tvec,dvec,type="b")
## plot(tvec,abs(dvec-d0),log="xy",col=ifelse(dvec<d0,"black","red"))


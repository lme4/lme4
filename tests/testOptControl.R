## https://github.com/lme4/lme4/issues/59
library(lme4)
dat <- read.csv(system.file("testdata","dat20101314.csv",package="lme4"))
library(lme4)

NMcopy <- lme4:::Nelder_Mead

fit <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
        control=lmerControl("NMcopy", check.numlev.gtreq.5="ignore"))

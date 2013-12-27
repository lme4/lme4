library(lme4)
n <- nrow(sleepstudy)
v <- rpois(n,1) + 1
w <- 1/v

fm1    <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE, weights=w)
fm2    <- lmer(Reaction ~ Days + (1    | Subject), sleepstudy, REML=FALSE, weights=w)

fm1.10 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE, weights=w*10)
fm2.10 <- lmer(Reaction ~ Days + (1    | Subject), sleepstudy, REML=FALSE, weights=w*10)

stopifnot(all.equal(anova(fm1, fm2)$Chisq[2],
                    anova(fm1.10, fm2.10)$Chisq[2]))

stopifnot(all.equal(AIC(fm1)-AIC(fm2),
                    AIC(fm1.10)-AIC(fm2.10)))

stopifnot(all.equal(BIC(fm1)-BIC(fm2),
                    BIC(fm1.10)-BIC(fm2.10)))

stopifnot(all.equal(logLik(fm1)-logLik(fm2),
                    logLik(fm1.10)-logLik(fm2.10)))

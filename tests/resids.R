library(lme4)

## raw residuals for LMMs
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
stopifnot(all.equal(residuals(fm1),sleepstudy$Reaction-fitted(fm1)))

r1 <- residuals(fm1,type="pearson")

## deviance/Pearson residuals for GLMMs
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             family = binomial, data = cbpp)
p <- fitted(gm1)
n <- cbpp$size
v <- n*p*(1-p)
obs_p <- cbpp$incidence/cbpp$size
rp <- residuals(gm1,"pearson")
rp1 <- (obs_p-p)/sqrt(p*(1-p))
rp2 <- rp1*n
## FIXME:: restore this test
## stopifnot(all.equal(rp,rp2))

r2 <- residuals(gm1,type="deviance")

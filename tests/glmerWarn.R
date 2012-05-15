library(lme4)
library(testthat)
m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = binomial, data = cbpp)
expect_that(m2 <- lmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = binomial, data = cbpp),
            gives_warning())

## FIXME: should glmer(..., [no family]) give a warning as well?
stopifnot(all.equal(m1,m2))

m3 <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy)
m4 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
m5 <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy, family=gaussian)
all.equal(fixef(m3),fixef(m5))
all.equal(m3,m4)
isTRUE(all.equal(m3,m5))

## would like m3==m5 != m4 ??
VarCorr(m4)
VarCorr(m5)  ## wrong??? is this the report or the
getME(m4,"theta")
getME(m5,"theta")

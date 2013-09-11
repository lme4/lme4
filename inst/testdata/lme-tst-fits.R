#### Saved fits for lme4 testing
####  ----------------------------------
if(FALSE) ### "Load" these by
    load(system.file("testdata/lme-tst-fits.rda", package="lme4", mustWork=TRUE))

library(lme4)
## intercept only in both fixed and random effects
fit_sleepstudy_0 <- lmer(Reaction ~ 1 + (1|Subject), sleepstudy)
## fixed slope, intercept-only RE
fit_sleepstudy_1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
## fixed slope, intercept & slope RE
fit_sleepstudy_2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
## fixed slope, independent intercept & slope RE
fit_sleepstudy_3 <- lmer(Reaction ~ Days + (1|Subject)+
                         (0+Days|Subject), sleepstudy)

save(list=ls(pattern="fit_sleepstudy"),file="lme-tst-fits.rda")

load("lme-tst-fits.rda")
cbpp$obs <- factor(seq(nrow(cbpp)))
fit_cbpp_0 <- glmer(cbind(incidence,size-incidence) ~ 1 + (1|herd),
                    cbpp, family=binomial)
fit_cbpp_1 <- update(fit_cbpp_0, . ~ . + period)
fit_cbpp_2 <- update(fit_cbpp_1, . ~ . + (1|obs))

save(list=ls(pattern="fit_"),file="lme-tst-fits.rda")

#### Saved fits for lme4 testing
####  ----------------------------------
fn <- system.file("testdata/lme-tst-fits.rda", package="lme4", mustWork=TRUE)
if(FALSE) ### "Load" these by
    load(fn)

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

save(list=ls(pattern="fit_"),file="lme-tst-fits.rda")

load("lme-tst-fits.rda")
cbpp$obs <- factor(seq(nrow(cbpp)))
## intercept-only fixed effect
fit_cbpp_0 <- glmer(cbind(incidence,size-incidence) ~ 1 + (1|herd),
                    cbpp, family=binomial)
## include fixed effect of period
fit_cbpp_1 <- update(fit_cbpp_0, . ~ . + period)
## include observation-level RE
fit_cbpp_2 <- update(fit_cbpp_1, . ~ . + (1|obs))
## specify formula by proportion/weights instead
fit_cbpp_3 <- update(fit_cbpp_0,
                     formula=incidence/size ~ period + (1 | herd),
                     weights=size)
save(list=ls(pattern="fit_"),file="lme-tst-fits.rda")

## an example with >20 fixed effects (for testing print.summary.merMod)
if (require(agridat)) {
    dat <- archbold.apple
    ## Define main plot and subplot
    dat <- transform(dat, rep=factor(rep),
                     spacing=factor(spacing), trt=factor(trt),
                     mp = factor(paste(row,spacing,sep="")),
                     sp = factor(paste(row,spacing,stock,sep="")))
    fit_agridat_archbold <- lmer(yield ~ -1 + trt + (1|rep/mp/sp), dat)
    save(list=ls(pattern="fit_"),file="lme-tst-fits.rda")
}

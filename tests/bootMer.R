library(lme4)
library(testthat)

mySumm <- function(.) {
  s <- sigma(.)
  c(beta =getME(., "beta"),
    sigma = s,
    sig01 = unname(s * getME(., "theta")))
}

fm1 <- lmer(Yield ~ 1|Batch, Dyestuff)
boo01 <- bootMer(fm1, mySumm, nsim = 10)
boo02 <- bootMer(fm1, mySumm, nsim = 10, use.u = TRUE)

## boo02 <- bootMer(fm1, mySumm, nsim = 500, use.u = TRUE)
if (require(boot)) {
    boot.ci(boo02,index=2,type="perc")
}

fm2 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
boo03 <- bootMer(fm2, mySumm, nsim = 10)
boo04 <- bootMer(fm2, mySumm, nsim = 10, use.u = TRUE)

if (lme4:::testLevel() > 1) {
    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)
    boo05 <- bootMer(gm1, mySumm, nsim = 10)
    boo06 <- bootMer(gm1, mySumm, nsim = 10, use.u = TRUE)

    cbpp$obs <- factor(seq(nrow(cbpp)))
    gm2 <- glmer(cbind(incidence, size - incidence) ~ period +
                 (1 | herd) +  (1|obs),
                 family = binomial, data = cbpp)
    boo03 <- bootMer(gm2, mySumm, nsim = 10)
    boo04 <- bootMer(gm2, mySumm, nsim = 10, use.u = TRUE)
}
load(system.file("testdata","culcita_dat.RData",package="lme4"))
cmod <- glmer(predation~ttt+(1|block),family=binomial,data=culcita_dat)
set.seed(101)
## FIXME: sensitive to step-halving PIRLS tests
## expect_warning(cc <- confint(cmod,method="boot",nsim=10,quiet=TRUE,
##              .progress="txt",PBargs=list(style=3)),"some bootstrap runs failed")

library(parallel)
if (detectCores()>1) {
    ## http://stackoverflow.com/questions/12983137/how-do-detect-if-travis-ci-or-not
    travis <- nchar(Sys.getenv("TRAVIS"))>0
    if(.Platform$OS.type != "windows" && !travis) {
        boo01P <- bootMer(fm1, mySumm, nsim = 10, parallel="multicore", ncpus=2)
    }

    ## works in Solaris from an interactive console but not ???
    ##   via R CMD BATCH
    if (Sys.info()["sysname"] != "SunOS")
        boo01P.snow <- bootMer(fm1, mySumm, nsim = 10, parallel="snow", ncpus=2)
}


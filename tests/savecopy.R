library(lme4)
fn <- "savecopy.rda"
unlink(fn)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
save("fm1",file="savecopy.rda")
rm(fm1)
load(fn)
## should work, but seems to be some problem with
##  comparing RC's after load&restore ... ?
## https://stat.ethz.ch/pipermail/r-devel/2016-August/072944.html ?
## all.equal(fm1@pp,fm1@pp$copy())

fm2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
## should work, but seems to be some problem with
##  comparing RC's after load&restore ... ?
## https://stat.ethz.ch/pipermail/r-devel/2016-August/072944.html ?
## all.equal(fm2@pp,fm2@pp$copy())

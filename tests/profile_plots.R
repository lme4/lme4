library(lme4)
fm1 <- lmer(Reaction~ Days + (Days|Subject), sleepstudy)
## system.time( tpr.fm1 <- profile(fm1, optimizer="Nelder_Mead") )  ## 20 seconds
## save("tpr.fm1",file="../../inst/testdata/tprfm1.RData")
load(system.file("testdata","tprfm1.RData",package="lme4"))
xyplot(tpr.fm1)
splom(tpr.fm1)
densityplot(tpr.fm1, main="densityplot( profile(lmer(..)) )") # does not work


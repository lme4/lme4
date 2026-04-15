load("glmer_basebatch.rda")
require("lme4")
#profile.gm1 <- profile(gm1,signames=FALSE) # does not run
profile.gm2 <- profile(gm2,signames=FALSE)
profile.gm3 <- profile(gm3,signames=FALSE)
#confint0.prof <- confint(profile.gm1)
confint.prof <- confint(profile.gm2)
confint3.prof <- confint(profile.gm3)
confint0.wald <-confint(gm1,method="Wald")
confint.wald <- confint(gm2,method="Wald")
confint3.wald <- confint(gm3,method="Wald")
save(list=c(ls(pattern="confint.*"),
            ls(pattern="profile.*")),
     file="glmer_profbatch.rda")
sessionInfo()

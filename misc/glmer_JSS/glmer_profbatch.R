load("glmer_basebatch.RData")
require("lme4")
profile.gm2 <- profile(gm2,signames=FALSE)
profile.gm3 <- profile(gm3,signames=FALSE)
confint.prof <- confint(profile.gm2)
confint3.prof <- confint(profile.gm3)
confint.wald <- confint(gm2,method="Wald")
confint3.wald <- confint(gm3,method="Wald")
save(list=c(ls(pattern="confint.*"),
            ls(pattern="profile.*")),
     file="glmer_profbatch.RData")
sessionInfo()

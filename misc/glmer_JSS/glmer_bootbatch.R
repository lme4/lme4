load("glmer_basebatch.RData")
require("lme4")
confint.boot <- confint(gm2,method="boot",seed=101,nsim=501)
confint3.boot <- confint(gm3,method="boot",seed=101,nsim=501)
save("confint.boot","confint3.boot",
     file="glmer_bootbatch.RData")
sessionInfo()

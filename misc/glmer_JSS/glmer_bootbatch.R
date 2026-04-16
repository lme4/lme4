load("glmer_basebatch.rda")
require("lme4")
confint0.boot <- confint(gm1,method="boot",seed=101,nsim=501)
confint.boot <- confint(gm2,method="boot",seed=101,nsim=501)
confint3.boot <- confint(gm3,method="boot",seed=101,nsim=501)
save("confint0.boot", "confint.boot","confint3.boot",
     file="glmer_bootbatch.rda")
sessionInfo()

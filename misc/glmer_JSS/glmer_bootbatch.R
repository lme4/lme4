require("lme4")
#load("glmer_basebatch.rda")
# trying to check origins of glmer_basebatch.rda
cbpp2 <- read.csv("cbpp2.csv")
cbpp2 <- transform(cbpp2,period=factor(period),
                   treatment=factor(treatment,
                                    levels=c("Partial/null","Complete","Unknown")))

gm1 <- glmer(incidence/size ~ period + treatment + avg_size + (1 | herd),
             family = binomial,
             data = cbpp2, weights = size)
cbpp2 <- transform(cbpp2,obs=factor(seq(nrow(cbpp2))))    
gm2 <- update(gm1,.~.+(1|obs))  ## herd and observation-level REs
gm3 <- update(gm1,.~.-(1|herd)+(1|obs))  ## observation-level REs only

confint0.boot2 <- confint(gm1,method="boot",seed=101,nsim=501)
confint.boot2 <- confint(gm2,method="boot",seed=101,nsim=501)
confint3.boot2 <- confint(gm3,method="boot",seed=101,nsim=501)

# this is me saving the version we suspect is correct
save("confint0.boot2", "confint.boot2","confint3.boot2",
     file="glmer_bootbatch_test.rda")

load("glmer_bootbatch_test.rda")
load("glmer_bootbatch.rda")

all.equal(confint0.boot, confint0.boot2)

# Commented this out in case I accidentally anything
#save("confint0.boot", "confint.boot","confint3.boot",
#     file="glmer_bootbatch.rda")
sessionInfo()

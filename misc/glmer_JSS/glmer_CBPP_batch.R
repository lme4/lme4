library("lme4")

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

confint0.boot2 <- confint(gm1,method="boot",seed=101,nsim=501,signames=FALSE)
confint.boot2 <- confint(gm2,method="boot",seed=101,nsim=501,signames=FALSE)
confint3.boot2 <- confint(gm3,method="boot",seed=101,nsim=501,signames=FALSE)

# Already saved below.
save("confint0.boot2", "confint.boot2","confint3.boot2",
     file="CBPP_bootbatch.rda")

profile.gm1 <- profile(gm1,signames=FALSE, devtol = 1e-2,
                       verbose = TRUE)
profile.gm2 <- profile(gm2,signames=FALSE, devtol = 1e-2)
profile.gm3 <- profile(gm3,signames=FALSE)
confint1.prof <- confint(profile.gm1)
#confint.prof <- confint(profile.gm2)
confint3.prof <- confint(profile.gm3)

confint0.wald <-confint(gm1,method="Wald",signames=FALSE)
confint.wald <- confint(gm2,method="Wald",signames=FALSE)
confint3.wald <- confint(gm3,method="Wald",signames=FALSE)
save(list=c(ls(pattern="confint.*"),
            ls(pattern="profile.*")),
     file="CBPP_profbatch.rda")
sessionInfo()


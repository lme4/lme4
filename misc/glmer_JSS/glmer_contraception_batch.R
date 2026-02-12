library(mlmRev)
library(lme4)
data(Contraception)
Contraception$ch <- factor(Contraception$livch != 0, labels = c("N","Y"))
# getting the models again...
cm1 <- glmer(use ~ age + I(age^2) + urban + livch + (1|district), 
             Contraception, binomial, nAGQ=0L)
cm2 <- glmer(use ~ age + I(age^2) + urban + ch + (1|district),
             Contraception, binomial, nAGQ=0L)
cm3 <- glmer(use ~ age*ch + I(age^2) + urban + (1|district),
             Contraception, binomial)
cm4 <- glmer(use ~ age*ch + I(age^2) + urban + (urban|district),
             Contraception, binomial)
cm5 <- glmer(use ~ age*ch + I(age^2) + urban + (1|urban:district) + (1|district),
             Contraception, binomial)
cm6 <- glmer(use ~ age*ch + I(age^2) + urban + (1|urban:district),
             Contraception, binomial)
# getting the bootstrap CI
confint.boot.cm1 <- confint(cm1,method="boot",seed=101,nsim=501)
confint.boot.cm2 <- confint(cm2,method="boot",seed=101,nsim=501)
# convergence failed, see if I can increase the tolerance later?
#927 warning(s): Model failed to converge with max|grad| = 0.00200872 (tol = 0.002, component 1)
#See ?lme4::convergence and ?lme4::troubleshooting. (and others)
confint.boot.cm3 <- confint(cm3,method="boot",seed=101,nsim=501)
confint.boot.cm4 <- confint(cm4,method="boot",seed=101,nsim=501)
confint.boot.cm5 <- confint(cm5,method="boot",seed=101,nsim=501)
confint.boot.cm6 <- confint(cm6,method="boot",seed=101,nsim=501)

# saving the info -> commented out for now in case a mistake is made
#save("confint.boot.cm1", "confint.boot.cm2",
#     "confint.boot.cm3", "confint.boot.cm4",
#     "confint.boot.cm5", "confint.boot.cm6",
#     file="glmer_bootbatch2.RData")

#load the file via the command below
#load(file="glmer_bootbatch2.RData")

## another file, likely will save it separately later.

# theta size mismatch??
#profile.cm1 <- profile(cm1,signames=FALSE)
#profile.cm2 <- profile(cm2,signames=FALSE)
profile.cm3 <- profile(cm3,signames=FALSE)
profile.cm4 <- profile(cm4,signames=FALSE)
profile.cm5 <- profile(cm5,signames=FALSE)
profile.cm6 <- profile(cm6,signames=FALSE)

#confint.prof.cm1 <- confint(profile.cm1)
#confint.prof.cm2 <- confint(profile.cm2)
confint.prof.cm3 <- confint(profile.cm3)
confint.prof.cm4 <- confint(profile.cm4)
confint.prof.cm5 <- confint(profile.cm5)
confint.prof.cm6 <- confint(profile.cm6)

# this runs quickly...
confint.wald.cm1 <- confint(cm1,method="Wald",seed=101,nsim=501)
confint.wald.cm2 <- confint(cm2,method="Wald",seed=101,nsim=501)
confint.wald.cm3 <- confint(cm3,method="Wald",seed=101,nsim=501)
confint.wald.cm4 <- confint(cm4,method="Wald",seed=101,nsim=501)
confint.wald.cm5 <- confint(cm5,method="Wald",seed=101,nsim=501)
confint.wald.cm6 <- confint(cm6,method="Wald",seed=101,nsim=501)

# verify the list pattern before running
#save(list=c(ls(pattern="confint.wald.*"),
#            ls(pattern="profile.*")),
#     file="glmer_profbatch2.RData")


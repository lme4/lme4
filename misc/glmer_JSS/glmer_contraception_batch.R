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
#            ls(pattern="profile.*"),
#            ls(pattern="confint.prof.*")),
#     file="glmer_profbatch2.RData")

gm.ar1 <- glmer(incidence/size ~ ar1(1 + herd | period),
                family = binomial,
                data = cbpp, weights = size)
cm.diag <- glmer(use ~ diag(1 + age|district), Contraception, binomial)
cm.hetdiag <- glmer(use ~ diag(1 + age|district, hom = FALSE), 
                    Contraception, binomial)
cm.cs <- glmer(use ~ cs(1 + age|district), Contraception, binomial)
cm.hetcs <- glmer(use ~ cs(1 + age|district, hom = FALSE), Contraception, binomial)

# verify the list pattern before running
#save("gm.ar1", "cm.diag", "cm.hetdiag", "cm.cs", "cm.hetcs",
#     file="glmer_cm_vc.RData")

##### SECOND BATCH
# issue with the plot: need to actually have the real-scaled values for comparisons.

data(Contraception)

# change the data to be 01 instead
Contraception$ch <- ifelse(Contraception$livch != 0, 1, 0)
Contraception$urban <- ifelse(Contraception$urban == "Y", 1, 0)

gelman_scale <- function(dataset, response) {
  new_dataset <- dataset
  for(var in setdiff(names(dataset), response)){
    x <- dataset[[var]]
    if (is.numeric(x) && length(unique(x)) > 2) {
      new_dataset[[var]] <- (x - mean(x, na.rm = TRUE)) / (2 * sd(x, na.rm = TRUE))
    } else if(all(sort(unique(x)) == c(0, 1))){
      new_dataset[[var]] <- x - mean(x, na.rm = TRUE)
    } else {
      new_dataset[[var]] <- x
    }
  }
  return(new_dataset)
}

Contraception2 <- gelman_scale(Contraception, "use")

#######

cm1_2 <- glmer(use ~ age + I(age^2) + urban + livch + (1|district), 
             data = Contraception2, binomial, nAGQ=0L)

cm2_2 <- glmer(use ~ age + I(age^2) + urban + ch + (1|district),
             data = Contraception2, binomial, nAGQ=0L)

cm3_2 <- glmer(use ~ age*ch + I(age^2) + urban + (1|district),
             data = Contraception2, binomial)

cm4_2 <- glmer(use ~ age*ch + I(age^2) + urban + (urban|district),
             data = Contraception2, binomial)

cm5_2 <- glmer(use ~ age*ch + I(age^2) + urban + (1|urban:district) + 
               (1|district), data = Contraception2, binomial)

cm6_2 <- glmer(use ~ age*ch + I(age^2) + urban + (1|urban:district),
             data = Contraception2, binomial)

confint.boot.cm1_2 <- confint(cm1_2,method="boot",seed=101,nsim=501)
confint.boot.cm2_2 <- confint(cm2_2,method="boot",seed=101,nsim=501)
confint.boot.cm3_2 <- confint(cm3_2,method="boot",seed=101,nsim=501)
#5 warning(s): Model failed to converge with max|grad| = 0.00492516 (tol = 0.002, component 1)
confint.boot.cm4_2 <- confint(cm4_2,method="boot",seed=101,nsim=501)
#433 warning(s): Model failed to converge with max|grad| = 0.00204504 (tol = 0.002, component 1)
#See ?lme4::convergence and ?lme4::troubleshooting. (and others)
confint.boot.cm5_2 <- confint(cm5_2,method="boot",seed=101,nsim=501)
#156 message(s): boundary (singular) fit: see help('isSingular')
#107 warning(s): Model failed to converge with max|grad| = 0.00205994 (tol = 0.002, component 1)
confint.boot.cm6_2 <- confint(cm6_2,method="boot",seed=101,nsim=501)
#7 warning(s): Model failed to converge with max|grad| = 0.0272862 (tol = 0.002, component 1)
#See ?lme4::convergence and ?lme4::troubleshooting. (and others)

profile.cm1_2 <- profile(cm1_2,signames=FALSE)
profile.cm2_2 <- profile(cm2_2,signames=FALSE)
profile.cm3_2 <- profile(cm3_2,signames=FALSE)
profile.cm4_2 <- profile(cm4_2,signames=FALSE)
profile.cm5_2 <- profile(cm5_2,signames=FALSE)
profile.cm6_2 <- profile(cm6_2,signames=FALSE)

#confint.prof.cm1_2 <- confint(profile.cm1_2)
#confint.prof.cm2_2 <- confint(profile.cm2_2)
confint.prof.cm3_2 <- confint(profile.cm3_2)
confint.prof.cm4_2 <- confint(profile.cm4_2)
#Warning messages:
#  1: In nextpar(mat, cc, i, delta, lowcut, upcut) :
#  unexpected decrease in profile: using minstep
confint.prof.cm5_2 <- confint(profile.cm5_2)
confint.prof.cm6_2 <- confint(profile.cm6_2)

# haven't run below yet.
confint.wald.cm1_2 <- confint(cm1_2,method="Wald",seed=101,nsim=501)
confint.wald.cm2_2 <- confint(cm2_2,method="Wald",seed=101,nsim=501)
confint.wald.cm3_2 <- confint(cm3_2,method="Wald",seed=101,nsim=501)
confint.wald.cm4_2 <- confint(cm4_2,method="Wald",seed=101,nsim=501)
confint.wald.cm5_2 <- confint(cm5_2,method="Wald",seed=101,nsim=501)
confint.wald.cm6_2 <- confint(cm6_2,method="Wald",seed=101,nsim=501)

save(list=c(ls(pattern="confint.prof.*"), 
            ls(pattern="confint.wald.*")),
     file="glmer_profbatch4.RData")

#save(list=c(ls(pattern="confint.boot.*")),
#     file="glmer_profbatch3.RData")

load("glmer_profbatch3.RData")
load("glmer_profbatch4.RData")

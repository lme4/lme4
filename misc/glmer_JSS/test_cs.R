data("Contraception", package = "mlmRev")
library("lme4")
Contraception$age_s <- drop(scale(Contraception$age))
##  still failing -- opt bad? 
cm.cs <- glmer(use ~ cs(1 + age|district),
               Contraception, binomial,
               start = list(theta=c(0.1, 0, 0.1)))
## problem in mkMerMod/upReCovs/setTheta
cm.hetcs <- glmer(use ~ cs(1 + age|district, hom = FALSE), Contraception, binomial,
                  start = list(theta = c(0.1, 0, 0.1)))
## these seem to be bad because both variances hit zero ... ???????
## (problem still arises with nAGQ=0)

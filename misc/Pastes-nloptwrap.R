library(lme4)
library(testthat)
library(lattice)

options(nwarnings = 5000)# instead of 50, and then use  summary(warnings())

###---- triggered from
"../tests/profile-tst.R"

## 1) "old" default  "bobyqa"
fmoB <- lmer(strength ~ 1 + (cask | batch), data=Pastes,
             control = lmerControl(optimizer = "bobyqa"))
(pfmoB <- profile(fmoB, which = "beta_", alphamax=.001))
xyplot(pfmoB)# nice and easy ..

## 2) new (potential) default "nloptwrap" which is  using
environment(nloptwrap)$defaultControl[["algorithm"]] ## == "NLOPT_LN_BOBYQA"

summary(
    fm <- lmer(strength ~ 1 + (cask | batch), data=Pastes,
               control = lmerControl(optimizer = "nloptwrap",
                                     calc.derivs= FALSE))
)

ls.str(environment(nloptwrap))# showing *its* defaults

pfm <- profile(fm, which = "beta_", alphamax=.001) # 197 warnings for "nloptwrap"
summary(warnings())
str(pfm) # only 3 rows, .zeta = c(0, NaN, Inf) !!!
try( xyplot(pfm) ) ## FIXME or rather the profiling or rather the "wrap on nloptr"

###---- triggered from
"../tests/testthat/test-methods.R"  ## has
## m1 <- lmer(strength ~ 1 + (cask|batch), Pastes)
## ----  but we reuse models fit above

## 1) "bobyqa" : "fine"
set.seed(47); cioB <- CI.boot(fmoB)
cioB # fine (all 6+1 = 7 "sigma"s are significant on 5% level)

## but
## 2) "nloptwrap" ("NLOPT_LN_BOBYQA")  fails badly:
set.seed(47)
try(ci <- CI.boot(fm)) # prints, then Error :
## [1] "All values of t are equal to  0.675470129506501 \n Cannot calculate confidence intervals"
## Error in dimnames(citab) <- list(names(bb[["t0"]]), format.perc(a, 3)) :
##   length of 'dimnames' [1] not equal to array extent

try(ci <- confint(fm,  method="boot", nsim=10, seed=101))

(lme4:::confint.merMod) # shows that bootMer() is used,
## then boot::boot.ci() [which gives the message "All ... t are equal .." + the error]





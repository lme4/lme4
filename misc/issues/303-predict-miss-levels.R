## https://github.com/lme4/lme4/issues/303#issuecomment-104406146

## Related comments from Russ Lenth:

## > Looking at the code in lme4 related to merMod objects, I noticed that
## > there did not appear to be any directly recoverable information on the
## > nature of rank deficiency when it exists. So I did a little test, and
## > sure enough, I verified that there's a problem...
## >
## > Load a test dataset and identify cases for one factor combination:


data(Oats, package = "nlme")
excl <-  with(Oats, which(nitro == .2 &
                          Variety == "Marvellous"))
excl
## [1] 10 22 34 46 58 70

## Fit a model, excluding the identified subset

library("lme4")
Oats.lmer <-  lmer(yield ~ Variety * factor(nitro) + (1|Block/Variety),
                   data = Oats, subset = -excl)
## "message":
## fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

## Get a summary showing estimability issue via `lsmeans`
library("lsmeans")
rg <- ref.grid(Oats.lmer)
(nd <- summary(rg))
## Variety     nitro prediction  SE       df
## Golden Rain   0.0  80.00000    8.839195 16.24
## Marvellous    0.0  86.66667    8.839195 16.24
## Victory       0.0  71.50000    8.839195 16.24
## Golden Rain   0.2  98.50000    8.839195 16.24
## Marvellous    0.2  NA          NA       NA
## Victory       0.2  89.66667    8.839195 16.24
## Golden Rain   0.4 114.66667    8.839195 16.24
## Marvellous    0.4 117.16667    8.839195 16.24
## Victory       0.4 110.83333    8.839195 16.24
## Golden Rain   0.6 124.83333    8.839195 16.24
## Marvellous    0.6 126.83333    8.839195 16.24
## Victory       0.6 118.50000    8.839195 16.24


## The `nd` result is a data frame with values of the fixed-effects
try(
predict(Oats.lmer, newdata = nd)
)## Error in X %*% fixef(object) : non-conformable arguments

## MM: error above no longer happens, but rather -- a NEW ERROR:
##----  Error in eval(expr, envir, enclos) : object 'Block' not found
## and
## and this error msg is correct :

## ===> Fix our "newdata":
nd2 <- rbind(cbind(nd, Block = as.factor("III")),
             cbind(nd, Block = as.factor("IV" )))

p2 <- predict(Oats.lmer, newdata = nd2)
## no problem here anymore

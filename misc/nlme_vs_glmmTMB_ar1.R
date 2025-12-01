library(glmmTMB)
library(nlme)

### INSTRUCTIONS:
# 1. glmmTMB, AR1, dispersion 0
# 2. glmmTMB, OU, dispersion 0
# 3. lme, corAR1
# 4. lme, corCAR1
# 5. lme, corEXP no nugget

### EXPECTATIONS:
# (i)   expect 1 and 3 to be the same
# (ii)  expect 2 and 5 to be the same
# (iii) expect 3 and 4 to be the same

# could be bad dataset; return to soy example
data("Soybean", package = "MEMSS")

# 1. glmmTMB, AR1, dispersion 0
Soybean$Timef <- factor(Soybean$Time, ordered = TRUE)
glmmTMB.ar1.disp0 <- glmmTMB(weight ~ Time + ar1(0 + Timef | Plot), 
                       Soybean, REML = FALSE, dispformula = ~0)

# 2. glmmTMB, ou, dispersion 0
Soybean$Timeff <- numFactor(Soybean$Time)
glmmTMB.ou.disp0 <- glmmTMB(weight ~ Time + ou(0 + Timeff | Plot), 
                          Soybean, REML = FALSE, dispformula = ~0)

# 3. nlme, corAR1
nlme.corAR1 <- gls(weight ~ Time,
                correlation = corAR1(form = ~ Time | Plot),
                data = Soybean, method = "ML")

# 4. nlme, corCAR1
nlme.corCAR1 <- gls(weight ~ Time,
                correlation = corCAR1(form = ~ Time | Plot),
                data = Soybean, method = "ML")

# 5. nlme, corEXP
nlme.corEXP <- gls(weight ~ Time,
                correlation = corExp(form = ~ Time | Plot, nugget = TRUE),
                data = Soybean, method = "ML")

# Do comparisons
mod_list <- list(glmmTMB.ar1.disp0 = glmmTMB.ar1.disp0,
                 glmmTMB.ou.disp0 = glmmTMB.ou.disp0,
                 nlme.corAR1 = nlme.corAR1,
                 nlme.corCAR1 = nlme.corCAR1,
                 nlme.corEXP = nlme.corEXP)

## Comparing certain values between them...
lliks <- sapply(mod_list, function(x) -1*c(logLik(x)))
lliks

sigma <- sapply(mod_list, function(x) sigma(x))
sigma

vcovs <- sapply(mod_list, function(x) vcov(x))
vcovs

fixefs <- sapply(mod_list, function(x){
  if(class(x) == "glmmTMB") fixef(x)$cond else coef(x)
  })
fixefs

## note: varcorrs is going to look fundamentally different; can ignore below..?
varcorrs <- sapply(mod_list, 
                function(x){
                  if(class(x) == "glmmTMB"){
                    VarCorr(x)$cond[[1]]
                  } else vcov(x)
                })
varcorrs

### Exploration Conclusion
# (i) 1 and 3 are NOT the same
# (ii) 2 and 5 are similar except for sigma??
# (ii) 3 and 4 are NOT the same


# Tested with simulated data below; doesn't work for nlme corAR1 for some reason,
# so ignore that.
simGroup <- function(g, n=6, phi=0.6) {
  x <- MASS::mvrnorm(mu = rep(0,n),
                     Sigma = phi^as.matrix(dist(1:n)) )  
  y <- x + rnorm(n)                              
  times <- factor(1:n)
  group <- factor(rep(g,n))
  data.frame(y, times, group)
}

set.seed(1)
dat <- do.call("rbind", lapply(1:20, simGroup))

# 1. glmmTMB, AR1, dispersion 0
glmmTMB.ar1.disp0 <- glmmTMB(y ~ times + ar1(0 + times | group), 
                             data = dat, dispformula=~0)

# 2. glmmTMB, OU, dispersion 0
dat$timesff <- numFactor(dat$times)
glmmTMB.ou.disp0 <- glmmTMB(y ~ times + ou(0 + timesff | group), 
                            data = dat, REML = FALSE, dispformula=~0)

# lme, corAR1
# DOES NOT WORK BELOW...
nlme.corAR1 <- gls(y ~ times,
                   correlation = corAR1(form = ~ times | group),
                   data = dat, method = "ML")

# lme, corCAR1
nlme.corCAR1 <- gls(y ~ times,
                    correlation = corCAR1(form = ~ times | group),
                    data = dat, method = "ML")

# lme, corEXP no nugget
nlme.corEXP <- gls(y ~ times,
                   correlation = corExp(form = ~ times | group, nugget = FALSE),
                   data = dat, method = "ML")


# Do comparisons
mod_list <- list(glmmTMB.ar1.disp0 = glmmTMB.ar1.disp0,
                 glmmTMB.ou.disp0 = glmmTMB.ou.disp0,
                 nlme.corCAR1 = nlme.corCAR1,
                 nlme.corEXP = nlme.corEXP)

## Comparing certain values between them...
lliks <- sapply(mod_list, function(x) -1*c(logLik(x)))
lliks
# they're all similar except nlme.corEXP

sigma <- sapply(mod_list, function(x) sigma(x))
sigma
# glmmTMB's are similar with each other, and same with nlme's

fixefs <- sapply(mod_list, function(x){
  if(class(x) == "glmmTMB") fixef(x)$cond else coef(x)
})
fixefs
# all relatively similar...

vcovs <- sapply(mod_list, function(x) vcov(x))
vcovs

### Exploration Conclusion
# (ii) 2 and 5 are similar except for sigma??

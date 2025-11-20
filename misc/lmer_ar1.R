library(lme4)
library(glmmTMB)
library(nlme)
library(broom.mixed)
library(tidyverse)

data("Soybean", package = "MEMSS")

library(ggplot2)
ggplot(Soybean, aes(Time, weight, group = Plot)) +
  geom_point() + geom_line() +
  ## scale_y_log10()
  scale_y_continuous(trans = "sqrt")

Soybean$Timef <- factor(Soybean$Time, ordered = TRUE)

soy.lme4 <- lmer(weight ~ Time + ar1(0 + Timef | Plot, hom = TRUE), 
            Soybean, REML = FALSE, 
            control = lmerControl(check.nobs.vs.nRE = "ignore"))

soy.glmmTMB <- glmmTMB(weight ~ Time + ar1(0 + Timef | Plot), 
                       Soybean, REML = FALSE)

## ornstein-uhlenbeck
Soybean$Timeff <- numFactor(Soybean$Time)
soy.glmmTMB.OU <- glmmTMB(weight ~ Time + ou(0 + Timeff | Plot), 
                       Soybean, REML = FALSE)

soy.nlme <- gls(weight ~ Time,
                correlation = corExp(form = ~ Time | Plot, nugget = TRUE),
                data = Soybean, method = "ML")

## Comparing lme4 to glmmTMB ########################
# likelihoods are similar; below returns TRUE
all.equal(-soy.glmmTMB$fit$objective, logLik(soy.lme4), 
          check.attributes = FALSE,
          check.class = FALSE)

mod_list <- list(glmmTMB_ar1 = soy.glmmTMB,
                 lme4 = soy.lme4,
                 nlme = soy.nlme,
                 glmmTMB_ou = soy.glmmTMB.OU)

lliks <- sapply(mod_list, function(x) -1*c(logLik(x)))
print(lliks)
lliks - lliks[1]

fixef_sum <- purrr::map_dfr(mod_list, tidy, .id = "model", effects = "fixed") |>
  mutate(across(model, ~ forcats::fct_inorder(factor(.)))) |>
  select(model, term, estimate, std.error) |>
  pivot_longer(c(estimate, std.error)) |>
  arrange(term, name, model) |>
  group_by(term, name)

ggplot(fixef_sum, aes(value, model)) +
  facet_wrap(~term + name, scales = "free") +
  geom_point()

all.equal(fixef(soy.glmmTMB)$cond, fixef(soy.lme4),
          check.attributes = FALSE,
          check.class = FALSE)

all.equal(vcov(soy.lme4), vcov(soy.glmmTMB)$cond, check.attributes = FALSE,
          check.class = FALSE)
# "Mean relative difference: 0.0001230642"
# not 5e-5, but should be ok?

## Comparing lme4 to nlme ###########################
# differences are relatively large...
all.equal(logLik(soy.lme4), logLik(soy.nlme), check.attributes = FALSE,
          check.class = FALSE)
# [1] "Mean relative difference: 0.012905"
all.equal(vcov(soy.lme4), vcov(soy.nlme), check.attributes = FALSE,
          check.class = FALSE)
# "Mean relative difference: 0.3676962"

#################################
## second version (try using REML!!!)

soy.lme4.REML <- lmer(weight ~ Timef + ar1(0 + Timef | Plot, hom = TRUE), 
            Soybean, control = lmerControl(check.nobs.vs.nRE = "ignore"))

soy.nlme.REML <- gls(weight ~ Timef,
                correlation = corExp(form = ~ Time | Plot, nugget = TRUE),
                data = Soybean)
# doesn't improve...
all.equal(logLik(soy.lme4.REML), logLik(soy.nlme.REML), 
          check.attributes = FALSE, check.class = FALSE)
# [1] "Mean relative difference: 0.01200982"
all.equal(vcov(soy.lme4.REML), vcov(soy.nlme.REML), 
          check.attributes = FALSE, check.class = FALSE)
# [1] "Mean relative difference: 0.2913537"

## TO DO
##  (1) why are OU and ar1 so different?
##  (2) comparing glmmTMB with nlme, I can turn off the nugget in glmmTMB by
##    using dispformula = ~0, then I can compare to (1) straight corAR1 (2)
##    corCAR1 without nugget
##  (3) simulate data like this with known true values, fit many times, and see
##   what the actual coverages are for the different methods ...

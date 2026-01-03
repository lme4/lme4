
## exploring https://stackoverflow.com/questions/37466771/using-profile-and-boot-method-within-confint-option-with-glmer-model
library(lme4)
library(tinyplot)
library(ggplot2); theme_set(theme_bw())
library(ggalt)
library(broom.mixed)

df2 <- data.frame(
  prop1 = c(0.46, 0.471, 0.458, 0.764, 0.742, 0.746, 0.569, 0.45,    0.491,    0.467, 0.464, 
            0.556, 0.584, 0.478, 0.456, 0.46, 0.493, 0.704, 0.693, 0.651), 
  prop2 = c(0.458, 0.438, 0.453, 0.731, 0.738, 0.722, 0.613, 0.498, 0.452, 0.451, 0.447,
            0.48, 0.576, 0.484, 0.473, 0.467, 0.467, 0.722, 0.707, 0.709),
  site = c(1,1,2,3,3,3,4,4,4,4,4,5,5,5,6,6,7,8,8,8)
)
df2$site <- factor(df2$site)

plt(prop2~prop1|site, data = df2, pch = 16)

ggplot(df2, aes(prop1, prop2)) +
  geom_point(aes(colour = site)) +
  geom_encircle(aes(fill = site), alpha = 0.3, colour = NA)

## there's not much variability to work with
## can't currently fit profiles for GLMMs with estimated scale parameters ...
## maybe finding local minimum anyway?

form <- prop2 ~ prop1 + (1|site)
fit_glmer <- glmer(form, data=df2, family=gaussian(link="logit"),
                   control = glmerControl(optimizer = "nloptwrap",
                                          optCtrl = list(itmax = 1000)))

bootFun <- function(x, signames = FALSE) {
  ## repeat useSc so this can be used upstream/outside if necessary
  useSc <- as.logical(x@devcomp$dims[["useSc"]])
  pp <- getME(x, "par")
  sig <- if (!useSc) NULL else sigma(x)
  ## hard-code profscale -- users can specify FUN if they want ...
  repars <- lme4:::convParToProfPar(pp, x, profscale = "sdcor" , sc = sig)
  ## FIXME: names are (may be?) still in the wrong order ...
  names(repars) <- lme4:::profnames(x, signames, useSc=useSc)
  c(repars, fixef(x))
}

## can we handle out-of-range stuff better?
bb <- bootMer(fit_glmer, FUN = bootFun, nsim = 2500, parallel = "multicore", ncpus = 8,
              seed = 101)
bci <- confint(bb)

bci_long <- as.data.frame.table(bb$t) |> setNames(c("junk", "par", "val"))
plt(~val, type = "histogram", data = bci_long, facet = ~par,
    facet.args = list(free = TRUE))

ss <- bci_long |> subset(par == "sd_(Intercept)|site")
hist(ss$val)

library(glmmTMB)
fit_glmmTMB <- glmmTMB(form,
                       data=df2,
                       family=gaussian(link="logit"))
bootFun2 <- function(x) {
  tt <- tidy(x)
  nm <- with(tt, paste(ifelse(is.na(group), "", group), term, sep = "."))
  setNames(tt$estimate, nm)
}

bb2 <- bootMer(fit_glmmTMB, FUN = bootFun2, nsim = 2500, parallel = "multicore", ncpus = 8,
               seed = 101)
confint(bb2)  ## site SD: approx ~ 0-0.085)

pp <- TMB:::tmbprofile(fit_glmmTMB$obj, name = 4, adaptive = FALSE, h = 1e-1)
pp$relval <- pp$value - min(pp$value)
plot(pp[,1], pp$relval, type = "b")
cutoff <- qchisq(0.95, df = 1)/2
abline(h = cutoff, col = 2, lty = 2)
exp(-2.09)

library(brms)
set.seed(101)
fit_brms <- brm(form, data=df2, family=gaussian(link="logit"), control = list(adapt_delta = 0.99))
tidy(fit_brms, conf.int = TRUE) ## 0.004 - 0.219

## conclusion: glmer is probably finding a local max (the best fit is probably
##  actually singular), which in turn throws off the PB
## glmmTMB gets a slightly higher log-likelihood (0.3 units)
## glmmTMB and brms get sufficiently similar CIs that I think I believe them
## brms still gets a couple of divergent transitions with adapt_delta = 0.9

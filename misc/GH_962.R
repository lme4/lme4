form <- out_bin_1 ~ ftime + trt + ar1(ftime + 0 | id_cluster) +
  (1 | id_individual)

form2 <- update(form, . ~ . - ar1(ftime + 0 | id_cluster) +
                        cs(ftime + 0 | id_cluster))

form3 <- update(form, . ~ . - ar1(ftime + 0 | id_cluster) +
                        (1|id_cluster/ftime))

library(glmmTMB)
library(lme4)
stopifnot(packageVersion("lme4") >= "2.0.2")
set.seed(101)
sfun <- function(n) factor(sample(1:n, replace=TRUE, size = 1465))
dd <- data.frame(id_cluster = sfun(7), id_individual = sfun(486),
                 ftime = sfun(4), trt = sfun(2))
np <- list(beta = c(0, rep(0.1, 3), 1), theta = rep(0, 3))
dd$out_bin_1 <- simulate_new(form[-2], newdata = dd,
             newparams = np,
             family = binomial)[[1]]

## devtools::load_all()
## debug(optimizeGlmer)
## looks OK after nAGQ0: environment(devfun)$pp$beta(1)
## but not after nAGQ1

## still not sure what's going on ...
## cross-check with other forms?
lme4_1 <- glmer(form, data = dd, family = binomial)
fixef(lme4_1)
lme4_1@optinfo$val

lme4_1B <- update(lme4_1, nAGQ=0)
fixef(lme4_1B)

lme4_2 <- update(lme4_1, form = form2)

lme4_2@optinfo$val

lme4_3 <- update(lme4_1, form = form3)
fixef(lme4_3)
lme4_3@optinfo$val

## only happens for 'new-style' (structured); something happens
## (or fails to happen) in the parameter mapping step?

glmmTMB_1 <- glmmTMB(form, data = dd, family = binomial)
c(logLik(lme4_1) - logLik(glmmTMB_1))
logLik(lme4_1)
logLik(glmmTMB_1)


dd$out_gauss <- simulate_new(form[-2], newdata = dd,
             newparams = c(np, list(betadisp=0)),
             family = gaussian)[[1]]

form2 <- update(form, out_gauss ~ .)
lme4_2 <- lmer(form2, data = dd, REML = FALSE)
glmmTMB_2 <- glmmTMB(form2, data = dd, family = gaussian)
fixef(lme4_2)

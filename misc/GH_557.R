## https://stats.stackexchange.com/questions/447029/parameter-recovery-gamma-distribution-and-model-with-a-random-intercept#447029
library(lme4)
library(GLMMadaptive)
library(glmmTMB)
library(broom.mixed)
library(dplyr)
library(ggplot2); theme_set(theme_bw())

simulate_gamma_mixed <- function (n,      # number of subjects
                                  seed = 1,
                                  betas = c(0.2, -0.1, 0.15),
                                  sigma_b = 0.2,
                                  scale = 0.15,
                                  K = 30) # number of measurements per subject
{
    if (!exists(".Random.seed", envir = .GlobalEnv)) 
        runif(1)
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
    # we construct a data frame with the design: 
    DF <- data.frame(id = rep(seq_len(n), each = K),
                     time = gl(3, 10, length = n * K,
                               labels = paste0("T", 0:2)))
    # design matrix for the fixed effects
    X <- model.matrix(~ time, data = DF)
    # we simulate random effects
    b <- rnorm(n, sd = sigma_b)
    # linear predictor
    eta <- c(X %*% betas + b[DF$id])
    # we simulate Gamma longitudinal data with log link
    DF$y <- rgamma(n * K, shape = exp(eta) / scale, scale = scale)
    DF
}

##########################################################################################
##########################################################################################

library("GLMMadaptive") # you need version >= 0.6-9
library("lme4")
library("glmmTMB")
library("broom.mixed")

M <- 100 # number of simulations
betas_glmer <- betas_mixmod <- betas_glmmTMB <- matrix(as.numeric(NA), M, 3)
sigma_b_glmer <- sigma_b_mixmod <- sigma_b_glmmTMB <- numeric(M)
pb <- txtProgressBar(max = M, style = 3)

DF_m <- simulate_gamma_mixed(n = 50, seed = 2020)
fm_0 <- glmer(y ~ time + (1 | id), data = DF_m, 
              family = Gamma(link = "log"), nAGQ = 15)

fm_1 <- glmer(y ~ time + (1 | id), data = DF_m, 
              family = Gamma(link = "log"), nAGQ = 15,
              start = list(fixef = c(0.2, -0.1, 0.15),
                           theta = 0.5))
## still converges badly
VarCorr(fm_1)
## "can't (yet) profile GLMMs with non-fixed scale parameters", sigh ...
dfun <- getME(fm_1, "devfun")
## profile by hand
p0 <- unlist(getME(fm_1, c("theta", "beta")))
stopifnot(all.equal(dfun(p0), c(-2*logLik(fm_1))))

pfun <- function(beta, theta) {
    dfun(c(theta, beta))
}
pfun(p0[-1], p0[1])
tvec_upr <- seq(p0[1], p0[1]+0.2, by = 0.01)
sigma(fm_1)*getME(fm_1, "theta")

## ugh, have to recapitulate sigma() calculation, from internals
get_sigma <- function(obj) {
    pp <- obj@pp$copy()
    resp <- obj@resp
    sqrLenU <- pp$sqrL(1)
    wrss    <- resp$wrss()
    n <- length(resp$y)
    pwrss   <- wrss + sqrLenU
    sqrt(pwrss/n)
}
stopifnot(identical(get_sigma(fm_1), sigma(fm_1)))

## need to multiply sigma by theta for each optimization to get
##  corresponding
prof_upr <- vapply(tvec_upr, \(t) nloptwrap(par = p0[-1],
                                fn = pfun,
                                theta = t,
                                lower = rep(-Inf, 3),
                                upper = rep(Inf, 3))$fval,
       FUN.VAL = numeric(1))
tvec_lwr <- seq(p0[1], 0.001, by = -0.01)
prof_lwr <- vapply(tvec_lwr, \(t) nloptwrap(par = p0[-1],
                                            fn = pfun,
                                            theta = t,
                                            lower = rep(-Inf, 3),
                                            upper = rep(Inf, 3))$fval,
       FUN.VAL = numeric(1))
prof <- c(rev(prof_lwr), prof_upr)
tvec <- c(rev(tvec_lwr), tvec_upr)
plot(tvec*sigma(fm_1), prof)  ## approximate scaling

tidy(fm_0, effects = "ran_pars") |>
    dplyr::filter(group=="id") |>
    dplyr::pull(estimate)

## GLMMadaptive tidy method needs to know about 'effects' ...
## tidy(gm_m, effects = "ran_pars")

for (m in seq_len(M)) {
    setTxtProgressBar(pb, m)
    DF_m <- simulate_gamma_mixed(n = 50, seed = 2020 + m)
    fm_m <- refit(fm_0, DF_m$y)
    betas_glmer[m, ] <- fixef(fm_m)
    sigma_b_glmer[m] <- sqrt(VarCorr(fm_m)$id[1, 1])
    gm_m <- mixed_model(y ~ time, random = ~ 1 | id, data = DF_m, 
                        family = Gamma(link = "log"), nAGQ = 15)
    betas_mixmod[m, ] <- fixef(gm_m)
    sigma_b_mixmod[m] <- sqrt(gm_m$D)
    ## Laplace approx only
    tm_m <- glmmTMB(y ~ time + (1 | id), data = DF_m, 
                    family = Gamma(link = "log"))
    betas_glmmTMB[m, ] <- fixef(tm_m)$cond
    sigma_b_glmmTMB[m] <- sqrt(VarCorr(tm_m)$cond$id[1, 1])
}
close(pb)

cbind(glmer = colMeans(betas_glmer),
      mixmod = colMeans(betas_mixmod),
      glmmTMB = colMeans(betas_glmmTMB))

c(glmer = mean(sigma_b_glmer),
  mixmod = mean(sigma_b_mixmod),
  glmmTMB = mean(sigma_b_glmmTMB))



## TO DO:
##  COMPARE PARAMETERIZATIONS: check inst/testdata/lme-tst-funs.R (for Gamma sim),
##      tests/glmmExt.R
##      gSim uses a default sigma_B of 1
##     note that sigma is the "the square root of the residual deviance per degree of freedom in more general models"
##     but this shouldn't affect theta results?
##   simplify simulation code (use simulate()) ?
##   think about scale specification -- seems odd but it's working ...
##   extend broom.mixed to handle everything?
##   consider alternative starting values for lme4?
##   check sigma to see if it is similarly/counterbalancingly off?
    
## simplified example
dd <- data.frame(f = factor(rep(1:30, each = 20)))
set.seed(101)
dd$y <- simulate(~ 1 + (1|f), family = Gamma(link = "log"),
                 newparams = list(beta = 0, theta = 1, sigma = 1),
                 seed = 101,
                 newdata = dd)[[1]]
g0 <- glmer(y ~ 1 + (1|f), family = Gamma(link = "log"), data = dd)
getME(g0, "theta")

set.seed(101)
tvals <- replicate(500,
{
    cat(".")
    suppressMessages(dd$y <- simulate(~ 1 + (1|f), family = Gamma(link = "log"),
                                      newparams = list(beta = 0, theta = 1, sigma = 1),
                                      newdata = dd)[[1]]
                     )
    g0 <- glmer(y ~ 1 + (1|f), family = Gamma(link = "log"), data = dd)
    getME(g0, "theta")
})

hist(tvals, breaks = 100)
## looks fine (slightly downward biased)


library(lme4)
set.seed(101)
d <- expand.grid(block=LETTERS[1:26],rep=1:100)
d$x <- runif(nrow(d))
reff_f <- rnorm(length(levels(d$block)),sd=1)
## need intercept large enough to avoid negative values
d$eta0 <- 4+3*d$x  ## version without random effects
d$eta <- d$eta0+reff_f[d$block]
## inverse link
d$mu <- 1/d$eta
d$y <- rgamma(nrow(d),scale=d$mu/2,shape=2)

(gm1 <- glmer(y ~ 1|block, d, Gamma, nAGQ=25L, compDev=FALSE))
## debug(lme4:::glmerPwrssUpdate)
## FIXME: hangs within optwrap, on second round of optimization:
if (FALSE) {
    gm2 <- glmer(y ~ 1|block, d, Gamma, nAGQ=25L, verbose=10)
    try({gm3 <- glmer(y ~ 1|block, d, Gamma, nAGQ=25L); print(gm1)})
}
## (NM) 1: f = inf at 1.05128 5.46369


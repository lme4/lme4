## original code was designed to detect segfaults/hangs from error handling

library(lme4)
set.seed(101)
d <- expand.grid(block = LETTERS[1:26],
                 rep = 1:100)
d$x <- runif(nrow(d))
reff_f <- rnorm(length(levels(d$block)),sd=1)
## need intercept large enough to avoid negative values
d$eta0 <- 4+3*d$x  ## version without random effects
d$eta <- d$eta0+reff_f[d$block]
## inverse link
d$mu <- 1/d$eta
d$y <- rgamma(nrow(d), scale=d$mu/2, shape=2)

if (.Platform$OS.type != "windows") {
gm0     <- glmer(y ~      1|block,  d, Gamma)
gm0.A25 <- glmer(y ~      1|block,  d, Gamma, nAGQ=25L)
gm1     <- glmer(y ~ x + (1|block), d, Gamma)
gm1.A25 <- glmer(y ~ x + (1|block), d, Gamma, nAGQ=25L)

## strange things happening for logLik  ==> AIC, etc for nAGQ ???
anova(gm0, gm1)
anova(gm0, gm0.A25)
anova(gm1, gm1.A25)

summary(gm1) # "fine"
summary(gm1.A25) # Inf logLik etc ?

}

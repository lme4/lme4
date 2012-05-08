## subsetted from glmmExt.R, for convenience
library(lme4)

set.seed(101)
d <- expand.grid(block=LETTERS[1:26], rep=1:100, KEEP.OUT.ATTRS = FALSE)
d$x <- runif(nrow(d))  ## sd=1
reff_f <- rnorm(length(levels(d$block)),sd=1)
## need intercept large enough to avoid negative values
d$eta0 <- 4+3*d$x  ## fixed effects only
d$eta <- d$eta0+reff_f[d$block]
dBc <- d
cc <- binomial(link="cloglog")
dBc$mu <- cc$linkinv(d$eta)
dBc$y <- rbinom(nrow(d),dBc$mu,size=1)

if (FALSE) {
    ## debug(lme4:::glmerPwrssUpdate)
    ## if we set compDev=FALSE we get
    ##   pdev eventually going to NaN in RglmerWrkIter/glmerPwrssUpdate
    gBc1 <- glmer(y ~ 1 + (1|block), data=dBc,
                  family=binomial(link="cloglog"), verbose= 3,
                  compDev=FALSE)
    gBc2 <- glmer(y ~ x + (1|block), data=dBc,
                  family=binomial(link="cloglog"), verbose= 3)
}





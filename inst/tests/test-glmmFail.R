library("testthat")
library("lme4")

set.seed(101)
d <- expand.grid(block=LETTERS[1:26], rep=1:100, KEEP.OUT.ATTRS = FALSE)
d$x <- runif(nrow(d))  ## sd=1
reff_f <- rnorm(length(levels(d$block)),sd=1)
## set intercept large enough to create a constant response
d$eta0 <- 4+3*d$x  ## fixed effects only
d$eta <- d$eta0+reff_f[d$block]
dBc <- d
cc <- binomial(link="cloglog")
dBc$mu <- cc$linkinv(d$eta)
dBc$y <- rbinom(nrow(d),dBc$mu,size=1)

context("glmer failure because of a constant response")
test_that("constant response", {
    expect_that(glmer(y ~ 1 + (1|block), data=dBc,
                      family=binomial(link="cloglog")),
                throws_error())
})


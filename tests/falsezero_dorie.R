## test of false zero problem reported by Vince Dorie
## (no longer occurs with current development lme4)
## https://github.com/lme4/lme4/issues/17
library(lme4)

sigma.eps <- 2
sigma.the <- 0.75
mu <- 2

n <- 5
J <- 10
g <- gl(J, n)

set.seed(1)

theta <- rnorm(J, 0, sigma.eps * sigma.the)
y <- rnorm(n * J, mu + theta[g], sigma.eps)
lmerFit <- lmer(y ~ 1 + (1 | g), REML = FALSE, verbose=TRUE)

y.bar <- mean(y)
y.bar.j <- sapply(1:J, function(j) mean(y[g == j]))
S.w <- sum((y - y.bar.j[g])^2)
S.b <- n * sum((y.bar.j - y.bar)^2)
R <- S.b / S.w

sigma.the.hat <- sqrt(max((n - 1) * R / n - 1 / n, 0))
stopifnot(all.equal(sigma.the.hat,lme4Sigma <- unname(getME(lmerFit,"theta")),
                    tolerance=2e-5))

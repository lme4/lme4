library("testthat")
library("lme4")

context("summarizing/printing models")
test_that("lmer", {
  set.seed(0)
  J <- 8
  n <- 10
  N <- J * n
  beta <- c(5, 2, 4)
  u <- matrix(rnorm(J * 3), J, 3)
  
  x.1 <- rnorm(N)
  x.2 <- rnorm(N)
  g <- rep(1:J, rep(n, J))
  
  y <- 1   * (beta[1] + u[g,1]) +
       x.1 * (beta[2] + u[g,2]) +
       x.2 * (beta[3] + u[g,3]) +
       rnorm(N)

  summary(lmer(y ~ x.1 + x.2 + (1 + x.1 | g)))
  summary(lmer(y ~ x.1 + x.2 + (1 + x.1 + x.2 | g)))
})

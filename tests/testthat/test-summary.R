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

  tmpf <- function(x) capture.output(print(summary(x),digits=1))
  tfun <- function(cc) {
      w <- grep("Fixed effects:",cc)
      cc[w:length(cc)]
  }
  cc1 <- tmpf(lmer(y ~ x.1 + x.2 + (1 + x.1 | g)))
  cc2 <- tmpf(lmer(y ~ x.1 + x.2 + (1 + x.1 + x.2 | g)))
  expect_equal(tfun(cc1),
               c("Fixed effects:",
                 "            Estimate Std. Error t value", 
                 "(Intercept)      5.4        0.5      12",
                 "x.1              1.9        0.4       5", 
                 "x.2              4.0        0.1      28",
                 "",
                 "Correlation of Fixed Effects:", 
                 "    (Intr) x.1   ",
                 "x.1 -0.019       ",
                 "x.2  0.029 -0.043"
                 ))

  expect_equal(tfun(cc2),
               c("Fixed effects:",
                 "            Estimate Std. Error t value", 
                 "(Intercept)      5.4        0.4      12",
                 "x.1              2.0        0.4       5", 
                 "x.2              4.0        0.3      15",
                 "",
                 "Correlation of Fixed Effects:",
                 "    (Intr) x.1   ",
                 "x.1 -0.069       ",
                 "x.2  0.136 -0.103"
                 ))
})

library("testthat")
try(detach("package:lmerTest"), silent = TRUE)
library("lme4")

testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1

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
      w <- grep("Fixed effects:", cc)
      cc[w:length(cc)]
  }
  C1 <- lmerControl(optimizer="nloptwrap",
                    optCtrl=list(xtol_abs=1e-6, ftol_abs=1e-6))
  m1 <- lmer(y ~ x.1 + x.2 + (1 + x.1 | g), control=C1)
  m2 <- lmer(y ~ x.1 + x.2 + (1 + x.1 + x.2 | g), control=C1)
  cc1 <- tmpf(m1)
  cc2 <- tmpf(m2)
  ## FIXME: correlation of fixed effects printed inconsistently.
  ## If (1) LME4_TEST_LEVEL == 100 *and* after running all of prior
  ##   tests, something (package load? options setting?) changes
  ##   so that the fixed-effect correlations are no longer printed
  ##   out, and this test fails
  ## would like to sort this out but realistically not sure it's worth it?
  t1 <- tfun(cc1)
  vv <- vcov(m1)
  ss <- sessionInfo()
  save(m1, vv, ss, file = sprintf("test-summary_testlevel_%d.rda", testLevel))
  expect_equal(t1,
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

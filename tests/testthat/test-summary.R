try(detach("package:lmerTest"), silent = TRUE)

testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1

#context("summarizing/printing models")
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

## Tests with regards to auto-scaling; making sure we get expected behaviours.

set.seed(1)
sleepstudy$var1 = runif(nrow(sleepstudy), 1e6, 1e7)

scf1 <- lmer(Reaction ~ var1 + Days + (Days | Subject), 
              control = lmerControl(autoscale = TRUE), sleepstudy)
scf2 <- suppressWarnings(lmer(Reaction ~ var1 + Days + (Days | Subject), 
             control = lmerControl(autoscale = FALSE), sleepstudy))

test_that("Ensuring we get the internal scale for X with getME", {
  res1 <- getME(scf1, "X")[180, ]
  valtest1 <- c(1, 1.456867279040504, 1.562340900990911)
  expect_equal(unname(res1), valtest1, tolerance =1e-10)
})

test_that("Ensuring we get the internal scale for beta with getME", {
  res2 <- getME(scf1, "beta")
  valtest2 <- c(298.507891666666751, 1.621716753196637, 30.161580776828778)
  expect_equal(res2, valtest2, tolerance =1e-10)
})

test_that("model.matrix should provide unscaled version at default", {
  res3 <- model.matrix.merMod(scf1)[180, ]
  valtest3 <- c(1, 9127734.503475949, 9)
  expect_equal(unname(res3), valtest3, tolerance =1e-6)
  
  res4 <- model.matrix.merMod(scf1, noScale = TRUE)[180, ]
  valtest4 <- c(1, 1.45686727904050395, 1.5623409009909106)
  expect_equal(unname(res4), valtest4, tolerance =1e-6)
})

test_that("fixef.merMod() should provide unscaled version at default", {
  res5 <- fixef.merMod(scf1)
  valtest5 <- fixef.merMod(scf2)
  expect_equal(res5, valtest5, tolerance =1e-10)
  
  res6 <- fixef.merMod(scf1, noScale = TRUE)
  valtest6 <- c(298.507891666666751, 1.621716753196637, 30.161580776828778) 
  expect_equal(unname(res6), valtest6, tolerance =1e-10)
})

test_that("vcov.merMod() should provide unscaled version at default", {
  res7 <- vcov.merMod(scf1)
  valtest7 <- vcov.merMod(scf2)
  expect_true(all.equal(unname(res7), unname(valtest7), tolerance =1e-6))
  
  res8 <- vcov.merMod(scf1, noScale = TRUE)
  valtest8 <- as(
    matrix(c(8.12278783645280e+01, -3.24525171154812e-15, 26.5884231100835215,
            -3.24525171154812e-15, 4.39230297872591e+00, 0.0344711252408276,
             2.65884231100835e+01, 3.44711252408276e-02, 19.7018459306493092), 
           nrow = 3, byrow = TRUE), "dpoMatrix")
  expect_true(all.equal(unname(res8), unname(valtest8), tolerance =1e-10))
})


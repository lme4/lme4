library("testthat")
library("lme4")

context("glmer.nb")
test_that("basic", {
   set.seed(101)
   dd <- expand.grid(f1 = factor(1:3),
                     f2 = LETTERS[1:2], g=1:9, rep=1:15,
                     KEEP.OUT.ATTRS=FALSE)
   mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
   dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)
   require("MASS")
   m.glm <- glm.nb(y ~ f1*f2, data=dd)
   library(lme4)
   m.nb <- glmer.nb(y ~ f1*f2 + (1|g), data=dd)
   expect_equal(fixef(m.nb),coef(m.glm),tol=1e-5)
})

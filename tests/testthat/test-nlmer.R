library("testthat")
library("lme4")

testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL")))
                 as.numeric(s) else 1

context("lower/upper bounds on nlmer models")
test_that("nlmer", {
    startvec <- c(Asym = 200, xmid = 725, scal = 350)
    nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
                 Orange, start = startvec,
                 control=nlmerControl(optCtrl=list(lower=c(0,200,-Inf,-Inf))))
    expect_equal(unname(fixef(nm1)[1]),200)
})

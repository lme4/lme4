library("testthat")
library("lme4")

context("Utilities (including *non*-exported ones")

test_that("namedList", {
    nList <- lme4:::namedList
    a <- b <- c <- 1
    expect_identical(nList(a,b,c),  list(a = 1, b = 1, c = 1))
    expect_identical(nList(a,b,d=c),list(a = 1, b = 1, d = 1))
    expect_identical(nList(a, d=pi, c), list(a = 1, d = pi, c = 1))
})

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

test_that("Var-Cov factor conversions", { ## from ../../R/vcconv.R
    mlist2vec <- lme4:::mlist2vec
    Cv_to_Vv <- lme4:::Cv_to_Vv
    Cv_to_Sv <- lme4:::Cv_to_Sv
    Sv_to_Cv <- lme4:::Sv_to_Cv
    Vv_to_Cv <- lme4:::Vv_to_Cv
    ##
    set.seed(1); cvec1 <- sample(10, 6)
    v1 <- Cv_to_Vv(cvec1)
    expect_equal(unname(v1), structure(c(9, 12, 15, 65, 34, 93), clen = 3))
    expect_equal(2, as.vector(Vv_to_Cv(Cv_to_Vv(2))))
    expect_equivalent(c(v1, 1),  Cv_to_Vv(cvec1, s=3) / 3^2)
    expect_equal(as.vector(ss1 <- Sv_to_Cv(Cv_to_Sv(cvec1))), cvec1)
    expect_equal(as.vector(vv1 <- Vv_to_Cv(Cv_to_Vv(cvec1))), cvec1)
    ## for length-1 matrices, Cv_to_Sv should be equivalent
    ##   to multiplying Cv by sigma and appending sigma ....
    clist2 <- list(matrix(1),matrix(2),matrix(3))
    cvec2 <- mlist2vec(clist2)
    expect_equal(cvec2, structure(1:3, clen = rep(1,3)), tolerance=0)
    expect_true(all((cvec3 <- Cv_to_Sv(cvec2, s=2)) == c(cvec2*2,2)))
    n3 <- length(cvec3)
    expect_equivalent(Sv_to_Cv(cvec3, n=rep(1,3), s=2), cvec3[-n3]/cvec3[n3])
})

test_that("nobar", {
    rr <- lme4:::RHSForm
    expect_equal(nobars(y~1+(1|g)),                      y~1)
    expect_equal(nobars(y~1|g),                          y~1)
    expect_equal(nobars(y~1+(1||g)),                     y~1)
    expect_equal(nobars(y~1||g),                         y~1)
    expect_equal(nobars(y~1+(x:z|g)),                    y~1)
    expect_equal(nobars(y~1+(x*z|g/h)),                  y~1)
    expect_equal(nobars(y~(1|g)+x+(x|h)),                y~x)
    expect_equal(nobars(y~(1|g)+x+(x+z|h)),              y~x)
    expect_equal(nobars(~1+(1|g)),                        ~1)
    expect_equal(nobars(~(1|g)),                          ~1)
    expect_equal(nobars(rr(y~1+(1|g))),                    1)
    expect_equal(nobars(rr(y~(1|g))),                      1)
})


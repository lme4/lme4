library("testthat")
library("lme4")

## use old (<=3.5.2) sample() algorithm if necessary
if ("sample.kind" %in% names(formals(RNGkind))) {
    suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
}

context("Utilities (including *non*-exported ones)")

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

## moved to lme4


test_that("getData", {
    ## test what happens when wrong version of 'data' is found in environment of formula ...
    f <- round(Reaction) ~ 1 + (1|Subject)
    g <- function() {
        data <- sleepstudy
        m1 <- glmer(f, data = data, family = poisson)
    }
    m1 <- g()
    expect_error(getData(m1), "object found is not a data frame or matrix")
    p1 <- suppressMessages(predict(m1, newparams = list(beta = 5, theta = 1), type = "response",
            re.form = ~ 1|Subject))
    expect_equal(unname(head(p1, 2)), c(444.407382298401, 444.407382298401))
})

library("testthat")

data(Dyestuff, package="lme4")
n     <- nrow(Dyestuff)
ones  <- rep.int(1, n)
zeros <- rep.int(0, n)
YY    <- Dyestuff$Yield
mYY   <- mean(YY)

context("lmerResp objects")
test_that("lmerResp", {
    mres  <- YY - mYY
    rr    <- lmerResp$new(y=YY)

    expect_that(rr$weights,                   equals(ones))
    expect_that(rr$sqrtrwt,                   equals(ones))
    expect_that(rr$sqrtXwt,                   equals(ones))
    expect_that(rr$offset,                    equals(zeros))
    expect_that(rr$mu,                        equals(zeros))
    expect_that(rr$wtres,                     equals(YY))
    expect_that(rr$wrss(),                    equals(sum(YY^2)))
    expect_that(rr$updateMu(rep.int(mYY, n)), equals(sum(mres^2)))
    expect_that(rr$REML,                      equals(0L))

    rr$REML <- 1L
    expect_that(rr$REML,                      equals(1L))
})

mlYY <- mean(log(YY))
gmeanYY <- exp(mlYY)                    # geometric mean

context("glmResp objects")
test_that("glmResp", {
    mres  <- YY - gmeanYY
    gmean <- rep.int(gmeanYY, n)
    rr    <- glmResp$new(family=poisson(), y=YY)
    
    expect_that(rr$weights,                    equals(ones))
    expect_that(rr$sqrtrwt,                    equals(ones))
    expect_that(rr$sqrtXwt,                    equals(ones))
    expect_that(rr$offset,                     equals(zeros))
    expect_that(rr$mu,                         equals(zeros))
    expect_that(rr$wtres,                      equals(YY))
    expect_that(rr$n,                          equals(ones))

    ## wrss() causes an update of mu which becomes ones, wtres also changes
    expect_that(rr$wrss(),                     equals(sum((YY-1)^2)))
    expect_that(rr$mu,                         equals(ones))
    expect_that(rr$wtres,                      equals(YY-ones))

    expect_that(rr$updateMu(rep.int(mlYY, n)), equals(sum(mres^2)))
    expect_that(rr$mu,                         equals(gmean))
    expect_that(rr$muEta(),                    equals(gmean))
    expect_that(rr$variance(),                 equals(gmean))

    rr$updateWts()
    expect_that(1/sqrt(rr$variance()),         equals(rr$sqrtrwt))
    expect_that(as.vector(rr$sqrtXwt),         equals(rr$sqrtrwt * rr$muEta()))
})


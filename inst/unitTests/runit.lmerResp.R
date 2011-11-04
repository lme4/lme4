set.seed(1)

test.lmerResp <- function() {
    n <- nrow(Dyestuff)
    YY <- Dyestuff$Yield
    mYY <- mean(YY)
    mres <- YY - mYY
    rr <- lmerResp$new(y=YY)
    checkEquals(rr$weights, rep.int(1, n))
    checkEquals(rr$sqrtrwt(), rep.int(1, n))
    checkEquals(rr$sqrtXwt(), array(rep.int(1, n), c(n, 1L)))
    checkEquals(rr$offset, rep.int(0, n))
    checkEquals(rr$fitted(), rep.int(0, n))
    checkEquals(rr$wtres(), YY)
    checkEquals(rr$wrss(), sum(YY * YY))
    checkEquals(rr$updateMu(rep.int(mYY, n)), sum(mres^2))
    checkEquals(rr$reml, 0L)
    rr$reml <- 1L
    checkEquals(rr$reml, 1L)
}

test.glmResp <- function() {
    n <- nrow(Dyestuff)
    YY <- Dyestuff$Yield
    mlYY <- mean(log(YY))
    gmeanYY <- exp(mlYY)                # geometric mean
    mres <- YY - gmeanYY
    rr <- glmResp$new(family=poisson(), y=YY)

    checkEquals(rr$weights, rep.int(1, n))
    checkEquals(rr$sqrtrwt(), rep.int(1, n))
    checkEquals(rr$sqrtXwt(), array(rep.int(1, n), c(n, 1L)))
    checkEquals(rr$offset, rep.int(0, n))
    checkEquals(rr$fitted(), rep.int(0, n))
    checkEquals(rr$wtres(), YY)
    checkEquals(rr$wrss(), sum(YY^2))
    checkEquals(rr$n, rep.int(1, n))
    checkEquals(rr$updateMu(rep.int(mlYY, n)), sum(mres^2))
    checkEquals(rr$fitted(), rep.int(gmeanYY, n))
    checkEquals(rr$muEta(), rep.int(gmeanYY, n))
    checkEquals(rr$variance(), rep.int(gmeanYY, n))
    rr$updateWts()
    checkEquals(1/sqrt(rr$variance()), rr$sqrtrwt())
    checkEquals(as.vector(rr$sqrtXwt()), rr$sqrtrwt() * rr$muEta())
}

set.seed(1)


test.lmerResp.twoarg <- function() {
    n <- nrow(Dyestuff)
    YY <- Dyestuff$Yield
    mYY <- mean(YY)
    mres <- YY - mYY
    rr <- lmerResp$new(1L, YY)
    checkEquals(rr$weights, rep.int(1, n))
    checkEquals(rr$sqrtrwt, rep.int(1, n))
    checkEquals(rr$sqrtXwt, rep.int(1, n))
    checkEquals(rr$offset,  rep.int(0, n))
    checkEquals(rr$mu,      rep.int(0, n))
    checkEquals(rr$wtres,   YY)
    checkEquals(rr$updateMu(rep.int(mYY, n)), sum(mres^2))
}

test.lmerResp.threearg <- function() {
    n <- nrow(Dyestuff)
    YY <- Dyestuff$Yield
    WW <- rep(1, n)
    ZZ <- rep(0, n)
    sqrtWW <- sqrt(WW)
    mYY <- mean(YY)
    mres <- YY - mYY
    rr <- lmerResp$new(1L, YY, WW)
    checkEquals(rr$weights, WW)
    checkEquals(rr$sqrtrwt, WW)
    checkEquals(rr$sqrtXwt, WW)
    checkEquals(rr$offset,  ZZ)
    checkEquals(rr$mu,      ZZ)
    checkEquals(rr$wtres,   YY)
    checkEquals(rr$updateMu(rep.int(mYY, n)), sum(mres^2))
}

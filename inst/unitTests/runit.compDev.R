test.compDev <- function() {            # cross-check compDev versus R evaluation
    fm1 <- lmer(Yield ~ 1|Batch, Dyestuff, doFit = FALSE)
    dd1c <- mkdevfun(fm1)
    dd1u <- mkdevfun(fm1, compDev = FALSE)
    checkEquals(dd1c(1), dd1u(1))
    checkEquals(bobyqa(1,dd1c,0), bobyqa(1,dd1u,0))
}

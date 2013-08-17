## various tests of random-effects identifiability:

## 1. from inst/tests/test-lmer.R:
## only two observations per individual, attempt to fit random slopes model
## correctly identified
library(testthat)
library(lme4)
lFormula(Reaction ~ 1 + Days + (1 + Days | Subject),
                      data = sleepstudy, subset = (Days == 1 | Days == 9))

## 2. from tests/lmer.R:
## false positive
lsDat <- data.frame(
    Operator = as.factor(rep(1:5, c(3,4,8,8,8))),
    Part = as.factor(
        c(2L, 3L, 5L,
          1L, 1L, 2L, 3L,
          1L, 1L, 2L, 2L, 3L, 3L, 4L, 5L,
          1L, 2L, 3L, 3L, 4L, 4L, 5L, 5L,
          1L, 2L, 2L, 3L, 3L, 4L, 5L, 5L)),
    y =
    c(0.34, -1.23, -2.46,
      -0.84, -1.57,-0.31, -0.18,
      -0.94, -0.81, 0.77, 0.4, -2.37, -2.78, 1.29, -0.95,
      -1.58, -2.06, -3.11,-3.2, -0.1, -0.49,-2.02, -0.75,
      1.71,  -0.85, -1.19, 0.13, 1.35, 1.92, 1.04,  1.08))
xtabs( ~ Operator + Part, data=lsDat) 
lf <- lFormula(y ~ (1|Part) + (1|Operator) + (1|Part:Operator), data = lsDat,
               control=lmerControl(check.nobs.vs.rankZ="ignore"))

## right now this is what checkRankZ does (without the 'tol' optional
## argument)
rfun <- function(Zt,tol) {
    rankMatrix(if (doTr) 
               t(Zt)
    else Zt, method = "qr", sval = numeric(min(d)),
               tol=tol)
}
image(lf$reTrms$Zt)
## exploring: drastically

tol <- 1e-3
Zt <- lf$reTrms$Zt
c(rankMatrix(Zt)) ## 21
c(rankMatrix(Zt,tol=tol)) ## 21
c(rankMatrix(Zt,method="qr",tol=tol)) ## 31
c(rankMatrix(t(Zt),method="qr",tol=tol)) ## 30

## 3. false negative: two identical observation-level random effects
cbppX <- transform(cbpp,obs=seq(nrow(cbpp)),obs2=seq(nrow(cbpp)))
gf <- glFormula(cbind(incidence, size - incidence) ~ period + (1 | herd) +
                 (1|obs) + (1|obs2),
                 data = cbppX, family = binomial)

image(gf$reTrms$Zt)

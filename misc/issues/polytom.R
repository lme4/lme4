## from
## library(polytomous)
## data(think)
## polytomous(Lexeme ~ Agent + Patient + (1|Register),
##      data=think, heuristic="poisson.reformulation")

load(system.file("testdata","polytom2.RData",package="lme4"))
library(lme4)
## library(lme4.0)
g1 <- glmer(formula=formula.poisson,data=data.poisson,family=poisson)
## from lme4.0:
## does work, but we get fixed effect parameters with abs > 18 ...
load(system.file("testdata","polytom3.RData",package="lme4"))
g2 <- glmer(formula=formula.poisson,data=data.poisson,family=poisson)
## Error: PIRLS step-halving failed to reduce deviance in pwrssUpdate
## In addition: Warning messages:
## 1: In pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, verbose) :
##   Cholmod warning 'not positive definite' at file:../Cholesky/t_cholmod_rowfac.c, line 431
## 2: In pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, verbose) :
##  Cholmod warning 'not positive definite' at file:../Cholesky/t_cholmod_rowfac.c, line 431



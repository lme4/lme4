## from
if(FALSE) {
    library(polytomous)
    data(think)
    system.time(
        fm.think <-
            polytomous(Lexeme ~ Agent + Patient + (1|Register),
                       data=think, heuristic="poisson.reformulation")
    ) ## 21 sec [lynne, 2014-11]
    ## --> 5 warnings
}

attach(system.file("testdata", "polytom2.RData", package="lme4"))
ls.str(2)
## formula.poisson, data.poisson

library(lme4)
## library(lme4.0)
system.time( g1 <- glmer(formula.poisson, data=data.poisson, family=poisson) )
## from lme4.0:
## does work, but we get fixed effect parameters with abs > 18 ...

attach(system.file("testdata", "polytom3.RData", package="lme4"))

system.time( g2 <- glmer(formula.poisson, data=data.poisson, family=poisson) )

## Error: PIRLS step-halving failed to reduce deviance in pwrssUpdate
## In addition: Warning messages:
## 1: In pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, verbose) :
##   Cholmod warning 'not positive definite' at file:../Cholesky/t_cholmod_rowfac.c, line 431
## 2: In pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, verbose) :
##  Cholmod warning 'not positive definite' at file:../Cholesky/t_cholmod_rowfac.c, line 431



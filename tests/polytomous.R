if (FALSE) {
    ## library(polytomous)
    data(think)
    think.polytomous.lmer1 <-  polytomous(Lexeme ~ Agent + Patient + (1|Register),
                                          data=think, heuristic="poisson.reformulation")
    save("formula.poisson","data.poisson",file="polytomous_test.RData")
}
load("polytomous_test.RData")
library(lme4)
if (FALSE) {
    ## infinite loop
    glmer(formula.poisson,data=data.poisson,family=poisson,verbose=10)
    ## Cholmod not positive definite -> infinite loop
    glmer(formula.poisson,data=data.poisson,family=poisson,verbose=10,optimizer="bobyqa")
    ## caught warning: maxfun < 10 * length(par)^2 is not recommended. -> infinite loop
}
## works but sloooow ....
if (FALSE) {
    try(g1 <- glmer(formula.poisson,data=data.poisson,family=poisson,
                    compDev=FALSE,verbose=1))
    ## runs for 2880 steps until:
    ## Error in pp$updateDecomp() : Downdated VtV is not positive definite
}
glmer(formula.poisson,data=data.poisson,family=poisson,
      compDev=FALSE,verbose=1,optimizer="bobyqa")
## caught warning: maxfun < 10 * length(par)^2 is not recommended.
## but runs to completion





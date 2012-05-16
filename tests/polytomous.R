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




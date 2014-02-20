library(lme4)

 ## setup
 ## library(polytomous)
 ##   data(think)
 ##   think.polytomous.lmer1 <-  polytomous(Lexeme ~ Agent + Patient + (1|Register),
 ##                                         data=think, heuristic="poisson.reformulation")
 ##   save("formula.poisson","data.poisson",file="polytomous_test.RData")

load(system.file("testdata","polytomous_test.RData",package="lme4"))

if (FALSE) {
    ## infinite loop
    glmer(formula.poisson,data=data.poisson,family=poisson,verbose=10)
    ## Cholmod not positive definite -> infinite loop
    glmer(formula.poisson,data=data.poisson,family=poisson,
          verbose=10,control=glmerControl(optimizer="bobyqa"))
    ## caught warning: maxfun < 10 * length(par)^2 is not recommended. -> infinite loop
}
## works but sloooow ....
if (FALSE) {
    try(g1 <- glmer(formula.poisson,data=data.poisson,family=poisson,
                    control=glmerControl(compDev=FALSE),verbose=1))
    ## runs for 2880 steps until:
    ## Error in pp$updateDecomp() : Downdated VtV is not positive definite
}

(testLevel <- lme4:::testLevel())
if (testLevel > 2) {
    glmer(formula.poisson,data=data.poisson,family=poisson,
          control=glmerControl(compDev=FALSE),optimizer="bobyqa")
    ## caught warning: maxfun < 10 * length(par)^2 is not recommended.
    ## but runs to completion
}





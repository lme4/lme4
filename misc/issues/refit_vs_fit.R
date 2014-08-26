## from Github #231

library(lme4)
full_mod1 <- glmer(TICKS ~ YEAR+scale(HEIGHT)+(1|LOCATION/BROOD/INDEX),
                   family="poisson", data=grouseticks)

set.seed(12345)
debug(.simulateFun)
simTICKS.tab <- simulate(full_mod1, nsim=2, seed=12345)
## ??? things seem to be much slower than I expected

r1 <- refit(full_mod1,simTICKS.tab[[1]])
r2 <- refit(full_mod1,simTICKS.tab[[2]],verbose=100) ## fails
debug(lme4:::refit.merMod)
load("devfunUp.RData")  ## contains ff, updated devfun from inside refit
ff2 <- update(full_mod1,devFunOnly=TRUE)
fitpars <- unlist(getME(full_mod1,c("theta","fixef")))
ff(fitpars)
ls(e1 <- environment(ff))
ls(e2 <- environment(ff2))
for (i in ls(environment(ff))) {
    print(i)
    print(all.equal(e1[[i]],e2[[i]]))
}
ls(e1)
ls(e2)

## next step: compare environments carefully, try to make them match
## as much as possible

r1A <- update(full_mod1,data=transform(grouseticks,TICKS=simTICKS.tab[[1]]))
r2A <- update(full_mod1,data=transform(grouseticks,TICKS=simTICKS.tab[[2]]),
              verbose=100)
## fit succeeds

r2B <- update(full_mod1,data=transform(grouseticks,TICKS=simTICKS.tab[[2]]),
              verbose=100,start=getME(full_mod1,c("theta","fixef")))
## skips nAGQ=0 stage

all.equal(fixef(r1),fixef(r1A),tol=2e-5)
all.equal(unlist(VarCorr(r1)),unlist(VarCorr(r1A)),tol=2e-5)

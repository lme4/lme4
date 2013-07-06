library(lme4)
(testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1)

## testing refit
## for each type of model, should be able to
##  (1) refit with same data and get the same answer,
##     at least structurally (small numerical differences
##     are probably unavoidable)
##  (2) refit with simulate()d data

getinfo <- function(x) {
  c(fixef(x),logLik(x))
}

## LMM
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm1R <- refit(fm1,sleepstudy$Reaction)
fm1S <- refit(fm1,simulate(fm1)[[1]])
getinfo(fm1)
getinfo(fm1R)
getinfo(fm1S)
stopifnot(all.equal(fm1,fm1R,tol=1.5e-5))
## sapply(slotNames(fm1),
##        function(x) isTRUE(all.equal(slot(fm1,x),slot(fm1R,x),tol=1.5e-5)))
## fm1@optinfo
## fm1R@optinfo

getinfo(refitML(fm1))

## binomial GLMM (two-column)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            cbpp, binomial)
gm1R <- refit(gm1,with(cbpp,cbind(incidence,size-incidence)))

## FIXME: testing all-zero responses
## this gives "pwrssUpdate did not converge in 30 iterations"
## not sure if it's pathological or not
if (FALSE) {
 sim1Z <- simulate(gm1)[[1]]
 sim1Z[4,] <- c(0,0)
 refit(gm1,sim1Z)
}

## FIXME: why is the difference this large?
stopifnot(all.equal(gm1,gm1R,tol=6e-4))
getinfo(gm1)
getinfo(gm1R)

## FIXME: still failing on Windows
## gm1S <- refit(gm1,simulate(gm1)[[1]])
## getinfo(gm1S)

## binomial GLMM (prob/weights)
gm2 <- glmer(incidence/size ~ period + (1 | herd), cbpp, binomial, weights=size)
gm2R <- refit(gm2,with(cbpp,incidence/size))
stopifnot(all.equal(gm2,gm2R,tol=6e-4))
getinfo(gm2)
getinfo(gm2R)
## FIXME: check on Windows
## gm2S <- refit(gm2,simulate(gm2)[[1]])
## getinfo(gm2S)

if (testLevel > 1) {
    ## Bernoulli GLMM (specified as factor)
    data(Contraception,package="mlmRev")
    gm3 <- glmer(use ~ urban+age+livch+(1|district), Contraception, binomial)
    gm3R <- refit(gm3,Contraception$use)
    gm3S <- refit(gm3,simulate(gm3)[[1]])
    stopifnot(all.equal(gm3,gm3R,tol=3e-4))
    getinfo(gm3)
    getinfo(gm3R)
    getinfo(gm3S)

    data(Mmmec,package="mlmRev")
    gm4 <- glmer(deaths ~ uvb + (1|region), data=Mmmec,
                 family=poisson,
                 offset = log(expected))
    gm4R <- refit(gm4,Mmmec$deaths)
    gm4S <- refit(gm4,simulate(gm4)[[1]])
    getinfo(gm4)
    getinfo(gm4R)
    getinfo(gm4S)

    stopifnot(all.equal(gm4,gm4R,tol=6e-5))
}


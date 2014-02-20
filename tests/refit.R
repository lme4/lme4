library(lme4)
library(testthat)

load(system.file("testdata","lme-tst-fits.rda",package="lme4"))

if (getRversion()>"3.0.0") {
    ## saved fits are not safe with old R versions

    ## testing refit
## for each type of model, should be able to
##  (1) refit with same data and get the same answer,
##     at least structurally (small numerical differences
##     are probably unavoidable)
##  (2) refit with simulate()d data

getinfo <- function(x) {
  c(fixef(x),logLik(x),unlist(ranef(x)),unlist(VarCorr(x)))
}

dropterms <- function(x) {
    attr(x@frame,"terms") <- NULL
    x
}

## LMM
fm1 <- fit_sleepstudy_2
fm1R <- refit(fm1,sleepstudy$Reaction)
fm1S <- refit(fm1,simulate(fm1)[[1]])

stopifnot(all.equal(getinfo(fm1),getinfo(fm1R),tolerance=3e-5))
## sapply(slotNames(fm1),
##        function(x) isTRUE(all.equal(slot(fm1,x),slot(fm1R,x),tolerance=1.5e-5)))
## fm1@optinfo
## fm1R@optinfo

getinfo(refitML(fm1))

## binomial GLMM (two-column)
gm1 <- fit_cbpp_1
gm1R <- refit(gm1, with(cbpp,cbind(incidence,size-incidence)))

## FIXME: testing all-zero responses
## this gives "pwrssUpdate did not converge in 30 iterations"
## not sure if it's pathological or not
if (FALSE) {
 sim1Z <- simulate(gm1)[[1]]
 sim1Z[4,] <- c(0,0)
 refit(gm1,sim1Z) # Error:  ... PIRLS ... failed ...
}

## FIXME: why is the difference this large?
## check components ...
stopifnot(all.equal(getinfo(gm1),getinfo(gm1R),tolerance=1e-4))

if(FALSE) {## FIXME: still failing on Windows
gm1S <- refit(gm1,simulate(gm1)[[1]])
getinfo(gm1S)
}

## binomial GLMM (prob/weights)
gm2 <- fit_cbpp_3
## glmer(incidence/size ~ period + (1 | herd), cbpp, binomial, weights=size)
gm2R <- refit(gm2,with(cbpp,incidence/size))
stopifnot(all.equal(getinfo(gm2),getinfo(gm2R),tolerance=6e-4))

## from Alexandra Kuznetsova
set.seed(101)
Y <- matrix(rnorm(1000),ncol=2)
d <- data.frame(y1=Y[,1],  x=rnorm(100), f=rep(1:10,10))
fit1 <- lmer(y1 ~ x+(1|f),data=d)
fit2 <- refit(fit1, newresp = Y[,2], rename.response=TRUE)
## check, but ignore terms attribute of model frame ...
expect_warning(refit(fit1, newresp = Y[,2], junk=TRUE))
if (isTRUE(all.equal(fit1,fit2))) stop("fit1 and fit2 should not be equal")
## hack number of function evaluations
u2 <- update(fit2)
fit2@optinfo$feval <- u2@optinfo$feval <-  NA
stopifnot(all.equal(d1 <- dropterms(fit2), d2 <- dropterms(u2)))

## for (i in slotNames(d1)) {
##     cat(i,"\n")
##    print(isTRUE(all.equal(slot(d1,i),slot(d2,i))))
## }
## FIXME: check on Windows
## gm2S <- refit(gm2,simulate(gm2)[[1]])
## getinfo(gm2S)

if (lme4:::testLevel() > 1) {
    ## Bernoulli GLMM (specified as factor)
    if (require("mlmRev")) {
        data(Contraception,package="mlmRev")
        gm3 <- glmer(use ~ urban+age+livch+(1|district),
                     Contraception, binomial)
        gm3R <- refit(gm3,Contraception$use)
        gm3S <- refit(gm3,simulate(gm3)[[1]])
        stopifnot(all.equal(getinfo(gm3),getinfo(gm3R),tolerance=3e-4))
        print(getinfo(gm3))
        print(getinfo(gm3R))
        print(getinfo(gm3S))

        data(Mmmec,package="mlmRev")
        gm4 <- glmer(deaths ~ uvb + (1|region), data=Mmmec,
                     family = poisson,
                     offset = log(expected))
        gm4R <- refit(gm4,Mmmec$deaths)
        if(FALSE) ## FIXME: following fails, not finding 'expected'
        gm4S <- refit(gm4,simulate(gm4)[[1]])
        print( getinfo(gm4) )
        print( getinfo(gm4R) )
        if(FALSE) ## FIXME (above)
        getinfo(gm4S)

        stopifnot(all.equal(getinfo(gm4),getinfo(gm4R),tolerance=6e-5))
    }
}

}

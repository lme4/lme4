## polytomous package example breaks with current release of lme4
## it's arguably a very dodgy example, but it is still worth trying
## to diagnose/investigate what's going on here ...
##
## setup code -- results are saved in a file, so you don't
## need to do this stuff.
if (FALSE) {
    library(polytomous)
    think.polytomous.lmer1 <-
        polytomous(Lexeme ~ Agent + Patient + (1|Register),
                   data=think, heuristic="poisson.reformulation")
    ## extract V, object from browser mode
    form <- formula(object)
    data.poisson <- model.frame(object)
    library(lme4.0)
    ## lme4.0 fit:
    m0 <- glmer(form,data=data.poisson,family=poisson)
    ## lme4 fit (vcov.merMod fails):
    m1 <- object
    save("m0","m1","V","form","data.poisson",file="polytomous_test.RData")
}
## start here ...
load("polytomous_test.RData")
## load both old & new versions so we can make comparisons
library(lme4)
library(lme4.0)
## full parameter set for lme4.0 and lme4 fits
m0parms <- c(getME(m0,"theta"),fixef(m0))
m1parms <- c(lme4::getME(m1,"theta"),fixef(m1))
length(m0parms)     ## 50 parameters
nrow(data.poisson)  ## 220 observations
## Number of obs: 220, groups: Observation, 55; Register, 2
## This is likely to be overfitted (N/(# fixed coefs) approx. 4,
## 10 is a sensible minimum ...)
## Also: should probably *not* be treating register (factor with only
## two levels) as a RE, at least not without regularization ...
## Compare parameter sets graphically -- not nearly identical, but
## at least vaguely in the same ballpark
## Some parameters <<0, in both old & new fits
##    -- suggests some categories are essentially all
##    zero counts (again, need for regularization)
par(las=1,bty="l")
plot(m0parms,m1parms,col=rep(2:1,c(2,48)))
abline(a=0,b=1)

## test parameter sets in 'new' deviance function
dd <- lme4::glmer(form,data.poisson,family=poisson,devFunOnly=TRUE)
dd(m0parms)
dd(m1parms)  ## m1parms are supposedly MUCH better ...

## but we have a problem with the variance-covariance matrix
try(vcov(m1))
## Error in vcov.merMod(m1) : 
##   Computed variance-covariance matrix problem: not a positive definite matrix
try(as(V,"dpoMatrix"))  ## ditto

## look for asymmetries
asymm <- function(x) {
    u <- x[upper.tri(x)]
    l <- t(x)[upper.tri(x)]
    sdiff <- abs(u-l)
    if (max(sdiff)==0) return(FALSE)
    c(max=max(sdiff),pos=which(x==u[which.max(sdiff)],arr.ind=TRUE))
}
image(V) ## doesn't look symmetric???  Is this matrix somehow
         ## confusing the image method for "dgeMatrix" objects?
V[47,47]
asymm(V)  ## symmetric
head(rev(e <- eigen(V)$values))  ## ... but one disgusting eigenvalue
range(e)
## https://stat.ethz.ch/pipermail/r-help/2004-May/051054.html
## trace("coerce",browser,sig=c("dgeMatrix","dpoMatrix"))
as(V,"dpoMatrix")
## guts of coerce method
as(V3 <- as(V2 <- as(V1 <- as(V, "symmetricMatrix"), "dMatrix"),
            "denseMatrix"), 
   "dpoMatrix")
as(V3,"dpoMatrix")
## trace("coerce",browser,sig=c("dsyMatrix","dpoMatrix"))
## again, extract the guts
.Call(Matrix:::dpoMatrix_chol,V3)
## the leading minor of order 47 is not pos def
## we *could* fall back on the RX-based computation:
image(vcov.merMod(m1,use.hessian=FALSE))


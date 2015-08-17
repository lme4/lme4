library("lme4")

set.seed(101)
dd <- expand.grid(f1 = factor(1:3),
                  f2 = LETTERS[1:2], g=1:9, rep=1:15,
          KEEP.OUT.ATTRS=FALSE)
mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

## mimic glmer.nb protocol

## basic Poisson fit
m.base <- glmer(y ~ f1*f2 + (1|g), data=dd, family=poisson)
all.equal(m.base@beta,(m.base.r <- refit(m.base))@beta)  ## OK

th <- lme4:::est_theta(m.base,limit=20,eps=1e-4,trace=FALSE)
th0 <- structure(0.482681268108477, SE = 0.0244825021248148)
th1 <- structure(0.482681277470945)
th2 <- 0.482681268108477
th3 <- 0.4826813
## NB update with raw number
m.numth1 <- update(m.base,family=negative.binomial(theta=0.4826813))
all.equal(m.numth1@beta,(m.numth1.r <- refit(m.numth1))@beta)  ## OK

## strip NB value
m.symth4 <- update(m.base,family=negative.binomial(theta=c(th)))
all.equal(m.symth4@beta,(m.symth4.r <- refit(m.symth4))@beta)  ## fails

## IDENTICAL numeric value to case #1 above
m.symth6 <- update(m.base,family=negative.binomial(theta=th3))
all.equal(m.symth6@beta,(m.symth6.r <- refit(m.symth6))@beta)  ## works

## standard NB update with computed theta from est_theta (incl SE attribute)
m.symth <- update(m.base,family=negative.binomial(theta=th))
all.equal(m.symth@beta,(m.symth.r <- refit(m.symth))@beta)  ## not OK

## NB update with equivalent value
m.symth2 <- update(m.base,family=negative.binomial(theta=th0))
all.equal(m.symth2@beta,(m.symth2.r <- refit(m.symth2))@beta)  ##

## NB update with theta value (stored as variable, no SE) only
##  (off by -1e-8)
m.symth3 <- update(m.base,family=negative.binomial(theta=th1))
all.equal(m.symth3@beta,(m.symth3.r <- refit(m.symth3))@beta)  ## works

## strip NB value (off by 5e-16)
m.symth5 <- update(m.base,family=negative.binomial(theta=th2))
all.equal(m.symth5@beta,(m.symth5.r <- refit(m.symth5))@beta)  ## fails

## Surprised by difference between specifying as variable vs.
##  number; could it be a deep-copying problem?

## if (FALSE) {  ## slow
    mfun <- function(x) isTRUE(all.equal(x,m.symth@beta))
    rmat <- matrix(ncol=10,nrow=length(fixef(m.symth)))
    for (i in 1:10) {
        cat("*",i)
        m.numthXX <- update(m.base,family=negative.binomial(theta=0.4826813+i*1e-10))
        tt <- system.time(rmat[,i] <- refit(m.numthXX)@beta)
        cat(" ",mfun(rmat[,i]),tt["elapsed"],"\n")
    }

    ## results highly unstable:
    ##
    apply(rmat,2,mfun)
    ##  Ubuntu 14.04 64-bit:
    ##  [1]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE  TRUE FALSE FALSE
    ##  Ubuntu 12.04 32-bit:
    ##  [1]  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE
    ##
    ## bad fits take longer (>30 seconds vs ~4 seconds on slow machine)
    ## pattern is repeatable (so probably not memory corruption??)
    ## numeric instability is surprising; does it happen on all platforms
    ## but in slightly different places on each?

## }
## try with different optimizer (red herring??)
## 
##  m4 <- update(m1,family=negative.binomial(theta=th),
##              control=glmerControl(optimizer="bobyqa"))
## all.equal(m4@beta,(m4B <- refit(m4))@beta)

refit.all <- ls(pattern="^m\\..*\\.r$")
refit.all <- mget(refit.all)
(rtab <- do.call(cbind,lapply(refit.all,fixef)))

sessionInfo()

#### diagnosing:
cmpFields <- function(x,y,comp=identical,...) {
    if (is.list(x)) {
        fields <- names(x)
    } else {
        fields <- names(x$getRefClass()$fields())
    } ## also need one for S4: have I done this already somewhere?
    res <- setNames(logical(length(fields)),fields)
    for (i in fields) {
        res[i] <- comp(x[[i]],y[[i]],...)
    }
    return(res)
} 
    
rr_i3 <- readRDS("respmod_i3.rds")
rr_i4 <- readRDS("respmod_i4.rds")
cmpFields(rr_i3,rr_i4)  ## only family is FALSE
all.equal(rr_i3$family,rr_i4$family)
cmpFields(rr_i3$family,rr_i4$family)
## variance, dev.resids, aic, validmu, simulate have different envs
cmpFields(rr_i3$family,rr_i4$family,ignore.environment=TRUE)

## don't think we've hit the problem yet
getRDA <- function(f,suffix) {
    L <- load(f)
    for (i in L) 
        assign(paste(i,suffix,sep="_"),
               get(i),.GlobalEnv)
}
newL <- new("dgCMatrix",i=pp$Lambdat@i,p=pp$Lambdat@p,
    x = rep(0,length(pp$Lambdat@x)), Dim=c(9L,9L))
pp$Lambdat <- newL
getRDA("after_GPU_i3.rda","i3")
getRDA("after_GPU_i4.rda","i4")
cmpFields(rr_i3,rr_i4)
cmpFields(rr_i3,rr_i4,comp=all.equal)
cmpFields(pp_i3,pp_i4,comp=all.equal)
cmpFields(pp_i3,pp_i4,comp=all.equal,tol=1e-3)
pp_i3$Lambdat  ## all zero -- does this collapse to nothing somewhere??
pp_i4$Lambdat  ## tiny but not zero
## would substituting one Lambdat for the other make a difference?
## but ... i3 works, i4 *doesn't*.  Should Lambdat, LamtUt really be all-zero?

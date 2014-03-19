## In both of these cases boundary fit (i.e. estimate of zero RE
## variance) is *incorrect*. (Nelder_Mead, restart_edge=FALSE) is the
## only case where we get stuck; either optimizer=bobyqa or
## restart_edge=TRUE (default) works

library(lme4)
library(testthat)

if(!dev.interactive(orNone=TRUE)) pdf("boundary_plots.pdf")

## Stephane Laurent:
dat <- read.csv(system.file("testdata","dat20101314.csv",package="lme4"))

fit <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=lmerControl(optimizer="Nelder_Mead"))
fit_b <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=lmerControl(optimizer="bobyqa",restart_edge=FALSE))
fit_c <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=lmerControl(optimizer="Nelder_Mead", restart_edge=FALSE,
              check.conv.hess="ignore"))
## final fit gives degenerate-Hessian warning
## FIXME: use fit_c with expect_warning() as a check on convergence tests
## tolerance=1e-5 seems OK in interactive use but not in R CMD check ... ??
stopifnot(all.equal(getME(fit,"theta"),getME(fit_b,"theta"),tolerance=2e-5))
stopifnot(all(getME(fit,"theta")>0))

## Manuel Koller

source(system.file("testdata","koller-data.R",package="lme4"))
ldata <- getData(13)
## old (backward compatible/buggy)
fm4  <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="Nelder_Mead",
                                  use.last.params=TRUE))

fm4b <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="Nelder_Mead",
                                  use.last.params=TRUE, restart_edge=FALSE,
                                  check.conv.hess="ignore",
                                  check.conv.grad="ignore"))
## FIXME: use as convergence test check
stopifnot(getME(fm4b,"theta")==0)
fm4c <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="bobyqa",
                                  use.last.params=TRUE))
stopifnot(all.equal(getME(fm4,"theta"),getME(fm4c,"theta"),tolerance=1e-4))
stopifnot(getME(fm4,"theta") > 0)

## new: doesn't get stuck at edge any more,  but gets stuck somewhere else ...
fm5 <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="Nelder_Mead",
				  check.conv.hess="ignore",
				  check.conv.grad="ignore"))
fm5b <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="Nelder_Mead",
                                  restart_edge=FALSE,
                                  check.conv.hess="ignore",
                                  check.conv.grad="ignore"))
fm5c <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="bobyqa"))
stopifnot(all.equal(unname(getME(fm5c,"theta")), 0.21067645, tolerance = 1e-7))
					#	 0.21067644264 [64-bit, lynne]

if (FALSE) {
    ## additional stuff for diagnosing Nelder-Mead problems.
    ## library(nloptr) call commented out to avoid R CMD check problems/needing to
    ##  Suggest: nloptr

    library(optimx)
    fm5d <- update(fm5,control=lmerControl(optimizer="optimx",
                       optCtrl=list(method="L-BFGS-B")))

    ## library(nloptr)
    defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-6,maxeval=1e5)
    nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
        for (n in names(defaultControl))
            if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
        res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
        with(res,list(par=solution,
                      fval=objective,
                      conv=if (status>0) 0 else status,
                      message=message))
    }
    fm5e <- update(fm5,control=lmerControl(optimizer="nloptwrap2"))

    mList <- setNames(list(fm4,fm4b,fm4c,fm5,fm5b,fm5c,fm5d,fm5e),
                      c("NM/uselast","NM/uselast/norestart","bobyqa/uselast",
                        "NM","NM/norestart","bobyqa","LBFGSB","nloptr/bobyqa"))
    pp <- profile(fm5c,which=1)
    dd <- as.data.frame(pp)
    par(las=1,bty="l")
    v <- sapply(mList,
                function(x) sqrt(VarCorr(x)[[1]]))
    (v2 <- sapply(mList,getME,"theta"))

    plot(.zeta^2~.sig01,data=dd,type="b")
    abline(v=v)
}

######################
library(lattice)
## testing boundary and near-boundary cases

tmpf <- function(i,...) {
    set.seed(i)
    d <- data.frame(x=rnorm(60),f=factor(rep(1:6,each=10)))
    d$y <- simulate(~x+(1|f),family=gaussian,newdata=d,
                    newparams=list(theta=0.01,beta=c(1,1),sigma=5))[[1]]
    lmer(y~x+(1|f),data=d,...)
}
sumf <- function(m) {
    unlist(VarCorr(m))[1]
}
if (FALSE) {
    ## figuring out which seeds will give boundary and
    ## near-boundary solutions
    mList <- lapply(1:201,tmpf)
    ss <- sapply(mList,sumf)+1e-50
    par(las=1,bty="l")
    hist(log(ss),col="gray",breaks=50)
    ## values lying on boundary
    which(log(ss)<(-40))   ## 5, 7-13, 15, 21, ...
    ## values close to boundary (if check.edge not set)
    which(log(ss)>(-40) & log(ss) <(-20))  ## 16, 44, 80, 86, 116, ...
}
## diagnostic plot
tmpplot <- function(i,FUN=tmpf) {
    dd <- FUN(i,devFunOnly=TRUE)
    x <- 10^seq(-10,-6.5,length=201)
    dvec <- sapply(x,dd)
    par(las=1,bty="l")
    plot(x,dvec-min(dvec)+1e-16,log="xy",type="b")
    abline(v=getME(FUN(i),"theta"),col=2)
}

## Case #1: boundary estimate with or without boundary.tol
m5 <- tmpf(5)
m5B <- tmpf(5,control=lmerControl(boundary.tol=0))
stopifnot(getME(m5,"theta")==0)
stopifnot(getME(m5B,"theta")==0)
p5 <- profile(m5)  ## bobyqa warnings but results look reasonable
xyplot(p5)
## reveals slight glitch (bottom row of plots doesn't look right)
expect_warning(splom(p5),"unreliable for singular fits")
p5B <- profile(m5,signames=FALSE)
expect_warning(splom(p5B),"unreliable for singular fits")

if(lme4:::testLevel() >= 2) { ## avoid failure to warn
    ## Case #2: near-boundary estimate, but boundary.tol can't fix it
    m16 <- tmpf(16)
    ## tmpplot(16)
    p16 <- profile(m16)  ## warning message (non-monotonic profile)
    if(FALSE) ## FIXME: *does* `warn' but that is not caught by expect_warning() , nor suppressWarning() ???
    if (!.Platform$OS.type=="windows") {
        expect_warning(xyplot(p16),"using linear interpolation")  ## warns about linear interpolation in profile for variable 1
        ## FIXME: don't know why this doesn't warn on windows ...
        ##        ... doesn't warn on mavericks either ...
    }
    ## (still not quite right)
    d16 <- as.data.frame(p16)
    xyplot(.zeta~.focal|.par,data=d16,type=c("p","l"),
           scales=list(x=list(relation="free")))
    try(splom(p16))  ## breaks
}

## bottom line:
##  * xyplot.thpr could still be improved
##  * most of the near-boundary cases are noisy and can't easily be
##    fixed

tmpf2 <- function(i,...) {
    set.seed(i)
    d <- data.frame(x=rnorm(60),f=factor(rep(1:6,each=10)),
                    w=rep(10,60))
    d$y <- simulate(~x+(1|f),family=binomial,
                    weights=d$w,newdata=d,
                    newparams=list(theta=0.01,beta=c(1,1)))[[1]]
    glmer(y~x+(1|f),data=d,family=binomial,weights=w,...)
}

if (FALSE) {
    ## figuring out which seeds will give boundary and
    ## near-boundary solutions
    mList <- lapply(1:201,tmpf2)
    ss <- sapply(mList,sumf)+1e-50
    par(las=1,bty="l")
    hist(log(ss),col="gray",breaks=50)
    ## values lying on boundary
    head(which(log(ss)<(-50)))   ## 1-5, 7 ...
    ## values close to boundary (if check.edge not set)
    which(log(ss)>(-50) & log(ss) <(-20))  ## 44, 46, 52, ...
}

##  m1 <- tmpf2(1)

## FIXME: doesn't work if we generate m1 via tmpf2(1) --
## some environment lookup problem ...

set.seed(1)
d <- data.frame(x=rnorm(60),f=factor(rep(1:6,each=10)),
                w=rep(10,60))
d$y <- simulate(~x+(1|f),family=binomial,
                weights=d$w,newdata=d,
                newparams=list(theta=0.01,beta=c(1,1)))[[1]]
m1 <- glmer(y~x+(1|f),data=d,family=binomial,weights=w)

p1 <- profile(m1)
xyplot(p1)
expect_warning(splom(p1),"splom is unreliable")


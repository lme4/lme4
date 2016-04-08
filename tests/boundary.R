## In both of these cases boundary fit (i.e. estimate of zero RE
## variance) is *incorrect*. (Nelder_Mead, restart_edge=FALSE) is the
## only case where we get stuck; either optimizer=bobyqa or
## restart_edge=TRUE (default) works

library(lme4)
library(testthat)

if(!dev.interactive(orNone=TRUE)) pdf("boundary_plots.pdf")

## Stephane Laurent:
dat <- read.csv(system.file("testdata","dat20101314.csv", package="lme4"))

fit   <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
	      control= lmerControl(optimizer="Nelder_Mead"))
fit_b <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
	      control= lmerControl(optimizer="bobyqa", restart_edge=FALSE))
fit_c <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
	      control= lmerControl(optimizer="Nelder_Mead", restart_edge=FALSE,
				   check.conv.hess="ignore"))
## final fit gives degenerate-Hessian warning
## FIXME: use fit_c with expect_warning() as a check on convergence tests
## tolerance=1e-5 seems OK in interactive use but not in R CMD check ... ??
stopifnot(all.equal(getME(fit,  "theta") -> th.f,
		    getME(fit_b,"theta"), tolerance=5e-5),
	  all(th.f > 0))

## Manuel Koller

source(system.file("testdata", "koller-data.R", package="lme4"))
ldata <- getData(13)
## old (backward compatible/buggy)
fm4  <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="Nelder_Mead",
                                  use.last.params=TRUE),
             start=list(theta=1))

fm4b <- lmer(y ~ (1|Var2), ldata,
             control = lmerControl(optimizer="Nelder_Mead", use.last.params=TRUE,
                                   restart_edge = FALSE,
                                   check.conv.hess="ignore", check.conv.grad="ignore"),
             start = list(theta=1))
## FIXME: use as convergence test check
stopifnot(getME(fm4b,"theta") == 0)
fm4c <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="bobyqa",
                                                      use.last.params=TRUE),
             start=list(theta=1))
stopifnot(all.equal(getME(fm4, "theta") -> th4,
		    getME(fm4c,"theta"), tolerance=1e-4),
	  th4 > 0)


## new: doesn't get stuck at edge any more,  but gets stuck somewhere else ...
fm5 <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="Nelder_Mead",
						     check.conv.hess="ignore",
						     check.conv.grad="ignore"),
	    start=list(theta=1))
fm5b <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="Nelder_Mead",
						      restart_edge=FALSE,
						      check.conv.hess="ignore",
						      check.conv.grad="ignore"),
	     start = list(theta = 1))
fm5c <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="bobyqa"),
	     start = list(theta = 1))
stopifnot(all.equal(unname(getME(fm5c,"theta")), 0.21067645, tolerance = 1e-7))
					#	 0.21067644264 [64-bit, lynne]

##{
    ## additional stuff for diagnosing Nelder-Mead problems.

    library(optimx)
    fm5d <- update(fm5,control=lmerControl(optimizer="optimx",
                       optCtrl=list(method="L-BFGS-B")))

    fm5e <- update(fm5, control=lmerControl(optimizer="nloptwrap"))

    mList <- setNames(list(fm4,fm4b,fm4c,fm5,fm5b,fm5c,fm5d,fm5e),
                      c("NM/uselast","NM/uselast/norestart","bobyqa/uselast",
                        "NM","NM/norestart","bobyqa","LBFGSB","nloptr/bobyqa"))
    pp <- profile(fm5c,which=1)
    dd <- as.data.frame(pp)
    par(las=1,bty="l")
    v <- sapply(mList,
                function(x) sqrt(VarCorr(x)[[1]]))
    plot(.zeta^2~.sig01, data=dd, type="b")
    abline(v=v)

    res <- cbind(VCorr  = sapply(mList, function(x) sqrt(VarCorr(x)[[1]])),
                 theta  = sapply(mList, getME,"theta"),
                 loglik = sapply(mList, logLik))
    res
    print(sessionInfo(), locale=FALSE)
##}

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
    mList <- lapply(1:201,tmpf) # [FIXME tons of messages "theta parameters vector not named"]
    ss <- sapply(mList,sumf)+1e-50
    par(las=1,bty="l")
    hist(log(ss),col="gray",breaks=50)
    ## values lying on boundary
    which(log(ss)<(-40))   ## 5, 7-13, 15, 21, ...
    ## values close to boundary (if check.edge not set)
    which(log(ss)>(-40) & log(ss) <(-20))  ## 16, 44, 80, 86, 116, ...
}
## diagnostic plot
tmpplot <- function(i, FUN=tmpf) {
    dd <- FUN(i, devFunOnly=TRUE)
    x <- 10^seq(-10,-6.5,length=201)
    dvec <- sapply(x,dd)
    op <- par(las=1,bty="l"); on.exit(par(op))
    plot(x,dvec-min(dvec)+1e-16, log="xy", type="b")
    r <- FUN(i)
    abline(v = getME(r,"theta"), col=2)
    invisible(r)
}

## Case #1: boundary estimate with or without boundary.tol
m5  <- tmpf(5)
m5B <- tmpf(5,control=lmerControl(boundary.tol=0))
stopifnot(getME(m5, "theta")==0,
          getME(m5B,"theta")==0)
p5 <- profile(m5)  ## bobyqa warnings but results look reasonable
xyplot(p5)
## reveals slight glitch (bottom row of plots doesn't look right)
expect_warning(splom(p5),"unreliable for singular fits")
p5B <- profile(m5, signames=FALSE) # -> bobyqa convergence warning (code 3)
expect_warning(splom(p5B), "unreliable for singular fits")

if(lme4:::testLevel() >= 2) { ## avoid failure to warn
    ## Case #2: near-boundary estimate, but boundary.tol can't fix it
    m16 <- tmpplot(16)
    ## sometimes[2014-11-11] fails (??) :
    p16 <- profile(m16)  ## warning message*s* (non-monotonic profile and more)
    plotOb <- xyplot(p16)
    ## NB: It's the print()ing of 'plotOb' which warns ==> need to do this explicitly:
    expect_warning(print(plotOb), ## warns about linear interpolation in profile for variable 1
                   "using linear interpolation")
    d16 <- as.data.frame(p16)
    xyplot(.zeta ~ .focal|.par, data=d16, type=c("p","l"),
           scales = list(x=list(relation="free")))
    try(splom(p16))  ## breaks when calling predict(.)
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


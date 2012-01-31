### suppressPackageStartupMessages(...)  as we have an *.Rout.save to Rdiff against
stopifnot(suppressPackageStartupMessages(require(lme4Eigen)))
options(show.signif.stars = FALSE)

source(system.file("test-tools.R", package = "Matrix"))# identical3() etc
all.EQ <- function(u,v, ...) all.equal.X(u, v, except = c("call", "frame"), ...)
S4_2list <- function(obj) {   # no longer used
   sn <- slotNames(obj)
   structure(lapply(sn, slot, object = obj), .Names = sn)
}
## Is now (2010-09-03) in Matrix' test-tools.R above
## showProc.time <- local({
##     pct <- proc.time()
##     function() { ## CPU elapsed __since last called__
## 	ot <- pct ; pct <<- proc.time()
## 	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
##     }
## })

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm1a <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE))
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
anova(fm1, fm2)

## Now works for glmer
fm1. <- glmer(Reaction ~ Days + (Days|Subject), sleepstudy)
## default family=gaussian -> automatically calls  lmer()
stopifnot(all.equal(fm1, fm1.))
## Test against previous version in lmer1 (using bobyqa for consistency)
#(fm1. <- lmer1(Reaction ~ Days + (Days|Subject), sleepstudy, opt = "bobyqa"))
#stopifnot(all.equal(fm1@devcomp$cmp['REML'], fm1.@devcomp$cmp['REML']),
#          all.equal(fixef(fm1), fixef(fm1.)),
#          all.equal(fm1@re@theta, fm1.@theta, tol = 1.e-7),
#          all.equal(ranef(fm1), ranef(fm1.)))

## compDev = FALSE no longer applies to lmer
## Test 'compDev = FALSE' (vs TRUE)
## fm1. <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
##              compDev = FALSE)#--> use R code (not C++) for deviance computation
## stopifnot(all.equal(fm1@devcomp$cmp['REML'], fm1.@devcomp$cmp['REML']),
##           all.equal(fixef(fm1), fixef(fm1.)),
##           all.equal(fm1@re@theta, fm1.@re@theta, tol = 1.e-7),
##           all.equal(ranef(fm1), ranef(fm1.), tol = 1.e-7))


stopifnot(all.equal(fixef(fm1), fixef(fm2), tol = 1.e-13),
          all.equal(unname(fixef(fm1)),
                    c(251.405104848485, 10.467285959595), tol = 1e-13),
	  all.equal(Matrix::cov2cor(vcov(fm1))["(Intercept)", "Days"],
		    -0.13755, tol=1e-4))

fm1ML <- refitML(fm1)
fm2ML <- refitML(fm2)
print(AIC(fm1ML)); print(AIC(fm2ML))
print(BIC(fm1ML)); print(BIC(fm2ML))

(fm3 <- lmer(Yield ~ 1|Batch, Dyestuff2))
stopifnot(all.equal(coef(summary(fm3)),
		    array(c(5.6656, 0.67838803150, 8.3515624346),
			  c(1,3), dimnames = list("(Intercept)",
				  c("Estimate", "Std. Error", "t value")))))
showProc.time() #

### {from ../man/lmer.Rd } --- compare lmer & lmer1 ---------------
(fmX1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm.1 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))

#(fmX2 <- lmer2(Reaction ~ Days + (Days|Subject), sleepstudy))
#(fm.2 <- lmer2(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
## check update(<mer>, <formula>):
fm.3 <- update(fmX1, . ~ Days + (1|Subject) + (0+Days|Subject))
stopifnot(all.equal(fm.1, fm.3))

fmX1s <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, sparseX=TRUE)
#fmX2s <- lmer2(Reaction ~ Days + (Days|Subject), sleepstudy, sparseX=TRUE)

showProc.time() #

for(nm in c("coef", "fixef", "ranef", "sigma",
	     "model.matrix", "model.frame" , "terms")) {
    cat(sprintf("%15s : ", nm))
    FUN <- get(nm)
    F.fmX1s <- FUN(fmX1s)
#    F.fmX2s <- FUN(fmX2s)
#    if(nm == "model.matrix") {
#        F.fmX1s <- as(F.fmX1s, "denseMatrix")
#        F.fmX2s <- as(F.fmX2s, "denseMatrix")
#	FF <- function(.) {r <- FUN(.); row.names(r) <- NULL
#			   as(r, "generalMatrix") }
#    } # else
    FF <- FUN
    stopifnot(
	      all.equal( FF(fmX1), F.fmX1s, tol =  1e-6)
#	      ,
#	      all.equal( FF(fmX2), F.fmX2s, tol = 1e-5)
#              ,
#	      all.equal( FF(fm.1), F.fmX2s, tol = 9e-6) ## these are different models
#              ,
#              all.equal(F.fmX2s,   F.fmX1s, tol = 6e-6)
#              ,
#              all.equal(FUN(fm.1), FUN(fm.2), tol = 6e-6)
              ,
              TRUE)
    cat("[Ok]\n")
}


## transformed vars should work[even if non-sensical as here;failed in 0.995-1]
fm2l <- lmer(log(Reaction) ~ log(Days+1) + (log(Days+1)|Subject),
             data = sleepstudy, REML = FALSE)
## no need for an expand method now : xfm2 <- expand(fm2)

stopifnot(dim(ranef(fm2l)[[1]]) == c(18, 2),
          is((c3 <- coef(fm3)), "coef.mer"),
          all(fixef(fm3) == c3$Batch),## <-- IFF  \hat{\sigma^2} == 0
          TRUE)

## generalized linear mixed model
## TODO: (1) move these to ./glmer-ex.R
## ----  (2) "rationalize" with ../man/cbpp.Rd
#m1e <- glmer1(cbind(incidence, size - incidence) ~ period + (1 | herd),
#              family = binomial, data = cbpp, doFit = FALSE)
## now
#bobyqa(m1e, control = list(iprint = 2L))
m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), cbpp, binomial, nAGQ=25L)
dput(unname(fixef(m1)))
dput(deviance(m1))
stopifnot(is((cm1 <- coef(m1)), "coef.mer"),
	  dim(cm1$herd) == c(15,4),
	  all.equal(fixef(m1), ##  these values are from an Ubuntu 11.10 amd64 system
                    c(-1.39922533406847, -0.991407294757321,
                      -1.12782184600404, -1.57946627431248),
		    tol = 1.e-5,
                    check.attr=FALSE),
          all.equal(deviance(m1), 100.010030538022, tol=1e-9)
	  )

## Simple example by Andrew Gelman (2006-01-10) ----
n.groups <- 10 ; n.reps <- 2
n <- length(group.id <- gl(n.groups, n.reps))
## simulate the varying parameters and the data:
set.seed(0)
a.group <- rnorm(n.groups, 1, 2)
y <- rnorm (n, a.group[group.id], 1)
## fit and summarize the model
fit.1 <- lmer (y ~ 1 + (1 | group.id))
coef (fit.1)
## check show( <"summary.mer"> ):
(sf1 <- summary(fit.1)) # --> now looks as for fit.1

stopifnot(all.equal(fixef(fit.1), c("(Intercept)" = 1.571312129)),
	  all.equal(unname(ranef(fit.1, drop=TRUE)[["group.id"]]),
		   c(1.8046888, -1.8097665, 1.6146451, 1.5408268, -0.1331995,
                     -3.3306655, -1.8259277, -0.8735145, -0.3591311,  3.3720441),
		    tol = 1e-5)
	  )


## ranef and coef
rr <- ranef(fm1)
stopifnot(is.list(rr), length(rr) == 1, class(rr[[1]]) == "data.frame")
print(plot(rr))
stopifnot(is(cc <- coef(fm1), "coef.mer"),
	  is.list(cc), length(cc) == 1, class(cc[[1]]) == "data.frame")
print(plot(cc))
rr <- ranef(fm2)
stopifnot(is.list(rr), length(rr) == 1, class(rr[[1]]) == "data.frame")
print(plot(rr))
stopifnot(is(cc <- coef(fm2), "coef.mer"),
	  is.list(cc), length(cc) == 1, class(cc[[1]]) == "data.frame")
print(plot(cc))

showProc.time() #

if (require('MASS', quietly = TRUE)) {
    bacteria$wk2 <- bacteria$week > 2
    contrasts(bacteria$trt) <-
        structure(contr.sdif(3),
                  dimnames = list(NULL, c("diag", "encourage")))
    print(fm5 <- glmer(y ~ trt + wk2 + (1|ID), bacteria, binomial, nAGQ=25L))
    ## used to fail with nlminb() : stuck at theta=1

    showProc.time() #

    stopifnot(
	      all.equal(logLik(fm5),
			## was	  -96.127838
			structure(-95.89706, nobs = 220L, nall = 220L,
				  df = 5L, REML = FALSE,
                                  class = "logLik"),
                        tol = 1e-5, check.attributes = FALSE)
	      ,
	      all.equal(fixef(fm5),
                        c("(Intercept)"= 2.85970407987798, "trtdiag"= -1.36896064622876,
                          "trtencourage"=0.579864265133904, "wk2TRUE"=-1.62687300090319),
                        tol = 1e-6)
	      )
}

## Invalid factor specification -- used to seg.fault:
set.seed(1)
dat <- within(data.frame(lagoon = factor(rep(1:4,each = 25)),
                         habitat = factor(rep(1:20, each = 5))),
          {
              y <- round(10*rnorm(100, m = 10*as.numeric(lagoon)))
          })

try(reg <- lmer(y ~ habitat + (1|habitat*lagoon), data = dat) # did seg.fault
    ) # now gives error                 ^- should be ":"
r1  <- lmer(y ~ 0+habitat + (1|habitat:lagoon), data = dat) # ok, but senseless
r1b <- lmer(y ~ 0+habitat + (1|habitat), data = dat) # same model, clearly indeterminable
## "TODO" :  summary(r1)  should ideally warn the user
stopifnot(all.equal(fixef(r1), fixef(r1b), tol= 1e-15),
          all.equal(ranef(r1), ranef(r1b), tol= 1e-15, check.attributes=FALSE))

## Use a more sensible model:
r2.0 <- lmer(y ~ 0+lagoon + (1|habitat:lagoon), data = dat) # ok
r2   <- lmer(y ~ 0+lagoon + (1|habitat), data = dat) # ok, and more clear
stopifnot(all.equal(fixef(r2), fixef(r2.0), tol= 1e-15),
          all.equal(ranef(r2), ranef(r2.0), tol= 1e-15, check.attributes=FALSE))
V2 <- vcov(r2)
assert.EQ.mat(V2, diag(x = 9.9833/3, nr = 4))
stopifnot(all.equal(unname(fixef(r2)) - (1:4)*100,
		    c(1.72, 0.28, 1.76, 0.8), tol = 1e-13))

## sparseX version should give same numbers:
r2.  <- lmer(y ~ 0+lagoon + (1|habitat), data = dat,
             sparseX = TRUE, verbose = TRUE)

## the summary() components we do want to compare 'dense X' vs 'sparse X':
nmsSumm <- c("methTitle", "devcomp", "logLik", "ngrps", "coefficients",
             "sigma", "REmat", "AICtab")
sr2  <- summary(r2)
sr2. <- summary(r2.)
sr2.$devcomp$dims['spFe'] <- 0L       # to allow for comparisons below
stopifnot(all.equal(sr2[nmsSumm], sr2.[nmsSumm], tol= 1e-14)
          , all.equal(ranef(r2), ranef(r2.), tol= 1e-14)
          , Matrix:::isDiagonal(vcov(r2.)) # ok
          , all.equal(Matrix::diag(vcov(r2.)), rep.int(V2[1,1], 4), tol= 1e-13)
#          , all(vcov(r2.)@factors$correlation == diag(4))  # not sure why this fails
          , TRUE)
r2.

## Failure to specify a random effects term - used to give an obscure message
## Ensure *NON*-translated message; works on Linux,... :
if(.Platform$OS.type == "unix") {
Sys.setlocale("LC_MESSAGES", "C")
tc <- tryCatch(
	       m2 <- glmer(incidence / size ~ period, weights = size,
			   family = binomial, data = cbpp)
	       , error = function(.) .)
stopifnot(inherits(tc, "error"),
	  identical(tc$message,
		    "No random effects terms specified in formula"))
}

### mcmcsamp() :
## From: Andrew Gelman <gelman@stat.columbia.edu>
## Date: Wed, 18 Jan 2006 22:00:53 -0500

if (FALSE) {  # mcmcsamp still needs work
    ## NB: Need to restore coda to the Suggests: field of DESCRIPTION
    ## file if this code block is reinstated.
    ## has.coda <- require(coda)
    ## if(!has.coda)
    ##     cat("'coda' package not available; some outputs will look suboptimal\n")

    ## Very simple example
    y <- 1:10
    group <- gl(2,5)
    (M1 <- lmer (y ~ 1 + (1 | group))) # works fine
    (r1 <- mcmcsamp (M1))              # dito
    r2 <- mcmcsamp (M1, saveb = TRUE)  # gave error in 0.99-* and 0.995-[12]
    (r10 <- mcmcsamp (M1, n = 10, saveb = TRUE))

    ## another one, still simple
    y <- (1:20)*pi
    x <- (1:20)^2
    group <- gl(2,10)
    M1 <- lmer (y ~ 1 | group)
    mcmcsamp (M1, n = 2, saveb=TRUE) # fine

    M2 <- lmer (y ~ 1 + x + (1 + x | group)) # false convergence
    ## should be identical (and is)
    M2 <- lmer (y ~ x + ( x | group))#  false convergence -> simulation doesn't work:
    if(FALSE) ## try(..) fails here (in R CMD check) [[why ??]]
        mcmcsamp (M2, saveb=TRUE)
    ## Error: inconsistent degrees of freedom and dimension ...

    ## mcmc for glmer:
    rG1k <- mcmcsamp(m1, n = 1000)
    summary(rG1k)
    rG2 <- mcmcsamp(m1, n = 3, verbose = TRUE)
}

## Spencer Graves' example (from a post to S-news, 2006-08-03) ----------------
## it should give an error, rather than silent non-sense:
tstDF <- data.frame(group = letters[1:5], y = 1:5)
assertError(## Now throws an error, as desired :
            lmer(y ~ 1 + (1|group), data = tstDF)
            )

showProc.time() #

## Wrong formula gave a seg.fault at times:
set.seed(2)# !
D <-  data.frame(y= rnorm(12,10), ff = gl(3,2,12),
                 x1=round(rnorm(12,3),1), x2=round(rnorm(12,7),1))
## NB: The first two are the same, having a length-3 R.E. with 3 x 3 vcov-matrix:
## --> do need CPU
m0 <- lmer(y ~ (x1 + x2)|ff, data = D)
m1 <- lmer(y ~ x1 + x2|ff  , data = D)
m2 <- lmer(y ~ x1 + (x2|ff), data = D)
m3 <- lmer(y ~ (x2|ff) + x1, data = D)
stopifnot(all.equal(ranef(m0), ranef(m1)),
          all.equal(ranef(m2), ranef(m3)),
          inherits(tryCatch(lmer(y ~ x2|ff + x1, data = D), error = function(e)e),
                   "error"))

showProc.time() #

## Reordering of grouping factors should not change the internal structure
#Pm1  <- lmer1(strength ~ (1|batch) + (1|sample), Pastes, doFit = FALSE)
#Pm2  <- lmer1(strength ~ (1|sample) + (1|batch), Pastes, doFit = FALSE)
#P2.1 <- lmer (strength ~ (1|batch) + (1|sample), Pastes, devFunOnly = TRUE)
#P2.2 <- lmer (strength ~ (1|sample) + (1|batch), Pastes, devFunOnly = TRUE)

## The environments of Pm1 and Pm2 should be identical except for
## "call" and "frame":
#stopifnot(## all.EQ(env(Pm1), env(Pm2)),
#	  all.EQ(S4_2list(P2.1),
#		 S4_2list(P2.2)))

## glmer - Modeling overdispersion as "mixture" aka
## ----- - *ONE* random effect *PER OBSERVATION" -- example inspired by Ben Bolker:

##' <description>
##'
##' <details>
##' @title
##' @param ng number of groups
##' @param nr number of "runs", i.e., observations per groups
##' @param sd standard deviations of group and "Individual" random effects,
##'    (\sigma_f, \sigma_I)
##' @param b  true beta (fixed effects)
##' @return a data frame (to be used in glmer()) with columns
##'    (x, f, obs, eta0, eta, mu, y), where y ~ Pois(lambda(x)),
##'                                   log(lambda(x_i)) = b_1 + b_2 * x + G_{f(i)} + I_i
##'    and G_k ~ N(0, \sigma_f);  I_i ~ N(0, \sigma_I)
##' @author Ben Bolker and Martin Maechler
set.seed(1)
rPoisGLMMi <- function(ng, nr, sd=c(f = 1, ind = 0.5), b=c(1,2))
{
  stopifnot(nr >= 1, ng >= 1,
            is.numeric(sd), names(sd) %in% c("f","ind"), sd >= 0)
  ntot <- nr*ng
  b.reff <- rnorm(ng,  sd= sd[["f"]])
  b.rind <- rnorm(ntot,sd= sd[["ind"]])
  x <- runif(ntot)
  within(data.frame(x,
                    f = factor(rep(LETTERS[1:ng], each=nr)),
                    obs = 1:ntot,
                    eta0 = cbind(1, x) %*% b),
     {
         eta <- eta0 + b.reff[f] + b.rind[obs]
         mu <- exp(eta)
         y <- rpois(ntot, lambda=mu)
     })
}
dd <- rPoisGLMMi(12, 20)
m0  <- glmer(y~x + (1|f),           family="poisson", data=dd)
(m1 <- glmer(y~x + (1|f) + (1|obs), family="poisson", data=dd))# must use Laplace
anova(m0, m1)

showProc.time()

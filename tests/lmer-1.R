### suppressPackageStartupMessages(...)  as we have an *.Rout.save to Rdiff against
stopifnot(suppressPackageStartupMessages(require(lme4)))
options(show.signif.stars = FALSE)

(fm1 <-  lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm1a <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE))
(fm2 <-  lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))

## transformed vars should work[even if non-sensical as here;failed in 0.995-1]
fm2l <- lmer(log(Reaction) ~ log(Days+1) + (log(Days+1)|Subject),
             data = sleepstudy, REML = FALSE)
xfm2 <- expand(fm2l)
stopifnot(is(fm1, "mer"), is(fm2l, "mer"),
          is(xfm2$P, "sparseMatrix"))

AIC(fm1); AIC(fm2)
BIC(fm1); BIC(fm2)
## not yet: if(getRversion() > "2.11.0") {
##  AIC(fm1, fm2)
##  BIC(fm1, fm2)
## }

## generalized linear mixed model
(m1 <- lmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = binomial, data = cbpp))
stopifnot(is(m1,"mer"), is((cm1 <- coef(m1)), "coef.mer"),
	  dim(cm1$herd) == c(15,4),
	  all.equal(fixef(m1),
		    c(-1.39853504914, -0.992334711,
		      -1.12867541477, -1.58037390498), check.attr=FALSE)
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
coef (fit.1)# failed in Matrix 0.99-6
(sf1 <- summary(fit.1)) # show() is as without summary()


## ranef and coef
rr <- ranef(fm1)
stopifnot(is.list(rr), length(rr) == 1, class(rr[[1]]) == "data.frame")
print(plot(rr))
cc <- coef(fm1)
stopifnot(is.list(cc), length(cc) == 1, class(cc[[1]]) == "data.frame")
print(plot(cc))
rr <- ranef(fm2)
stopifnot(is.list(rr), length(rr) == 1, class(rr[[1]]) == "data.frame")
print(plot(rr))
cc <- coef(fm2)
stopifnot(is.list(cc), length(cc) == 1, class(cc[[1]]) == "data.frame")
print(plot(cc))

if (require('MASS', quietly = TRUE)) {
    bacteria$wk2 <- bacteria$week > 2
    contrasts(bacteria$trt) <-
        structure(contr.sdif(3),
                  dimnames = list(NULL, c("diag", "encourage")))
    print(fm5 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial))
    print(fm6 <- lmer(y ~ trt + wk2 + (1|ID), bacteria, binomial))
}

## Invalid factor specification -- used to seg.fault:
set.seed(1)
dat <- data.frame(y = round(10*rnorm(100)), lagoon = factor(rep(1:4,each = 25)),
                  habitat = factor(rep(1:20, each = 5)))
r1  <- lmer(y ~ habitat + (1|habitat:lagoon), data = dat) # ok

try(
    reg <- lmer(y ~ habitat + (1|habitat*lagoon), data = dat) # did seg.fault
    ) # now gives error                 ^- should be ":"


## Failure to specify a random effects term - used to give an obscure message
try(
m2 <- lmer(incidence / size ~ period, weights = size,
            family = binomial, data = cbpp)
)

## needed as MacOSX (i386) gave differing results in some cases,
## (even on *repeating* the *same* computations!) :
if(Sys.info()[["sysname"]] == "Darwin") identical <- all.equal ## the horror..!

### mcmcsamp() :
## From: Andrew Gelman <gelman@stat.columbia.edu>
## Date: Wed, 18 Jan 2006 22:00:53 -0500

    has.coda <- require(coda)
    if(!has.coda)
        cat("'coda' package not available; some outputs will look suboptimal\n")

    ## Very simple example
    y <- 1:10
    group <- gl(2,5)
    (M1 <- lmer (y ~ 1 + (1 | group))) # works fine
    set.seed(25)
    (r1 <- mcmcsamp (M1))              # dito
    r2 <- mcmcsamp (M1, saveb = TRUE)  # gave error in 0.99-* and 0.995-[12]
    (r10 <- mcmcsamp (M1, n = 10, saveb = TRUE))

    ## another one, still simple
    y <- (1:20)*pi
    x <- (1:20)^2
    group <- gl(2,10)
    M1 <- lmer (y ~ 1 | group)
    mcmcsamp (M1, n = 2, saveb=TRUE) # fine

    M2. <- lmer (y ~ 1 + x + (1 + x | group)) # had false convergence
    ## convergence now ok (but ranef corr. = -1; fixef = -.996 :
    summary(M2.)
    M2  <- lmer (y ~ x + ( x | group))
    ## should be identical (and is .. well, not on all versions on Mac OSX):
    stopifnot(identical(fixef(M2), fixef(M2.)),
	      identical(ranef(M2), ranef(M2.)),
	      identical(resid(M2), resid(M2.)))

    mcmcsamp (M2, saveb=TRUE)

if (FALSE) {  # mcmcsamp for  glmer  not yet available

    ## mcmc for glmer:
    rG1k <- mcmcsamp(m1, n = 1000)
    summary(rG1k)
    rG2 <- mcmcsamp(m1, n = 3, verbose = TRUE)
}

## Spencer Graves' example (from a post to S-news, 2006-08-03): ----------------
## FIXME?
tstDF <- data.frame(group = letters[1:5], y = 1:5)
var(tstDF$y) # == 2.5
## Now throws an error
try(f.oops <- lmer(y ~ 1 + (1|group), data = tstDF))
##  summary(f.oops) ## or print(Matrix:::formatVC(VarCorr(f.oops)), quote = FALSE)
## ...
##   Groups   Name        Variance Std.Dev.
##   group    (Intercept) 1.81818  1.34840
##   Residual             0.68182  0.82572
## ...
##
##SG>	 This is ... silly, because there are zero degrees of freedom
##SG> to distinguish "group" from Residual.  It is comforting that the sum of
##SG> the variances sum to the variance of "y", ......
##SG>	 However, I would prefer to have the multilevel software catch this
##SG> case and optionally return an error or drop the redundant group
##SG> with a warning.

## Wrong formula gave a seg.fault at times:
D <-  data.frame(y= rnorm(20,10), ff = gl(4,5),
                 x1=rnorm(20,3), x2=rnorm(20,7),
                 x3=rnorm(20,1))
m0 <- lmer(y ~ (x1 + x2)|ff, data = D)
m1 <- lmer(y ~ x1 + x2|ff  , data = D)
m2 <- lmer(y ~ x1 + (x2|ff), data = D)
m3 <- lmer(y ~ (x2|ff) + x1, data = D)
stopifnot(identical(ranef(m0), ranef(m1)),
	  identical(ranef(m2), ranef(m3)))
stopifnot(inherits(tryCatch(lmer(y ~ x2|ff + x1, data = D), error = function(e)e),
		   "error"))

## Check the use of offset
om2 <- lmer(y ~ x1 + (x2|ff), data = D, offset = x3)
om3 <- lmer(y ~ x1 + (x2|ff) + offset(x3), data = D)

stopifnot(identical(ranef(om2), ranef(om3)),
          identical(deviance(om2), deviance(om3)))
if (identical(TRUE, all.equal(fixef(m2), fixef(om2))))
    stop("offset does not change the fixed effects")

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
(m1 <- glmer(y~x + (1|f) + (1|obs), family="poisson", data=dd))
anova(m0, m1)

## in case Mac useRs source() this file:
if(Sys.info()[["sysname"]] == "Darwin") rm(identical)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

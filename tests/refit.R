#### Testing  refit()
#### ----------------

library(lme4)
set.seed(101)

## for each type of model, should be able to
##  (1) refit with same data and get the same answer,
##     at least structurally (small numerical differences
##     are probably unavoidable)
##  (2) refit with simulate()d data

getinfo <- function(x) {
  c(fixef(x), logLik(x), unlist(ranef(x)), unlist(VarCorr(x)))
}

dropterms <- function(x) {
    attr(x@frame,"terms") <- NULL
    x
}

if (getRversion() >= "3.0.0") {
    attach(system.file("testdata", "lme-tst-fits.rda", package="lme4"))
} else {
    ## saved fits are not safe with old R versions; just re-compute ("cheat"!):

    fit_sleepstudy_2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

    cbpp$obs <- factor(seq(nrow(cbpp)))
    ## intercept-only fixed effect
    fit_cbpp_0 <- glmer(cbind(incidence, size-incidence) ~ 1 + (1|herd),
                        cbpp, family=binomial)
    ## include fixed effect of period
    fit_cbpp_1 <- update(fit_cbpp_0, . ~ . + period)
    if(FALSE) ## include observation-level RE
        fit_cbpp_2 <- update(fit_cbpp_1, . ~ . + (1|obs))
    ## specify formula by proportion/weights instead
    fit_cbpp_3 <- update(fit_cbpp_1, incidence/size ~ period + (1 | herd), weights = size)
}


## LMM
fm1 <- fit_sleepstudy_2
fm1R <- refit(fm1, sleepstudy$Reaction)
fm1S <- refit(fm1, simulate(fm1)[[1]])

stopifnot(all.equal(getinfo(fm1 ),
                    getinfo(fm1R), tolerance = 6e-3),
          all.equal(getinfo(fm1 ),
                    getinfo(fm1S), tolerance = 0.5) # <- simulate()
          )

if(FALSE) { ## show all differences
    sapply(slotNames(fm1), function(.)
        all.equal( slot(fm1,.), slot(fm1R,.), tolerance=0))
}

if (getRversion() >= "3.4.0") {
    ## differences: FALSE for resp, theta, u, devcomp, pp, optinfo?
    ## FIXME: this isn't actually tested in any way ...
    sapply(slotNames(fm1),
           function(.) isTRUE(all.equal( slot(fm1,.), slot(fm1R,.), tolerance= 1.5e-5)))
    str(fm1 @ optinfo)
    str(fm1R@ optinfo)
}

fm1ML <- refitML(fm1)
stopifnot(
    all.equal(getinfo(fm1), getinfo(fm1ML), tolerance=0.05)# 0.029998
)


## binomial GLMM (two-column)
gm1 <- fit_cbpp_1
gm1R <- refit(gm1, with(cbpp, cbind(incidence,size-incidence)))

 sim1Z <- simulate(gm1)[[1]]
 sim1Z[4,] <- c(0,0)
 (gm1. <- refit(gm1, sim1Z)) # earlier gave Error:  ... PIRLS ... failed ...

all.equal(getinfo(gm1), getinfo(gm1R), tolerance=0) # to see it --> 5.52e-4
# because glmer() uses Laplace approx. (? -- still, have *same* y !)
stopifnot(all.equal(getinfo(gm1), getinfo(gm1R), tolerance = 1e-4))

gm1S <- refit(gm1, simulate(gm1)[[1]])
          all.equal(getinfo(gm1), getinfo(gm1S), tolerance=0) # to see:
stopifnot(all.equal(getinfo(gm1), getinfo(gm1S), tolerance = 0.4))

## binomial GLMM (prob/weights)
formula(gm2 <- fit_cbpp_3)
## glmer(incidence/size ~ period + (1 | herd), cbpp, binomial, weights=size)
gm2R <- refit(gm2,  with(cbpp, incidence/size))
          all.equal(getinfo(gm2), getinfo(gm2R), tolerance= 0)
stopifnot(all.equal(getinfo(gm2), getinfo(gm2R), tolerance= 6e-4))

## FIXME: check on Windows == 2015-06: be brave
gm2S <- refit(gm2, simulate(gm2)[[1]])
          all.equal(getinfo(gm2), getinfo(gm2S), tolerance=0)# 0.17 .. upto 0.28
stopifnot(all.equal(getinfo(gm2), getinfo(gm2S), tolerance=0.40))

## from Alexandra Kuznetsova
set.seed(101)
Y <- matrix(rnorm(1000),ncol=2)
d <- data.frame(y1=Y[,1],  x=rnorm(100), f=rep(1:10,10))
fit1 <- lmer(y1 ~ x+(1|f),data=d)
fit2 <- refit(fit1, newresp = Y[,2], rename.response=TRUE)
## check, but ignore terms attribute of model frame ...
tools::assertWarning(refit(fit1, newresp = Y[,2], junk=TRUE))
if (isTRUE(all.equal(fit1,fit2))) stop("fit1 and fit2 should not be equal")
## hack number of function evaluations
u2 <- update(fit2)
fit2@optinfo$feval <- u2@optinfo$feval <-  NA

d1 <- dropterms(fit2)
d2 <- dropterms( u2 )
## They are not "all equal", but mostly :
for (i in slotNames(d1)) {
    ae <- all.equal(slot(d1,i), slot(d2,i))
    cat(sprintf("%10s: %s\n", i, if(isTRUE(ae)) "all.equal"
		else paste(ae, collapse="\n  ")))
}
          all.equal(getinfo(d1), getinfo(d2),  tolerance = 0)# -> 0.00126
stopifnot(all.equal(getinfo(d1), getinfo(d2),  tolerance = 0.005))


## Bernoulli GLMM (specified as factor)
if (requireNamespace("mlmRev")) {
    data(Contraception, package="mlmRev")
    gm3 <- glmer(use ~ urban + age + livch + (1|district),
                 Contraception, binomial)
    gm3R <- refit(gm3, Contraception$use)
    gm3S <- refit(gm3, simulate(gm3)[[1]])
    stopifnot(all.equal(getinfo(gm3 ),
                        getinfo(gm3R), tolerance = 1e-5),# 64b_Lx: 7.99e-7
              all.equal(getinfo(gm3 ),
                        getinfo(gm3S), tolerance = 0.05) # <- simulated data
              )
    cat("gm3: glmer(..):\n"        ); print(getinfo(gm3))
    cat("gm3R: refit(*, y):\n"     ); print(getinfo(gm3R))
    cat("gm3S: refit(*, sim.()):\n"); print(getinfo(gm3S))

    data(Mmmec, package="mlmRev")
    if (lme4:::testLevel() > 1) {
        gm4 <- glmer(deaths ~ uvb + (1|region), data=Mmmec,
                     family = poisson,
                     offset = log(expected))
        ## FIXME: Fails to converge (with larger maxit: "downdate .. not pos.def..")
        try( gm4R <- refit(gm4, Mmmec $ deaths) )
        try( gm4S <- refit(gm4, simulate(gm4)[[1]]) )
        if(FALSE) { ## FIXME (above)
            cat("gm4R: refit(*,y):\n" ); print( getinfo(gm4R) )
            cat("gm4S: refit(*,y):\n" ); print( getinfo(gm4S) )
            stopifnot(all.equal(getinfo(gm4),getinfo(gm4R),tolerance=6e-5))
        }
    }
}

## ----------------------------------------------------------------------
## issue: #231, http://ms.mcmaster.ca/~bolker/misc/boot_reset.html
## commits: 1a34cd0, e33d698, 53ce966, 7dbfff1, 73aa1bb, a693ba9, 8dc8cf0
## ----------------------------------------------------------------------

formGrouse <- TICKS ~ YEAR + scale(HEIGHT) + (1 | BROOD) + (1 | INDEX) + (1 | LOCATION)
gmGrouse <- glmer(formGrouse, family = "poisson", data = grouseticks)
set.seed(105)
simTICKS <- simulate(gmGrouse)[[1]]
newdata <- transform(grouseticks, TICKS = simTICKS)
gmGrouseUpdate <-  update(gmGrouse, data = newdata)
gmGrouseRefit  <-  refit(gmGrouse, newresp = simTICKS)

## compute and print tolerances
all.equal(bet.U <- fixef(gmGrouseUpdate),
          bet.R <- fixef(gmGrouseRefit), tolerance = 0)
all.equal(th.U <- getME(gmGrouseUpdate, "theta"),
          th.R <- getME(gmGrouseRefit,  "theta"), tolerance = 0)
all.equal(dev.U <- deviance(gmGrouseUpdate),
          dev.R <- deviance(gmGrouseRefit),	tolerance = 0)
stopifnot(
    all.equal(bet.U, bet.R, tolerance = 6e-5), # saw 1.0e-5
    all.equal( th.U,  th.R, tolerance = 4e-5), # saw 1.2e-5
    all.equal(dev.U, dev.R, tolerance = 2e-5)) # saw 4.6e-6

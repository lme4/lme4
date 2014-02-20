stopifnot(require(lme4))

(testLevel <- lme4:::testLevel())

## "MEMSS" is just 'Suggest' -- must still work, when it's missing:
if (suppressWarnings(!require(MEMSS,quietly=TRUE)) ||
    (data(ergoStool, package="MEMSS") != "ergoStool")) {

    cat("'ergoStool' data from package 'MEMSS' is not available --> skipping test\n")
} else {

fm1   <-  lmer (effort ~ Type + (1|Subject), data = ergoStool)
fm1.s  <- lmer (effort ~ Type + (1|Subject), data = ergoStool, sparseX=TRUE)
## was segfaulting with sparseX (a while upto 2010-04-06)

fe1   <- fixef(fm1)
fe1.s <- fixef(fm1.s)

print(s1.d <- summary(fm1))
print(s1.s <- summary(fm1.s))
stopifnot(
	  all.equal(fe1, fe1.s, tolerance= 1e-12)
	  ,
	  all.equal(se1.d <- coef(s1.d)[,"Std. Error"],
		    se1.s <- coef(s1.s)[,"Std. Error"])#, tol = 1e-10
	  ,
	  all.equal(V.d <- vcov(fm1),
		    V.s <- vcov(fm1.s))#, tol = 1e-9
	  ,
	  all.equal(Matrix::diag(V.d), unname(se1.d)^2, tolerance= 1e-12)
	  ,
	  all.equal(unname(se1.d),
		    if(fm1@optinfo$optimizer == "Nelder_Mead")
		    c(0.57601196, rep(0.51868399, 3)) else ## "bobyqa"
		    c(0.57601226, rep(0.51868384, 3)),
		    tolerance = 1e-6) ## std.err.: no extreme accuracy needed
          )

}## if( ergoStool is available from pkg MEMSS )

### -------------------------- a "large" example -------------------------
str(InstEval)

if (FALSE) { # sparse X is not currently implemented, so forget about this:

system.time(## works with 'sparseX'; d has 1128 levels
fm7 <- lmer(y ~ d + service + studage + lectage + (1|s),
             data = InstEval, sparseX=TRUE, verbose=1L, REML=FALSE)
)
system.time(sfm7 <- summary(fm7))
fm7 # takes a while as it computes summary() again !

range(t.fm7 <- coef(sfm7)[,"t value"])## -10.94173  10.61535 for REML, -11.03438  10.70103 for ML

m.t.7 <- mean(abs(t.fm7), trim = .01)
#stopifnot(all.equal(m.t.7, 1.55326395545110, tolerance = 1.e-9)) ##REML value
stopifnot(all.equal(m.t.7, 1.56642013605506, tolerance = 1.e-6)) ## ML

hist.t <- cut(t.fm7, floor(min(t.fm7)) : ceiling(max(t.fm7)))
cbind(table(hist.t))

}# fixed effect 'd' -- with 'sparseX' only --------------------------------

if(testLevel <= 1) { cat('Time elapsed: ', proc.time(),'\n'); q("no") }

## ELSE : (testLevel > 1) :
library(lattice)
source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
##--> all.equal(), isOptimized(), ...

system.time(
    fm8.N <- lmer(y ~ service * dept + studage + lectage +
                  (1|s) + (1|d), InstEval, REML=FALSE,
                  control=lmerControl("Nelder_Mead"), verbose = 1L)
    ) ## 62   sec [MM@lynne; 2013-11]
      ## 59.5 sec [nb-mm3;   2013-12-31]
system.time(
    fm8.B <- lmer(y ~ service * dept + studage + lectage +
                  (1|s) + (1|d), InstEval, REML=FALSE,
                  control=lmerControl("bobyqa"), verbose = 2L)
    )
      ## 34.1 sec [nb-mm3; 2013-12-31]

stopifnot(isOptimized(fm8.N), isOptimized(fm8.B))
all.equal(fm8.B, fm8.N, tolerance=0)
## "Mean relative difference: 3.31 e-06"   [nb-mm3; 2013-12-31]

str(baseOpti(fm8.N))
str(baseOpti(fm8.B))

(sm8 <- summary(fm8.B))
str(r8 <- ranef(fm8.B))
sapply(r8, summary)
r.m8 <- cov2cor(vcov(sm8))
Matrix::image(r.m8, main="cor(<fixed eff[ lmer(*,InstEval) ]>)")

if(testLevel <= 2) { cat('Time elapsed: ', proc.time(),'\n'); q("no") }

## ELSE: testLevel > 2

## Clearly smaller X, but more RE pars
## ==> currently considerably slower than the above

system.time(
    fm9 <-
    lmer(y ~ studage + lectage +
         (1|s) + (1|d) + (1|dept:service) + (1|dept),
         InstEval, verbose = 1L, REML=FALSE)
    ) ## 410 secs [MM@lynne; 2013-11]

fm9
(sm9 <- summary(fm9))
rr <- ranef(fm9, condVar = TRUE) ## ~ 10 secs

sapply(rr, summary)
qqr <- qqmath(rr, strip=FALSE)
qqr$d
qqr$s
dotplot(rr,strip=FALSE)$`dept:service`

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

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
	  all.equal(fe1, fe1.s, tol= 1e-12
	  ),
	  all.equal(se1.d <- coef(s1.d)[,"Std. Error"],
		    se1.s <- coef(s1.s)[,"Std. Error"]#, tol = 1e-10
	  ),
	  all.equal(V.d <- vcov(fm1),
		    V.s <- vcov(fm1.s)#, tol = 1e-9
	  ),
	  all.equal(Matrix::diag(V.d), unname(se1.d)^2, tol= 1e-12)
	  ,
	  all.equal(unname(se1.d),
                    c(0.576011960125099, rep.int(0.518683987372292,3L)),
                    tol = 1.e-8)
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
#stopifnot(all.equal(m.t.7, 1.55326395545110, tol = 1.e-9)) ##REML value
stopifnot(all.equal(m.t.7, 1.56642013605506, tol = 1.e-6)) ## ML

hist.t <- cut(t.fm7, floor(min(t.fm7)) : ceiling(max(t.fm7)))
cbind(table(hist.t))

}# fixed effect 'd' -- with 'sparseX' only --------------------------------

if(testLevel <= 1) { cat('Time elapsed: ', proc.time(),'\n'); q("no") }

## ELSE : (testLevel > 1) :
library(lattice)

system.time(
    fm8 <- lmer(y ~ service * dept + studage + lectage +
                (1|s) + (1|d), InstEval, verbose = 1L, REML=FALSE)
    ) ## 62 secs [MM@lynne; 2013-11]
fm8
sm8 <- summary(fm8)
str(r8 <- ranef(fm8))
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

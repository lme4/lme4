require(lme4)
source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

## Check that quasi families throw an error
assertError(lmer(cbind(incidence, size - incidence) ~ period + (1|herd),
		 data = cbpp, family = quasibinomial))
assertError(lmer(incidence ~ period + (1|herd),
		 data = cbpp, family = quasipoisson))
assertError(lmer(incidence ~ period + (1|herd),
		 data = cbpp, family = quasi))

## check bug found by Kevin Buhr
set.seed(7)
n <- 10
X <- data.frame(y=runif(n), x=rnorm(n), z=sample(c("A","B"), n, TRUE))
fm <- lmer(log(y) ~ x | z, data=X)
## gave error inside  model.frame()
stopifnot(all.equal(unname(fixef(fm)), -0.8345, tol=.01))

sstudy9 <- subset(sleepstudy, Days == 1 | Days == 9)
try({## This "did work" in lme4.0 and nlme -- FIXME ??
 m1 <- lmer(Reaction ~ 1 + Days + (1 + Days | Subject), data = sstudy9)
 if (FALSE) {
     image(reTrms$Zt)
     testRankZt <- function(object) {
         rankMatrix(cBind(reTrms$Zt,1))<nrow(reTrms$Zt)
     }
 }
 ## -> Error in ptr() : Downdated VtV is not positive definite
 ## FIXME?(2): More helpful error message, or rank-checking diagnostics?
 ## how would we check the rank here?

 print(sm1 <- summary(m1))
 fm1 <- fitted(m1)
})
showProc.time()

## check working of Matrix methods on  vcov(.) etc ----------------------
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
V  <- vcov(fm)
V1 <- vcov(fm1)
TOL <- 0 # to show the differences below
TOL <- 1e-5 # for the check
stopifnot(
	  all.equal(diag(V), 0.176078, tol = TOL) # 64b: 2.4e-8
	  ,
	  all.equal(as.numeric(chol(V)), 0.4196165, tol = TOL)	# 64b: 3.2e-8
	  ,
	  all.equal(diag(V1), c(46.574978, 2.389469), tol = TOL)# 64b: 9.8e-9
	  , dim(C1 <- chol(V1)) == c(2,2) ,
	  all.equal(as.numeric(C1),
		    c(6.82458627, 0, -0.2126260, 1.5310973), tol=TOL)# 64b: 1.6e-9
          ,
          dim(chol(crossprod(getME(fm1, "Z")))) == 36
	  , TRUE)
## printing
signif(chol(crossprod(getME(fm,"Z"))), 4)# -> simple 4 x 4 sparse

showProc.time() #

## From: Stephane Laurent
## To:   r-sig-mixed-models@..
## "crash with the latest update of lme4"
##
## .. example for which lmer() crashes with the last update of lme4 ...{R-forge},
## .. but not with version CRAN version (0.999999-0)
lsDat <- data.frame(
                  Operator = as.factor(rep(1:5, c(3,4,8,8,8))),
                  Part = as.factor(
                  c(2L, 3L, 5L,
                    1L, 1L, 2L, 3L,
                    1L, 1L, 2L, 2L, 3L, 3L, 4L, 5L,
                    1L, 2L, 3L, 3L, 4L, 4L, 5L, 5L,
                    1L, 2L, 2L, 3L, 3L, 4L, 5L, 5L)),
                  y =
                  c(0.34, -1.23, -2.46,
                    -0.84, -1.57,-0.31, -0.18,
                    -0.94, -0.81, 0.77, 0.4, -2.37, -2.78, 1.29, -0.95,
                    -1.58, -2.06, -3.11,-3.2, -0.1, -0.49,-2.02, -0.75,
                    1.71,  -0.85, -1.19, 0.13, 1.35, 1.92, 1.04,  1.08))
xtabs( ~ Operator + Part, data=lsDat) # --> 4 empty cells, quite a few with only one obs.:
##         Part
## Operator 1 2 3 4 5
##        1 0 1 1 0 1
##        2 2 1 1 0 0
##        3 2 2 2 1 1
##        4 1 1 2 2 2
##        5 1 2 2 1 2

old_lme4_vcov <- 0.254166752231642 ## lme4.0_0.999999-1

fm3 <- suppressWarnings(lmer(y ~ (1|Part) + (1|Operator) + (1|Part:Operator),
            data = lsDat))

## FIXME: should suppress 'caught warning' message in optwrap!
## --> *many* warnings   Cholmod warning 'not positive definite' ..

## FIXME: failing in vcov.merMod <- mkVcov <- as(V,"dpoMatrix")
## (fm3@pp$unsc() contains NaN, function catches "not a positive definite matrix"
##  and converts to "Computed variance-covariance matrix is not positive definite")
## fm3 # gave an error {no longer does}
## currently (2012-06-25) full of  NA/NaN  -- old lme4 does
showProc.time()

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

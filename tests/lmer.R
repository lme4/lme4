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
 ## -> Error in ptr() : Downdated VtV is not positive definite
 ## FIXME?(2): More helpful error message
 print(sm1 <- summary(m1))
 fm1 <- fitted(m1)
})

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



cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

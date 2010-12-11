require(lme4)
source(system.file("test-tools.R", package = "Matrix"))# identical3() etc

## Check that quasi families throw an error
assertError(lmer(cbind(incidence, size - incidence) ~ period + (1|herd),
		 data = cbpp, family = quasibinomial))
assertError(lmer(incidence ~ period + (1|herd),
		 data = cbpp, family = quasipoisson))
assertError(lmer(incidence ~ period + (1|herd),
		 data = cbpp, family = quasi))

## check bug found by Kevin Buhner
set.seed(7)
n <- 10
X <- data.frame(y=runif(n), x=rnorm(n), z=sample(c("A","B"), n, TRUE))
fm <- lmer(log(y) ~ x | z, data=X)
## gave error inside  model.frame()
stopifnot(all.equal(unname(fixef(fm)), -0.8345, tol=.01))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

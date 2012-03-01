## generalized linear mixed model
stopifnot(suppressPackageStartupMessages(require(lme4Eigen)))
options(show.signif.stars = FALSE)
source(system.file("test-tools.R", package = "Matrix"))# identical3(), showProc.time, etc


## TODO: (1) move these to ./glmer-ex.R [DONE]
## ----  (2) "rationalize" with ../man/cbpp.Rd
#m1e <- glmer1(cbind(incidence, size - incidence) ~ period + (1 | herd),
#              family = binomial, data = cbpp, doFit = FALSE)
## now
#bobyqa(m1e, control = list(iprint = 2L))

m0 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = binomial, data = cbpp, verbose = 2L)

m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            optimizer="bobyqa",
            family = binomial, data = cbpp, verbose = 2L,
            control = list(rhobeg=0.2, rhoend=2e-7), tolPwrss=1e-8)

m2 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            family = binomial, data = cbpp, nAGQ=25L)

## loosened tolerance on parameters
stopifnot(is((cm1 <- coef(m2)), "coef.mer"),
	  dim(cm1$herd) == c(15,4),
	  all.equal(fixef(m2),
### lme4a [from an Ubuntu 11.10 amd64 system]
                    ### c(-1.39922533406847, -0.991407294757321,
                    ###   -1.12782184600404, -1.57946627431248),
                    c(-1.3766013, -1.0058773,
                      -1.1430128, -1.5922817),
		    tol = 5.e-4,
                    check.attr=FALSE),
##        all.equal(deviance(m2), 100.010030538022, tol=1e-9)
          all.equal(deviance(m2), 101.11977669, tol=1e-9)
)

stopifnot(is((cm1 <- coef(m1)), "coef.mer"),
	  dim(cm1$herd) == c(15,4),
	  all.equal(fixef(m1),
                    ##  these values are those of "old-lme4":
		    ## c(-1.39853504914, -0.992334711,
		    ##   -1.12867541477, -1.58037390498),
                    ## lme4Eigen[r 1636], 64-bit ubuntu 11.10:
                    c(-1.3788385, -1.0589543,
                      -1.1936382, -1.6306271),
		    tol = 1.e-3,
                    check.attr=FALSE)
	  )


## Deviance for the new algorithm is lower, eventually we should change the previous test
##stopifnot(deviance(m1) <= deviance(m1e))

showProc.time() #

## FIXME -- non-convegence!!
if (FALSE && require('MASS', quietly = TRUE)) {
    bacteria$wk2 <- bacteria$week > 2
    contrasts(bacteria$trt) <-
        structure(contr.sdif(3),
                  dimnames = list(NULL, c("diag", "encourage")))
    print(fm5 <- glmer(y ~ trt + wk2 + (1|ID), bacteria, binomial))
    ## again *fails* (lme4Eigen[r 1636], 64-bit ubuntu 11.10)
    ## used to fail with nlminb() : stuck at theta=1

    showProc.time() #

    stopifnot(
	      all.equal(logLik(fm5),
			## was	  -96.127838
			structure(-96.13069, nobs = 220L, nall = 220L,
				  df = 5L, REML = FALSE,
                                  class = "logLik"),
                        tol = 1e-5, check.attributes = FALSE)
	      ,
	      all.equal(fixef(fm5),
			## was		 2.834218798		 -1.367099481
			c("(Intercept)"= 2.831609490, "trtdiag"= -1.366722631,
			  ## now	 0.5842291915,		 -1.599148773
			  "trtencourage"=0.5840147802, "wk2TRUE"=-1.598591346), tol = 1e-4)
	      )
}

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
set.seed(1)
dd <- rPoisGLMMi(12, 20)
m0  <- glmer(y~x + (1|f),           family="poisson", data=dd)
(m1 <- glmer(y~x + (1|f) + (1|obs), family="poisson", data=dd))
(a01 <- anova(m0, m1))

stopifnot(all.equal(a01$Chisq[2], 554.334056, tol=1e-6),
	  all.equal(a01$logLik, c(-1073.77193, -796.604902), tol=1e-6),
          a01$ Df == 3:4,
	  a01$`Chi Df`[2] == 1)


showProc.time()

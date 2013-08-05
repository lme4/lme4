####  Utility Functions for lme4 Testing
####  ----------------------------------
if(FALSE) ### "Load" these by
    source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
## e.g. from ../../tests/glmmWeights.R

##' General GLMM simulation  for the model   y ~ x + (1|block)
##' i.e. (2-way notation)   Y_{ij} = \beta_1 + \beta_2 x_{ij} +  B_j,
##'                                  i = 1..nperblk,  j = 1..nblk
gSim <- function(nblk=26,
                 nperblk=100,
                 sigma=1,    ## st.dev. of ran.eff  B
                 beta=c(4,3),
                 x = runif(n),## sd = sqrt(1/12) = 0.2887; allow also x=c(0,1)
                 shape=2,    ## shape parameter for Gamma
                 nbinom=10,  ## N for binomial trials
                 family=Gamma())
{
    stopifnot(nblk <= 50000)# some sanity, may increase (but remain "finite"!)
    ## ch.set: a potentially large set of "letters", as level-labels for 'block':
    nc <- length(ch.set <- c(LETTERS, letters,
                             paste0(LETTERS,LETTERS), paste0(LETTERS,letters)))
    while(nblk > nc) nc <-
	length(ch.set <- c(paste0(ch.set, LETTERS),
			   paste0(ch.set, letters)))
    stopifnot(1 <= nblk, nblk <= length(ch.set))
    d <- expand.grid(block = ch.set[1:nblk],
                     rep= 1:nperblk, KEEP.OUT.ATTRS=FALSE)
    stopifnot(nblk == length(levels(d$block)),
              (n <- nrow(d)) == nblk * nperblk, length(x) == n, length(beta) == 2)
    d$ x <- x
    reff_f <- rnorm(nblk, sd=sigma)
    ## need intercept large enough to avoid negative values
    d$ eta0 <- beta[1] + beta[2]*x  ## fixed effects only
    d$ eta  <- d$eta0 + reff_f[d$block]
    d$ mu   <- family$linkinv(d$eta)
    d$ y <- switch(family$family,
                   "Gamma" = rgamma(n,scale=d$mu/shape,shape=shape),
                   "poisson" = rpois(n,d$mu),
		   "binomial"= {
		       z <- rbinom(n, prob=d$mu, size=nbinom)
		       if (nbinom==1) z else cbind(succ = z, fail = nbinom-z)
		   },
                   stop("Family ", family$family, " not supported here"))
    d
}

.simpleCap <- function(x) {
  paste(toupper(substring(x, 1,1)), substring(x, 2),
        sep="", collapse=" ")
}

### bootMer() --- <==>	(TODO: semi-*)parametric bootstrap
### -------------------------------------------------------
## Doc: show how  this is equivalent - but faster than
##		boot(*, R = nsim, sim = "parametric", ran.gen = simulate(.,1,.), mle = x)
## --> return a "boot" object -- so we could use boot.ci() etc
## TODO: also allow "semi-parametric" model-based bootstrap:
##    resampling the (centered!) residuals (for use.u=TRUE) or for use.u=FALSE,
##    *both* the centered u's + centered residuals
##    instead of using	rnorm()

##' Model-based (Semi-)Parametric Bootstrap for Mixed Models
##' 
##' Perform model-based (Semi-)parametric bootstrap for mixed models.
##' 
##' The semi-parametric variant is not yet implemented, and we only
##' provide a method for \code{\link{lmer}}  and \code{\link{glmer}} results.
##' 
##' The working name for bootMer() was \dQuote{simulestimate()}, as it is an
##' extension of \code{\link{simulate}}, but we rather want to emphasize its
##' potential for valid inference.
##' 
##' In each of the \code{nsim} simulations --- that is what the
##' \emph{parametric} bootstrap does, both \dQuote{\emph{spherized}} random
##' effects \eqn{u} and the i.i.d errors \eqn{\epsilon} are generated, using
##' \code{\link{rnorm}()} with parameters corresponding to the fitted model
##' \code{x}.
##' 
##' @param x fitted \code{*lmer()} model, see \code{\link{lmer}},
##'     \code{\link{glmer}}, etc.
##' @param FUN a \code{\link{function}(x)}, computating the \emph{statistic} of
##'     interest, which must be a numeric vector, possibly named.
##' @param nsim number of simulations, positive integer; the bootstrap \eqn{B}
##'     (or \eqn{R}).
##' @param seed optional argument to \code{\link{set.seed}}.
##' @param use.u logical, indicating, if the spherized random effects should be
##'     simulated / bootstrapped as well.  If \code{FALSE}, they are not changed,
##'     and all inference is conditional on these.
##' @param verbose logical indicating if progress should print output
##' @param control an optional \code{\link{list}}, to be passed to the minimizer
##'     (of the log-likelihood, or RE likelihood).
##' @return an object of S3 \code{\link{class}} \code{"boot"}, compatible with
##'     \pkg{boot} package's \code{boot()} result.
##' @seealso For inference, including confidence intervals,
##'     \code{\link{profile-methods}}.
##' 
##'     \code{\link[boot]{boot}()}, and then \code{\link[boot]{boot.ci}} from
##'     package \pkg{boot}.
##' @references Davison, A.C. and Hinkley, D.V. (1997) \emph{Bootstrap Methods
##'     and Their Application}.  Cambridge University Press.
##' @keywords models htest
##' @examples
##' fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
##' ## see ?"profile-methods"
##' mySumm <- function(.) { s <- sigma(.)
##'     c(beta =getME(., "beta"), sigma = s, sig01 = s * getME(., "theta")) }
##' (t0 <- mySumm(fm01ML)) # just three parameters
##' 
##' \dontrun{%%--- fails for now --- FIXME
##' 
##' ## 3.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
##' system.time( boo01 <- bootMer(fm01ML, mySumm, nsim = 100) )
##' 
##' ## to "look" at it
##' if(need.boot <- is.na(match("package:boot", search())))
##'   require("boot")## is a recommended package, i.e. *must* be there
##' boo01 # look at estimated bias for sig01 (-9.1, also when nsim = 1000)
##' 
##' ## ------ Bootstrap-based confidence intervals ------------
##' 
##' (bCI.1 <- boot.ci(boo01, index=1, type=c("norm", "basic", "perc")))# beta
##' 
##' ## Sigma - original scale, and interval calculation on log-scale:
##' (bCI.2  <- boot.ci(boo01, index=2, type=c("norm", "basic", "perc")))
##' (bCI.2l <- boot.ci(boo01, index=2, type=c("norm", "basic", "perc"),
##'                    h = log, hdot = function(.) 1/., hinv = exp))
##' 
##' (bCI.3 <- boot.ci(boo01, index=3, type=c("norm", "basic", "perc")))# sig01
##' }
##' @export
bootMer2 <- function(x, FUN, nsim = 1, seed = NULL, use.u = FALSE,
                    type=c("parametric","semiparametric"),
		    verbose = FALSE, control = list(),
                    .progress="none", PBargs=list()) {
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  if (.progress!="none") {  ## progress bar
    pbfun <- get(paste(.progress,"ProgressBar",sep=""))
    setpbfun <- get(paste("set",.simpleCap(.progress),"ProgressBar",sep=""))
    pb <- do.call(pbfun,PBargs)
  }
    FUN <- match.fun(FUN)
    type <- match.arg(type)
    if(!is.null(seed)) set.seed(seed)
    else if(!exists(".Random.seed", envir = .GlobalEnv))
	runif(1) # initialize the RNG if necessary

    mc <- match.call()
    t0 <- FUN(x)
    if (!is.numeric(t0))
	stop("bootMer currently only handles functions that return numeric vectors")

    mle <- list(beta = getME(x,"beta"), theta = getME(x,"theta"))
    if (lme4Eigen:::isLMM(x)) mle <- c(mle,list(sigma = sigma(x)))
    ## FIXME: what about GLMMs with scale parameters??
    ## FIXME: remove prefix when incorporated in package

    if (type=="parametric") {
      ss <- simulate(x, nsim=nsim, use.u=use.u)
    } else {
      if (use.u) {
        ## FIXME: does this work for GLMMs???
        ss <- replicate(nsim,fitted(x)+sample(residuals(x,"response")),
                        simplify=FALSE)
      } else {
        stop("not implemented")
      }
    }
    t.star <- matrix(t0, nrow = length(t0), ncol = nsim)
    for(i in 1:nsim) {
      if (.progress!="none") { setpbfun(pb,i/nsim) }
      foo <- try(FUN(refit(x,ss[[i]])))
      if(verbose) { cat(sprintf("%5d :",i)); str(foo) }
      t.star[,i] <- if (inherits(foo, "error")) NA else foo
    }
  if (.progress!="none") { close(pb) }
    rownames(t.star) <- names(t0)

    ## boot() ends with the equivalent of
    ## structure(list(t0 = t0, t = t.star, R = R, data = data, seed = seed,
    ##		      statistic = statistic, sim = sim, call = call,
    ##		      ran.gen = ran.gen, mle = mle),
    ##		 class = "boot")
    structure(list(t0 = t0, t = t(t.star), R = nsim, data = x@frame,
		   seed = .Random.seed,
		   statistic = FUN, sim = "parametric", call = mc,
		   ## these two are dummies
		   ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
	      class = "boot")
}## {bootMer}

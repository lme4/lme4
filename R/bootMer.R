.simpleCap <- function(x) {
  paste0(toupper(substr(x, 1,1)), substr(x, 2, 1000000L), collapse=" ")
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
##  BUT see:
## @article{morris_blups_2002,
##	title = {The {BLUPs} are not "best" when it comes to bootstrapping},
##	volume = {56},
##	issn = {0167-7152},
##	url = {http://www.sciencedirect.com/science/article/pii/S016771520200041X},
##	doi = {10.1016/S0167-7152(02)00041-X},
##	journal = {Statistics \& Probability Letters},
##	author = {Morris, Jeffrey S},
##	year = {2002},
## }
## for an indication of why this is not necessarily a good idea!

##' Model-based (Semi-)Parametric Bootstrap for Mixed Models
##'
##' Perform model-based (Semi-)parametric bootstrap for mixed models.
##'
##' The semi-parametric variant is not yet implemented, and we only
##' provide a method for \code{\link{lmer}}  and \code{\link{glmer}} results.
##'
##' The working name for bootMer() was \dQuote{simulestimate()}, as it is an
##' extension of \code{\link{simulate}}, but we want to emphasize its
##' potential for valid inference.
##'
##' @param x fitted \code{*lmer()} model, see \code{\link{lmer}},
##'     \code{\link{glmer}}, etc.
##' @param FUN a \code{\link{function}(x)}, computating the \emph{statistic} of
##'     interest, which must be a numeric vector, possibly named.
##' @param nsim number of simulations, positive integer; the bootstrap \eqn{B}
##'     (or \eqn{R}).
##' @param seed optional argument to \code{\link{set.seed}}.
##' @param use.u logical, indicating whether the spherical random effects
##' should be simulated / bootstrapped as well.  If \code{TRUE}, they are not
##' changed, and all inference is conditional on these values. If \code{FALSE},
##' new normal deviates are drawn (see Details).
##' @param type character string specifying the type of bootstrap,
##'    \code{"parametric"} or \code{"semiparametric"}; partial matching is allowed.
##' @param verbose logical indicating if progress should print output
##' @param .progress character string - type of progress bar to
##'     display.  Default is \code{"none"}; the function will look for a relevant \code{*ProgressBar} function, so \code{"txt"} will work in general; \code{"tk"} is available if the \code{tcltk} package is loaded; or \code{"win"} on Windows systems.
##' @param PBargs a list of additional arguments to the progress bar function (the package authors like \code{list(style=3)}).
##' @param parallel The type of parallel operation to be used (if any).  If missing, the
##' default is taken from the option \code{"boot.parallel"} (and if that
##' is not set, \code{"no"}).
##' @param ncpus integer: number of processes to be used in parallel operation:
##' typically one would chose this to be the number of available CPUs.
##' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use if
##' \code{parallel = "snow"}.  If not supplied, a cluster on the
##' local machine is created for the duration of the \code{boot} call.

##' @return an object of S3 \code{\link{class}} \code{"boot"}, compatible with
##'     \pkg{boot} package's \code{\link[boot]{boot}()} result.
##' @seealso
##' \itemize{
##' \item For inference, including confidence intervals,
##'     \code{\link{profile-methods}}.
##' \item \code{\link[boot]{boot}()}, especially for information about parallel operation, and then \code{\link[boot]{boot.ci}} from
##'     package \pkg{boot}.
##' @author Martin Maechler and Ben Bolker; parallel interface adapted from code by Brian Ripley
##' }
##' @details
##' \itemize{
##' \item If \code{use.u} is \code{FALSE} and \code{type} is \code{"parametric"}, each simulation generates
##' new values of both the \dQuote{\emph{spherical}} random
##' effects \eqn{u} and the i.i.d. errors \eqn{\epsilon}, using
##' \code{\link{rnorm}()} with parameters corresponding to the fitted model
##' \code{x}.
##' \item If \code{use.u} is \code{TRUE} and \code{type=="parametric"}, only the i.i.d. errors
##' (or, for GLMMs, response values drawn from the appropriate distributions)
##' are resampled, with the values of \eqn{u} staying fixed at their
##' estimated values.
##' \item If \code{use.u} is \code{TRUE} and \code{type=="semiparametric"},
##' the i.i.d. errors are sampled from the distribution of (response) residuals.
##' (For GLMMs, the resulting sample will no longer have the same
##' properties as the original sample, and the method may not make sense;
##' a warning is generated.)
##' The semiparametric bootstrap is currently an experimental feature, and therefore
##' may not be stable.
##' \item The case where \code{use.u} is \code{FALSE} and \code{type=="semiparametric"} is not implemented; Morris (2002) suggests that resampling from the estimated values of \eqn{u} is not good practice.
##' }
##' @references
##' Davison, A.C. and Hinkley, D.V. (1997)
##' \emph{Bootstrap Methods and Their Application}.  Cambridge University Press.
##'
##' Morris, J. S. (2002).
##' The BLUPs Are Not \sQuote{best} When It Comes to Bootstrapping.
##' \emph{Statistics & Probability Letters} \bold{56}(4): 425--430.
##' doi:10.1016/S0167-7152(02)00041-X.
##' @keywords models htest
##' @examples
##' fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
##' ## see ?"profile-methods"
##' mySumm <- function(.) { s <- sigma(.)
##'     c(beta =getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) }
##' (t0 <- mySumm(fm01ML)) # just three parameters
##' ## alternatively:
##' mySumm2 <- function(.) {
##' c(beta=fixef(.),sigma=sigma(.),sig01=unlist(VarCorr(.)))
##' }
##'
##' set.seed(101)
##' ## 3.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
##' system.time( boo01 <- bootMer(fm01ML, mySumm, nsim = 100) )
##'
##' ## to "look" at it
##' require("boot") ## a recommended package, i.e. *must* be there
##' boo01
##' ## note large estimated bias for sig01
##' ## (~30% low, decreases _slightly_ for nsim = 1000)
##'
##' ## extract the bootstrapped values as a data frame ...
##' head(as.data.frame(boo01))
##' 
##' ## ------ Bootstrap-based confidence intervals ------------
##'
##' ## intercept
##' (bCI.1 <- boot.ci(boo01, index=1, type=c("norm", "basic", "perc")))# beta
##'
##' ## Residual standard deviation - original scale:
##' (bCI.2  <- boot.ci(boo01, index=2, type=c("norm", "basic", "perc")))
##' ## Residual SD - transform to log scale:
##' (bCI.2l <- boot.ci(boo01, index=2, type=c("norm", "basic", "perc"),
##'                    h = log, hdot = function(.) 1/., hinv = exp))
##'
##' ## Among-batch variance:
##' (bCI.3 <- boot.ci(boo01, index=3, type=c("norm", "basic", "perc")))# sig01
##'
##' ## Graphical examination:
##' plot(boo01,index=3)
##'
##' ## Check stored values from a longer (1000-replicate) run:
##' load(system.file("testdata","boo01L.RData",package="lme4"))
##' plot(boo01L,index=3)
##' 
##' @export
bootMer <- function(x, FUN, nsim = 1, seed = NULL, use.u = FALSE,
                    type=c("parametric","semiparametric"),
		    verbose = FALSE,
                    .progress="none", PBargs=list(),
                    parallel = c("no", "multicore", "snow"),
                    ncpus = getOption("boot.ncpus", 1L), cl = NULL)
{
    stopifnot((nsim <- as.integer(nsim[1])) > 0)
    if (.progress!="none") { ## progress bar
        pbfun <- get(paste0(.progress,"ProgressBar"))
        setpbfun <- get(paste0("set",.simpleCap(.progress),"ProgressBar"))
        pb <- do.call(pbfun,PBargs)
    }
    if (missing(parallel)) parallel <- getOption("boot.parallel", "no")
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    if (parallel != "no" && ncpus > 1L) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
        if (!have_mc && !have_snow) ncpus <- 1L
    }
    do_parallel <- (ncpus > 1L && (have_mc || have_snow))
    if (do_parallel & .progress!="none")
        message("progress bar disabled for parallel operations")

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
    if (isLMM(x)) mle <- c(mle,list(sigma = sigma(x)))
    ## FIXME: what about GLMMs with scale parameters??
    ## FIXME: remove prefix when incorporated in package

    if (type=="parametric") {
        ss <- simulate(x, nsim=nsim, use.u=use.u, na.action=na.exclude)
    } else {
        if (use.u) {
            if (isGLMM(x)) warning("semiparametric bootstrapping is questionable for GLMMs")
            ss <- replicate(nsim,fitted(x)+sample(residuals(x,"response")),
                            simplify=FALSE)
        } else {
            stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
        }
    }
    
    # define ffun as a closure containing the referenced variables 
    # in its scope to avoid explicit clusterExport statement
    # in the PSOCKcluster case 
    ffun <- local({
      FUN 
      refit 
      x 
      ss 
      verbose 
      do_parallel
      length.t0 <- length(t0)
      function(i) {
        foo <- try(FUN(refit(x,ss[[i]])),silent=TRUE)
        if(verbose) { cat(sprintf("%5d :",i)); str(foo) }
        if (!do_parallel && .progress!="none") { setpbfun(pb,i/nsim) }
        if (inherits(foo, "try-error")) rep(NA, length.t0) else foo
    }})

    simvec <- seq_len(nsim)
     res <- if (do_parallel) {
        if (have_mc) {
            parallel::mclapply(simvec, ffun, mc.cores = ncpus)
        } else if (have_snow) {
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                ## explicit export of the lme4 namespace since most FUNs will probably
                ## use some of them
                parallel::clusterExport(cl, varlist=getNamespaceExports("lme4"))
                if(RNGkind()[1L] == "L'Ecuyer-CMRG")
                    parallel::clusterSetRNGStream(cl)
                res <- parallel::parLapply(cl, simvec, ffun)
                parallel::stopCluster(cl)
                res
            } else parallel::parLapply(cl, simvec, ffun)
        }
    } else lapply(simvec, ffun)

    t.star <- do.call(cbind,res)
    rownames(t.star) <- names(t0)
    if ((numFail <- sum(apply(is.na(t.star),2,all)))>0) {
        warning("some bootstrap runs failed (",numFail,"/",nsim,")")
    }
    ## boot() ends with the equivalent of
    ## structure(list(t0 = t0, t = t.star, R = R, data = data, seed = seed,
    ##		      statistic = statistic, sim = sim, call = call,
    ##		      ran.gen = ran.gen, mle = mle),
    ##		 class = "boot")
    s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = x@frame,
		   seed = .Random.seed,
		   statistic = FUN, sim = "parametric", call = mc,
		   ## these two are dummies
		   ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
	      class = "boot")
    attr(s,"bootFail") <- numFail
    s
} ## {bootMer}

##' @S3method as.data.frame boot
as.data.frame.boot <- function(x,...) {
  as.data.frame(x$t)
}


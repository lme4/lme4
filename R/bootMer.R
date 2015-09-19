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

##' Perform model-based (Semi-)parametric bootstrap for mixed models.
##'
##' The semi-parametric variant is not yet implemented, and we only
##' provide a method for \code{\link{lmer}}  and \code{\link{glmer}} results.
bootMer <- function(x, FUN, nsim = 1, seed = NULL,
                    use.u = FALSE, re.form = NA,
                    type = c("parametric","semiparametric"),
		    verbose = FALSE,
                    .progress = "none", PBargs=list(),
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
    do_parallel <- (parallel != "no" && ncpus > 1L)
    if (do_parallel) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
	if (!(have_mc || have_snow))
	    do_parallel <- FALSE # (only for "windows")
    }
    if (do_parallel && .progress != "none")
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
        argList <- list(x, nsim=nsim, na.action=na.exclude)
        if (!missing(re.form)) {
            argList <- c(argList,list(re.form=re.form))
        } else {
            argList <- c(argList,list(use.u=use.u))
        }
        ss <- do.call(simulate,argList)
    } else {
        if (!missing(re.form))
            stop(paste(sQuote("re.form")),
                 "cannot be used with semiparametric bootstrapping")
        if (use.u) {
            if (isGLMM(x)) warning("semiparametric bootstrapping is questionable for GLMMs")
            ss <- replicate(nsim,fitted(x)+sample(residuals(x,"response"),
                                                  replace=TRUE),
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
        ret <- tryCatch(FUN(refit(x,ss[[i]])), error=function(e)e)
        if(verbose) { cat(sprintf("%5d :",i)); str(ret) }
        if (!do_parallel && .progress!="none") { setpbfun(pb,i/nsim) }
        if (inherits(ret, "error"))
            structure(rep(NA, length.t0), "fail.msgs" = ret$message)
        else
            ret
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
    if ((numFail <- sum(bad.runs <- apply(is.na(t.star),2,all)))>0) {
        warning("some bootstrap runs failed (",numFail,"/",nsim,")")
        fail.msgs <- vapply(res[bad.runs],FUN=attr,FUN.VALUE=character(1),
                            "fail.msgs")
    } else fail.msgs <- NULL
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
    attr(s,"boot.fail.msgs") <- fail.msgs
    s
} ## {bootMer}

##' @S3method as.data.frame boot
as.data.frame.boot <- function(x,...) {
  as.data.frame(x$t)
}


lmer <- function(formula, data, REML = TRUE, sparseX = FALSE,
		  control = list(), start = NULL,
		  verbose = 0L, doFit = TRUE,
		  subset, weights, na.action, offset,
		  contrasts = NULL, devFunOnly=FALSE, ...)
{
    mf <- mc <- match.call()
    ## '...' handling up front, safe-guarding against typos ("familiy") :
    if(length(l... <- list(...))) {
	if (!is.null(l...$family)) {  # call glmer if family specified
	    mc[[1]] <- as.name("glmer")
	    return( eval(mc, parent.frame()) )
	}
	## Check for method argument which is no longer used
	if (!is.null(method <- l...$method)) {
	    msg <- paste("Argument", sQuote("method"), "is deprecated.")
	    if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
		warning(msg)
		l... <- l...[names(l...) != "method"]
	    } else stop(msg)
	}
	if(length(l...))
	    warning("extra arguments ", paste(names(l...), sep=", "),
		    " are disregarded")
    }

    stopifnot(length(formula <- as.formula(formula)) == 3)
    if (missing(data)) data <- environment(formula)
					# evaluate and install the model frame :
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
	       names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substituted "|" by "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
					# random effects and terms modules
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    if (any(unlist(lapply(reTrms$flist, nlevels)) >= nrow(fr)))
	stop("number of levels of each grouping factor must be less than number of obs")
    ## fixed-effects model matrix X - remove random parts from formula:
    form <- formula
    form[[3]] <- if(is.null(nb <- nobars(form[[3]]))) 1 else nb
    X <- model.Matrix(form, fr, contrasts, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    p <- ncol(X)
    pp <- new("merPredD", X=X, Zt=reTrms$Zt, Lambdat=reTrms$Lambdat, Lind=reTrms$Lind,
	      theta=reTrms$theta + 1e-06) # ensure all entries of
					# theta are nonzero because
					# Eigen tends to prune sparse matrices
                                        # (Probably not needed any more).
    resp <- mkRespMod2(fr)
    if (REML) resp$REML <- p

    devfun <- mkdevfun(pp, resp)
    if (devFunOnly) return(devfun)
					# one evaluation to ensure all values are set
    opt <- list(fval=devfun(reTrms$theta))

    if (doFit) {			# optimize estimates
	if (verbose) control$iprint <- as.integer(verbose)
	opt <- bobyqa(reTrms$theta, devfun, reTrms$lower, control = control)
    }
    LL <- pp$Lambdat
    LL@x <- pp$theta[pp$Lind]
    stopifnot(validObject(LL))
    pp$Lambdat <- LL
    sqrLenU <- pp$sqrL(1.)
    wrss <- resp$wrss()
    pwrss <- wrss + sqrLenU
    n <- nrow(fr)

    dims <- c(N=n, n=n, nmp=n-p, nth=length(pp$theta), p=p, q=nrow(reTrms$Zt),
	      nAGQ=NA_integer_, useSc=1L, reTrms=length(reTrms$cnms),
	      spFe=0L, REML=resp$REML, GLMM=0L, NLMM=0L)
    cmp <- c(ldL2=pp$ldL2(), ldRX2=pp$ldRX2(), wrss=wrss,
	      ussq=sqrLenU, pwrss=pwrss,
	     drsum=NA, dev=if(REML)NA else opt$fval, REML=if(REML)opt$fval else NA,
	     sigmaML=sqrt(pwrss/n), sigmaREML=sqrt(pwrss/(n-p)))

    new("lmerMod", call=mc, frame=fr, flist=reTrms$flist, cnms=reTrms$cnms,
	theta=pp$theta, beta=pp$delb(), u=pp$delu(), lower=reTrms$lower,
	devcomp=list(cmp=cmp, dims=dims), pp=pp, resp=resp)
}## { lmer }

mkdevfun <- function(pp, resp) {
    if (is(resp, "lmerResp"))
	return (function(theta) .Call(lmerDeviance, pp$ptr, resp$ptr, theta))
    stop("unknown response type: ", class(resp))
}


glmer <- function(formula, data, family = gaussian, sparseX = FALSE,
                  control = list(), start = NULL, verbose = 0L, nAGQ = 1L,
                  doFit = TRUE, compDev=TRUE, subset, weights, na.action, offset,
                  contrasts = NULL, mustart, etastart, devFunOnly = FALSE,
                  tolPwrss = 0.000001, ...)
{
    verbose <- as.integer(verbose)
    mf <- mc <- match.call()
    if (missing(family)) { ## divert using lmer()
	mc[[1]] <- as.name("lmer")
	return(eval(mc, parent.frame()))
    }
### '...' handling up front, safe-guarding against typos ("familiy") :
    if(length(l... <- list(...))) {
	## Check for invalid specifications
	if (!is.null(method <- list(...)$method)) {
	    msg <- paste("Argument", sQuote("method"),
			 "is deprecated.\nUse", sQuote("nAGQ"),
			 "to choose AGQ.  PQL is not available.")
	    if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
		warning(msg)
		l... <- l...[names(l...) != "method"]
	    } else stop(msg)
	}
	if(length(l...))
	    warning("extra arguments ", paste(names(l...), sep=", "),
		    " are disregarded")
    }
    if(is.character(family))
	family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    if (family$family %in% c("quasibinomial", "quasipoisson", "quasi"))
	stop('"quasi" families cannot be used in glmer')
    nAGQ <- as.integer(nAGQ)[1]
    if (nAGQ > 1L) warning("nAGQ > 1 has not been implemented, using Laplace")
    stopifnot(length(formula <- as.formula(formula)) == 3)
    if (missing(data)) data <- environment(formula)
					# evaluate and install the model frame :
    m <- match(c("data", "subset", "weights", "na.action", "offset",
		 "mustart", "etastart"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substitute "|" for "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
					# random-effects module
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    ## fixed-effects model matrix X - remove random parts from formula:
    form <- formula
    form[[3]] <- if(is.null(nb <- nobars(form[[3]]))) 1 else nb
    X <- model.Matrix(form, fr, contrasts, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    p <- ncol(X)
    pp <- merPredD$new(X=X, Zt=reTrms$Zt, Lambdat=reTrms$Lambdat, Lind=reTrms$Lind, theta=reTrms$theta)
					# response module
    resp <- mkRespMod2(fr, family=family)
					# initial step from working response
    if (compDev) {
	.Call(glmerWrkIter, pp$ptr, resp$ptr)
	lapply(1:3, function(n).Call(glmerPwrssUpdate, pp$ptr, resp$ptr, verbose, FALSE, tolPwrss))
    } else {
        pp$updateXwts(resp$sqrtWrkWt())
        pp$updateDecomp()
        pp$updateRes(resp$wrkResp())
        pp$solve()
        resp$updateMu(pp$linPred(1))	# full increment
        resp$updateWts()
        pp$installPars(1)
        lapply(1:3, function(n) pwrssUpdate(pp, resp, verbose, tol=tolPwrss))
    }

    u0 <- pp$u0
    beta0 <- pp$beta0

    opt <- list(fval=resp$Laplace(pp$ldL2(), pp$ldRX2(), pp$sqrL(0.)))

    if (doFit || devFunOnly) {			# optimize estimates
        rho <- as.environment(list(u0=pp$u0, beta0=pp$beta0, pp=pp, resp=resp,
                                   verbose=verbose, control=control, tolPwrss=tolPwrss))
        parent.env(rho) <- parent.frame()
        devfun <- if (compDev) {
            function(theta)
                .Call(glmerLaplace, pp$ptr, resp$ptr,
                      theta, u0, beta0, verbose, FALSE, tolPwrss)
        } else {
            function(theta) {
                pp$u0 <- u0
                pp$beta0 <- beta0
                pp$theta <- theta
                pwrssUpdate(pp, resp, verbose, tol=tolPwrss)
                resp$Laplace(pp$ldL2(), pp$ldRX2(), pp$sqrL(0))
            }
        }
        environment(devfun) <- rho
        if (devFunOnly) return(devfun)
	control$iprint <- min(verbose, 3L)
	opt <- bobyqa(pp$theta, devfun, lower=reTrms$lower, control=control)
	if (nAGQ == 1L) {
            rho$u0 <- pp$u0
            rho$dpars <- seq_along(pp$theta)
	    if(!is.numeric(control$rhobeg)) control$rhobeg <- 0.0002
	    if(!is.numeric(control$rhoend)) control$rhoend <- 2e-7
            rho$control <- control
            devfunb <- if (compDev) {
                function(pars)
                    .Call(glmerLaplace, pp$ptr, resp$ptr, pars[dpars],
                          u0, pars[-dpars], verbose, TRUE, tolPwrss)
            } else {
                function(pars) {
                    pp$u0 <- u0
                    pp$theta <- pars[dpars]
                    pp$beta0 <- pars[-dpars]
                    pwrssUpdate(pp, resp, verbose, uOnly=TRUE, tol=tolPwrss)
                    resp$Laplace(pp$ldL2(), pp$ldRX2(), pp$sqrL(0))
                }
            }
            environment(devfunb) <- rho
            opt <- bobyqa(c(pp$theta, pp$beta0), devfunb,
                          lower=c(reTrms$lower, rep.int(-Inf, length(pp$beta0))),
                          control=control)
        }
    }
    LL <- pp$Lambdat
    LL@x <- pp$theta[pp$Lind]
    stopifnot(validObject(LL))
    pp$Lambdat <- LL
    sqrLenU <- pp$sqrL(0.)
    wrss <- resp$wrss()
    pwrss <- wrss + sqrLenU
    n <- nrow(fr)

    dims <- c(N=n, n=n, nmp=n-p, nth=length(pp$theta), p=p, q=nrow(reTrms$Zt),
	      nAGQ=nAGQ, useSc=0L, reTrms=length(reTrms$cnms),
	      spFe=0L, REML=0L, GLMM=1L, NLMM=0L)
    cmp <- c(ldL2=pp$ldL2(), ldRX2=pp$ldRX2(), wrss=wrss,
             ussq=sqrLenU, pwrss=pwrss,
	     drsum=resp$resDev(), dev=opt$fval, REML=NA,
	     sigmaML=NA, sigmaREML=NA)

    new("glmerMod", call=mc, frame=fr, flist=reTrms$flist, cnms=reTrms$cnms,
	theta=pp$theta, beta=pp$beta0, u=pp$u0, lower=reTrms$lower,
	devcomp=list(cmp=cmp, dims=dims), pp=pp, resp=resp)
}## {glmer}

##' Create an lmerResp, glmerResp or (later) nlmerResp instance
##'
##' @title Create a [ng]lmerResp instance
##' @param fr a model frame
##' @param family the optional glm family (glmRespMod only)
##' @param nlenv the nonlinear model evaluation environment (nlsRespMod only)
##' @param nlmod the nonlinear model function (nlsRespMod only)
##' @return a lmerResp (or glmerResp) instance
mkRespMod2 <- function(fr, family = NULL, nlenv = NULL, nlmod = NULL) {
    n <- nrow(fr)
    y <- model.response(fr)
    if(length(dim(y)) == 1) {
	## avoid problems with 1D arrays, but keep names
	nm <- rownames(y)
	dim(y) <- NULL
	if(!is.null(nm)) names(y) <- nm
    }
    yy <- y
    if (is.factor(yy)) yy <- as.numeric(yy != levels(yy)[1])
    if (is.matrix(yy) && ncol(yy) == 2L) yy <- yy[,1]/(yy[,1] + yy[,2])
    ans <- new("lmerResp", y = yy)
    if (!is.null(offset <- model.offset(fr))) {
	if (length(offset) == 1L) offset <- rep.int(offset, n)
	stopifnot(length(offset) == n)
	ans$offset <- unname(offset)
    }
    if (!is.null(weights <- model.weights(fr))) {
	stopifnot(length(weights) == n, all(weights >= 0))
	ans$weights <- unname(weights)
    }
    if (is.null(family)) return(ans)
    rho <- new.env()
    rho$etastart <- model.extract(fr, "etastart")
    rho$mustart <- model.extract(fr, "mustart")
    rho$weights <- ans$weights
    rho$offset <- ans$offset
    rho$nobs <- n
    if (is.null(rho$y)) rho$y <- ans$y
    eval(family$initialize, rho)
    family$initialize <- NULL	    # remove clutter from str output

    ans <- new("glmerResp", family, rho$y)
    ans$weights <- rho$weights
    ans$offset <- rho$offset
    ans$updateMu(family$linkfun(unname(rho$mustart)))
    ans
}

###' Create Z, Lambda, Lind, etc.
##'
##' @param bars a list of parsed random-effects terms
##' @param fr a model frame in which to evaluate these terms
##' @param s Number of parameters in the nonlinear mean function (nlmer only)
##'
##' @return a list of Zt, Lambdat, Lind, theta, lower, flist and cnms
mkReTrms <- function(bars, fr, s = 1L) {
    if (!length(bars))
	stop("No random effects terms specified in formula")
    stopifnot(is.list(bars), all(sapply(bars, is.language)),
	      inherits(fr, "data.frame"))
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))

    ## auxiliary {named, for easier inspection}:
    mkBlist <- function(x) {
	ff <- eval(substitute(factor(fac), list(fac = x[[3]])), fr)
	if (all(is.na(ff)))
	    stop("Invalid grouping factor specification, ",
		 deparse(x[[3]]))
	nl <- length(levels(ff))
	mm <- model.matrix(eval(substitute( ~ foo,
					   list(foo = x[[2]]))), fr)
	nc <- ncol(mm)
	nseq <- seq_len(nc)
	sm <- as(ff, "sparseMatrix")
	if (nc	> 1)
	    sm <- do.call(rBind, lapply(nseq, function(i) sm))
	sm@x[] <- t(mm[])
	## When nc > 1 switch the order of the rows of sm
	## so the random effects for the same level of the
	## grouping factor are adjacent.
	if (nc > 1)
	    sm <- sm[as.vector(matrix(seq_len(nc * nl),
				      nc = nl, byrow = TRUE)),]
	list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
    }
    blist <- lapply(bars, mkBlist)
    nl <- unlist(lapply(blist, "[[", "nl")) # no. of levels per term

    ## order terms stably by decreasing number of levels in the factor
    if (any(diff(nl)) > 0) {
	ord <- rev(order(nl))
	blist <- blist[ord]
	nl <- nl[ord]
    }
    Zt <- do.call(rBind, lapply(blist, "[[", "sm"))
    q <- nrow(Zt)

    ## Create and install Lambdat, Lind, etc.  This must be done after
    ## any potential reordering of the terms.
    cnms <- lapply(blist, "[[", "cnms")
    nc <- sapply(cnms, length)		# no. of columns per term
    nth <- as.integer((nc * (nc+1))/2)	# no. of parameters per term
    nb <- nc * nl			# no. of random effects per term
    stopifnot(sum(nb) == q)
    boff <- cumsum(c(0L, nb))		# offsets into b
    thoff <- cumsum(c(0L, nth))		# offsets into theta
### FIXME: should this be done with cBind and avoid the transpose
### operator?  In other words should Lambdat be generated directly
### instead of generating Lambda first then transposing?
    Lambdat <-
	t(do.call(sparseMatrix,
		  do.call(rBind,
			  lapply(seq_along(blist), function(i)
			     {
				 mm <- matrix(seq_len(nb[i]), nc = nc[i],
					      byrow = TRUE)
				 dd <- diag(nc[i])
				 ltri <- lower.tri(dd, diag = TRUE)
				 ii <- row(dd)[ltri]
				 jj <- col(dd)[ltri]
				 dd[cbind(ii, jj)] <- seq_along(ii)
				 data.frame(i = as.vector(mm[, ii]) + boff[i],
					    j = as.vector(mm[, jj]) + boff[i],
					    x = as.double(rep.int(seq_along(ii),
					    rep.int(nl[i], length(ii))) +
					    thoff[i]))
			     }))))
    thet <- numeric(sum(nth))
    ll <- list(Zt = Zt, theta = thet, Lind = as.integer(Lambdat@x))
    ## lower bounds on theta elements are 0 if on diagonal, else -Inf
    ll$lower <- -Inf * (thet + 1)
    ll$lower[unique(diag(Lambdat))] <- 0
    ll$theta[] <- is.finite(ll$lower) # initial values of theta are 0 off-diagonal, 1 on
    Lambdat@x[] <- ll$theta[ll$Lind]  # initialize elements of Lambdat
    ll$Lambdat <- Lambdat
					# massage the factor list
    fl <- lapply(blist, "[[", "ff")
					# check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
	fl <- fl[match(ufn, fnms)]
	asgn <- match(fnms, ufn)
    } else asgn <- seq_along(fl)
    names(fl) <- ufn
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn
    ll$flist <- fl
    ll$cnms <- cnms
    ll
} ## {mkReTrms2}

##' Determine a step factor that will reduce the pwrss
##'
##' The penalized, weighted residual sum of squares (pwrss) is the sum
##' of the weighted residual sum of squares from the resp module and
##' the squared length of u from the predictor module.  The predictor module
##' contains a base value and an increment for the coefficients.
##' @title Determine a step factor
##' @param pp predictor module
##' @param resp response module
##' @param verbose logical value determining verbose output
##' @return NULL if successful
stepFac <- function(pp, resp, verbose, maxSteps = 10) {
    stopifnot(is.numeric(maxSteps), maxSteps >= 2)
    pwrss0 <- resp$wrss() + pp$sqrL(0)
    for (fac in 2^(-(0:maxSteps))) {
	wrss <- resp$updateMu(pp$linPred(fac))
	pwrss1 <- wrss + pp$sqrL(fac)
	if (verbose > 3L)
	    cat(sprintf("pwrss0=%10g, diff=%10g, fac=%6.4f\n",
			pwrss0, pwrss0 - pwrss1, fac))
	if (pwrss1 <= pwrss0) {
	    pp$installPars(fac)
	    return(NULL)
	}
    }
    stop("step factor reduced below 0.001 without reducing pwrss")
}

pwrssUpdate <- function(pp, resp, verbose, uOnly=FALSE, tol, maxSteps = 10) {
    stopifnot(is.numeric(tol), tol > 0)
    repeat {
	resp$updateMu(pp$linPred(0))
	resp$updateWts()
        pp$updateXwts(resp$sqrtXwt())
        pp$updateDecomp()
        pp$updateRes(resp$wtres())
        if (uOnly) pp$solveU() else pp$solve()
	if ((pp$CcNumer())/(resp$wrss() + pp$sqrL(0)) < tol)
	    break
	stepFac(pp, resp, verbose, maxSteps=maxSteps)
    }
}

## setMethod("show", signature("lmerResp"), function(object)
##       {
##           with(object,
##                print(head(cbind(weights, offset, mu, y, sqrtrwt, wtres, sqrtXwt))))
##       })

## setMethod("show", signature("glmerResp"), function(object)
##       {
##           with(object,
##                print(head(cbind(weights, offset, eta, mu, y, muEta,
##                                 variance, sqrtrwt, wtres, sqrtXwt,
##                                 sqrtWrkWt, wrkResids, wrkResp))))
##       })

## setMethod("show", signature("Rcpp_reModule"), function(object)
##	 {
##	     with(object, print(head(cbind(u0, incr, u))))
##	 })


## setMethod("show", signature("Rcpp_deFeMod"), function(object)
##	 {
##	     with(object, print(cbind(coef0, incr, coef, Vtr)))
##	 })

## create a deviance evaluation function that uses the sigma parameters
## df2 <- function(dd) {
##     stopifnot(is.function(dd),
## 	      length(formals(dd)) == 1L,
## 	      is((rem <- (rho <- environment(dd))$rem), "Rcpp_reModule"),
## 	      is((fem <- rho$fem), "Rcpp_deFeMod"),
## 	      is((resp <- rho$resp), "Rcpp_lmerResp"),
## 	      all((lower <- rem$lower) == 0))
##     Lind <- rem$Lind
##     n <- length(resp$y)
##     function(pars) {
## 	sigma <- pars[1]
## 	sigsq <- sigma * sigma
## 	sigmas <- pars[-1]
## 	theta <- sigmas/sigma
## 	rem$theta <- theta
## 	resp$updateMu(numeric(n))
## 	solveBetaU(rem, fem, resp$sqrtXwt, resp$wtres)
## 	resp$updateMu(rem$linPred1(1) + fem$linPred1(1))
## 	n * log(2*pi*sigsq) + (resp$wrss + rem$sqrLenU)/sigsq + rem$ldL2
##     }
## }

simulate.merMod <- function(object, nsim = 1, seed = NULL, use.u = FALSE, ...)
{
    stopifnot((nsim <- as.integer(nsim[1])) > 0,
	      is(object, "merMod"),
	      ## i.e. not yet for glmer etc:
	      is(object@resp, "lmerResp"))
    if(!is.null(seed)) set.seed(seed)
    if(!exists(".Random.seed", envir = .GlobalEnv))
	runif(1) # initialize the RNG if necessary

    n <- nrow(X <- object@X)
    ## result will be matrix  n x nsim :
    as.vector(X %*% object@beta) +  # fixed-effect contribution
	sigma(object) * (## random-effects contribution:
			 if(use.u) {
			     object@ u
			 } else {
			     U <- object@Z %*% object@Lambda
			     q <- ncol(U)
			     as(U %*% matrix(rnorm(q * nsim), nc = nsim), "matrix")
			 }
			 ## residual contribution:
			 + matrix(rnorm(n * nsim), nc = nsim))
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

##' <description>
##'  Perform model-based (Semi-)parametric bootstrap for mixed models;
##'  The working name for bootMer() was simulestimate(), but this is serious..
##' <details>
##'  ...
##' @title Model-based (Semi-)Parametric Bootstrap for Mixed Models
##' @param x fitted *lmer() model
##' @param FUN a function(x) computing the *statistic* of interest,
##' which must be a numeric vector, possibly named.
##' @param nsim number of simulations, positive integer; the bootstrap
##' 'B' (or 'R').
##' @param seed optional argument to \code{\link{set.seed}}.
##' @param use.u
##' @param verbose logical indicating if progress should print output
##' @param control
##' @return an object of S3 class "boot" {compatible with boot package}
##' @author Martin Maechler
bootMer <- function(x, FUN, nsim = 1, seed = NULL, use.u = FALSE,
		    verbose = FALSE, control = list())
{
    stopifnot((nsim <- as.integer(nsim[1])) > 0,
	      is(x, "merMod"), is(x@resp, "lmerResp"))
    FUN <- match.fun(FUN)
    if(!is.null(seed)) set.seed(seed)
    else if(!exists(".Random.seed", envir = .GlobalEnv))
	runif(1) # initialize the RNG if necessary

    mc <- match.call()
    t0 <- FUN(x)
    if (!is.numeric(t0))
	stop("bootMer currently only handles functions that return numeric vectors")

    ## simplistic approach {original for old-lme4 by DB was much smarter}

    n <- nrow(X <- getME(x, "X"))
    if(use.u) {
	u <- getME(x,"u")
    } else {
	U <- getME(x,"Z") %*% getME(x, "Lambda")
	q <- ncol(U)
    }
    ##    Zt <- x@Z

    beta <- getME(x, "beta")
    X.beta <- as.vector(X %*% beta) # fixed-effect contribution
    sigm.x <- sigma(x)

    ## Here, and below ("optimize"/"bobyqa") using the "logic" of lmer() itself:
## lmer..Update <- if(is(x, "lmerSp")) lmerSpUpdate else lmerDeUpdate
#    devfun <- mkdevfun(x)
##    oneD <- length(x@re@theta) < 2
    theta0 <- getME(x,"theta")
    ## just for the "boot" result -- TODOmaybe drop
    mle <- list(beta = beta, theta = theta0, sigma = sigm.x)

    t.star <- matrix(t0, nr = length(t0), nc = nsim)
    resp <- x@resp
    for(i in 1:nsim) {
	y <- {
	    X.beta + sigm.x *
		((if(use.u) u else as.vector(U %*% rnorm(q))) + rnorm(n))
	    ##	    random effects  contribution	    +	  Error
	}
	x @ resp <- new(Class=class(resp), REML=resp$REML, y=y, offset=resp$offset, weights=resp$weights)

	## if (oneD) { # use optimize
	##     d0 <- devfun(0)
	##     opt <- optimize(devfun, c(0, 10))
	##     ##		       -------- <<< arbitrary
	##     ## FIXME ?! if optimal theta > 0, optimize will *not* warn!
	##     if (d0 <= opt$objective) ## prefer theta == 0 when close
	##	   devfun(0) # -> theta	 := 0  and update the rest
	## } else {
	opt <- bobyqa(theta0, mkdevfun(resp, x@pp), x@lower, control = control)
##	  xx <- updateMod(x, opt$par, opt$fval)
	    ## FIXME: also here, prefer \hat\sigma^2 == 0 (exactly)
##	  }
	foo <- tryCatch(FUN(xx), error = function(e)e)
	if(verbose) { cat(sprintf("%5d :",i)); str(foo) }
	t.star[,i] <- if (inherits(foo, "error")) NA else foo
    }
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

##' Fit a nonlinear mixed-effects model
##'
##' @param formula a nonlinear mixed model formula (see detailed documentation)
##' @param data an optional data frame containing the variables named in
##'    \code{formula}.	By default the variables are taken from the
##'    environment from which \code{nlmer} is called.
##' @param family
##' @param start starting estimates for the nonlinear model
##'    parameters, as a named numeric vector
##' @param verbose integer scalar passed to nlminb.  If negative then
##'    diagnostic output from the PIRLS (penalized iteratively
##'    reweighted least squares) step is also provided.
##' @param nAGQ number of adaptive Gauss-Hermite quadrature points to use
##' @param doFit logical scalar.  If FALSE the optimization
##'    environment is returned. Otherwise the parameters are estimated
##'    and an object of S4 class "mer" is returned.
##' @param subset further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param weights  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param na.action  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param mustart
##' @param etastart
##' @param sparseX
##' @param contrasts  further model specifications as in
##'    \code{\link[stats]{lm}}; see there for details.
##' @param control a list of control parameters passed to bobyqa.
##' @param ...

##' @return an object of S4 class "merMod"
nlmer <- function(formula, data, family = gaussian, start = NULL,
		  verbose = 0L, nAGQ = 1L, doFit = TRUE,
		  subset, weights, na.action, mustart, etastart,
		  sparseX = FALSE, contrasts = NULL, control = list(), ...)
{
    if (!missing(family)) stop("code not yet written")
    mf <- mc <- match.call()
    m <- match(c("data", "subset", "weights", "na.action",
		 "offset", "etastart", "mustart"),
	       names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
					# check the formula
    formula <- as.formula(formula)
    if (length(formula) < 3) stop("formula must be a 3-part formula")
    nlform <- as.formula(formula[[2]])
					# Really do need to check twice
    if (length(nlform) < 3) stop("formula must be a 3-part formula")
    nlmod <- as.call(nlform[[3]])
					# check for parameter names in start
    if (is.numeric(start)) start <- list(nlpars = start)
    stopifnot((s <- length(pnames <- names(start$nlpars))) > 0,
	      is.numeric(start$nlpars))
    if (!all(pnames %in% (anms <- all.vars(nlmod))))
	stop("not all parameter names are used in the nonlinear model expression")
    fr.form <- nlform
## FIXME: This should be changed to use subbars and subnms.
    fr.form[[3]] <-
	parse(text = paste(setdiff(all.vars(formula), pnames),
			 collapse = ' + '))[[1]]
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())

    ## First create nlenv.  For this the nlpar columns are numeric
    for (nm in pnames) fr[[nm]] <- start$nlpars[[nm]]
    nlenv <- new.env()	# inherit from this environment (or environment(formula)?)
    lapply(all.vars(nlmod),
	   function(nm) assign(nm, fr[[nm]], envir = nlenv))

    ## Second, extend the frame and convert the nlpar columns to indicators
    n <- nrow(fr)
    frE <- do.call(rbind, lapply(1:s, function(i) fr)) # rbind s copies of the frame
    for (nm in pnames) # convert these variables in fr to indicators
	frE[[nm]] <- as.numeric(rep(nm == pnames, each = n))
					# random-effects module
    reTrms <- mkReTrms(findbars(formula[[3]]), frE, s = s)
    dcmp <- updateDcmp(reTrms, .dcmp())
    dcmp$dims["nAGQ"] <- as.integer(nAGQ)[1]

    fe.form <- nlform
    fe.form[[3]] <- formula[[3]]
    feMod <- mkFeModule(fe.form, frE, contrasts, reTrms, sparseX = sparseX)
					# should this check be in mkFeModule?
    p <- length(feMod@coef)
    if ((qrX <- qr(feMod@X))$rank < p)
	stop(gettextf("rank of X = %d < ncol(X) = %d", qrX$rank, p))
    feMod@coef <- qr.coef(qrX, unlist(lapply(pnames, get, envir = nlenv)))
    respMod <- mkRespMod2(fr, nlenv = nlenv, nlmod = nlmod)
    ans <- new("merMod",
	       call = mc,
	       devcomp = updateDcmp(respMod, updateDcmp(feMod, dcmp)),
	       frame = fr, re = reTrms, fe = feMod, resp = respMod)
    if (!doFit) return(ans)
    PIRLSest(ans, verbose, control, nAGQ)
}


## Methods for the merMod class
fixef.merMod <- function(object, ...)
    structure(object@beta, names = dimnames(object@pp$X)[[2]])

formula.merMod <- function(x, ...) formula(x@call, ...)

##' Extract the random effects.
##'
##' Extract the conditional modes, which for a linear mixed model are
##' also the conditional means, of the random effects, given the
##' observed responses.	 These also depend on the model parameters.
##'
##' @param object an object that inherits from the \code{\linkS4class{mer}} class
##' @param postVar logical scalar - should the posterior variance be returned
##' @param drop logical scalar - drop dimensions of single extent
##' @param whichel - vector of names of factors for which to return results

##' @return a named list of arrays or vectors, aligned to the factor list

ranef.merMod <- function(object, postVar = FALSE, drop = FALSE,
			 whichel = names(ans), ...)
{
    ans <- as.vector(crossprod(object@pp$Lambdat, object@u))
    if (!is.null(object@flist)) {
	## evaluate the list of matrices
	levs <- lapply(fl <- object@flist, levels)
	asgn <- attr(fl, "assign")
	cnms <- object@cnms
	nc <- sapply(cnms, length)
	nb <- nc * (nl <- unlist(lapply(levs, length))[asgn])
	nbseq <- rep.int(seq_along(nb), nb)
	ml <- split(ans, nbseq)
	for (i in seq_along(ml))
	    ml[[i]] <- matrix(ml[[i]], nc = nc[i], byrow = TRUE,
			      dimnames = list(NULL, cnms[[i]]))
	## create a list of data frames corresponding to factors
	ans <- lapply(seq_along(fl),
		      function(i)
		      data.frame(do.call(cbind, ml[asgn == i]),
				 row.names = levs[[i]],
				 check.names = FALSE))
	names(ans) <- names(fl)
					# process whichel
	stopifnot(is(whichel, "character"))
	whchL <- names(ans) %in% whichel
	ans <- ans[whchL]

	if (postVar) {
	    vv <- .Call(reTrmsCondVar, re, sigma(object))
	    for (i in seq_along(ans))
		attr(ans[[i]], "postVar") <- vv[[i]]
	}
	if (drop)
	    ans <- lapply(ans, function(el)
		      {
			  if (ncol(el) > 1) return(el)
			  pv <- drop(attr(el, "postVar"))
			  el <- drop(as.matrix(el))
			  if (!is.null(pv))
			      attr(el, "postVar") <- pv
			  el
		      })
	class(ans) <- "ranef.mer"
    }
    ans
}

setMethod("sigma", "merMod", function(object, ...)
      {
	  dc <- object@devcomp
	  dd <- dc$dims
	  if(dd[["useSc"]])
	      dc$cmp[[ifelse(dd[["REML"]], "sigmaREML", "sigmaML")]] else 1.
      })

model.matrix.merMod <- function(object, ...) object@X

terms.merMod <- function(x, ...) attr(x@frame, "terms")

model.frame.merMod <- function(formula, ...) formula@frame

deviance.merMod <- function(object, REML = NULL, ...) {
    if (!missing(REML)) stop("REML argument not supported")
    object@devcomp$cmp[["dev"]]
}

logLik.merMod <- function(object, REML = NULL, ...)
{
    if (!missing(REML)) stop("REML argument not supported")
    dc <- object@devcomp
    dims <- dc$dims
    val <- - dc$cmp[["dev"]]/2
    attr(val, "nall") <- attr(val, "nobs") <- nrow(object@frame)
    attr(val, "df") <- length(object@beta) + length(object@theta) + dims[["useSc"]]
    class(val) <- "logLik"
    val
}

update.merMod <- function(object, formula., ..., evaluate = TRUE)
{
    if (is.null(call <- getCall(object)))
	stop("object should contain a 'call' component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
	call$formula <- update.formula(formula(object), formula.)
    if (length(extras) > 0) {
	existing <- !is.na(match(names(extras), names(call)))
	for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	if (any(!existing)) {
	    call <- c(as.list(call), extras[!existing])
	    call <- as.call(call)
	}
    }
    if (evaluate)
	eval(call, parent.frame())
    else call
}

#update.merMod <- updateMer

## This is modeled a bit after	print.summary.lm :
## Prints *both*  'mer' and 'merenv' - as it uses summary(x) mainly
printMerenv <- function(x, digits = max(3, getOption("digits") - 3),
			correlation = NULL, symbolic.cor = FALSE,
			signif.stars = getOption("show.signif.stars"), ...)
{
    so <- summary(x)
    cat(sprintf("%s ['%s']\n",so$methTitle, class(x)))
    if (!is.null(f <- so$family)) {
	cat(" Family:", f)
        if (!(is.null(ll <- so$link))) cat(" (", ll, ")")
        cat("\n")
    }
    if (!is.null(cc <- so$call$formula))
	cat("Formula:", deparse(cc),"\n")
    if (!is.null(cc <- so$call$data))
	cat("   Data:", deparse(cc), "\n")
    if (!is.null(cc <- so$call$subset))
	cat(" Subset:", deparse(asOneSidedFormula(cc)[[2]]),"\n")
    cat("\n")
    tab <- so$AICtab
    if (length(tab) == 1 && names(tab) == "REML")
	cat("REML criterion at convergence:", round(tab, 4), "\n")
    else print(round(so$AICtab, 4))
    cat("\nRandom effects:\n")
    print(formatVC(so$varcor, digits = digits, useScale = so$useScale),
	  quote = FALSE, digits = digits, ...)

    ngrps <- so$ngrps
    cat(sprintf("Number of obs: %d, groups: ", so$devcomp$dims[["n"]]))
    cat(paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "))
    cat("\n")
    p <- nrow(so$coefficients)
    if (p > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(so$coefficients, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
	if(!is.logical(correlation)) { # default
	    correlation <- p <= 20
	    if(!correlation) {
		nam <- deparse(substitute(x)) # << TODO: improve if this is called from show()
		cat(sprintf(paste("\nCorrelation matrix not shown by default, as p = %d > 20.",
				  "Use print(%s, correlation=TRUE)  or",
				  "    vcov(%s)	 if you need it\n", sep="\n"),
			    p, nam, nam))
	    }
	}
	if(correlation) {
	    if(is.null(VC <- so$vcov)) VC <- vcov(x)
	    corF <- VC@factors$correlation
	    if (is.null(corF)) {
		cat("\nCorrelation of Fixed Effets is not available\n")
	    }
	    else {
		p <- ncol(corF)
		if (p > 1) {
		    rn <- rownames(so$coefficients)
		    rns <- abbreviate(rn, minlen=11)
		    cat("\nCorrelation of Fixed Effects:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			corf <- as(corF, "matrix")
			dimnames(corf) <- list(rns,
					       abbreviate(rn, minlength=1, strict=TRUE))
			print(symnum(corf))
		    }
		    else {
			corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				       nc = p,
				       dimnames = list(rns, abbreviate(rn, minlen=6)))
			corf[!lower.tri(corf)] <- ""
			print(corf[-1, -p, drop=FALSE], quote = FALSE)
		    }
		}
	    }
	}
    }
    invisible(x)
}## printMerenv()

print.merMod <- printMerenv
setMethod("show",  "merMod", function(object) printMerenv(object))
fitted.merMod <- function(object, ...) {object <- object@resp; NextMethod()}
residuals.merMod <- function(object, type = c("deviance", "pearson",
				     "working", "response", "partial"), ...)
{
    object <- object@resp
    NextMethod()
}

print.summary.mer <- printMerenv
#setMethod("show",  "summary.mer", function(object) printMerenv(object))


## coef() method for all kinds of "mer", "*merMod", ... objects
## ------  should work with fixef() + ranef()  alone
coefMer <- function(object, ...)
{
    if (length(list(...)))
	warning(paste('arguments named "',
		      paste(names(list(...)), collapse = ", "),
		      '" ignored', sep = ''))
    fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
    ref <- ranef(object)
    val <- lapply(ref, function(x)
		  fef[rep.int(1L, nrow(x)),,drop = FALSE])
    for (i in seq(a = val)) {
	refi <- ref[[i]]
	row.names(val[[i]]) <- row.names(refi)
	nmsi <- colnames(refi)
	if (!all(nmsi %in% names(fef)))
	    stop("unable to align random and fixed effects")
	for (nm in nmsi) val[[i]][[nm]] <- val[[i]][[nm]] + refi[,nm]
    }
    class(val) <- "coef.mer"
    val
} ##  {coefMer}

coef.merMod <- coefMer

## FIXME: Do we really need a separate devcomp extractor?  I suppose it can't hurt.
devcomp <- function(x, ...) UseMethod("devcomp")
devcomp.merMod <- function(x, ...) x@devcomp
##setMethod("devcomp", "merMod", function(x, ...) x@devcomp)

setMethod("getL", "merMod", function(x) {
    .Deprecated("getME(., \"L\")")
    getME(x, "L")
})

##' "Generalized Extractor" -- providing compatibility between lme4* versions -- lme4Eigen ---
##' @param object [ng]lmer() fit
##' @param name character string
##' @return the corresponding "part" of the [gn]?lmer()-Fit
##' @note
##'
getME <- function(object,
		  name = c("X", "Z","Zt", "u",
		  "Gp",
		  "L", "Lambda", "Lambdat",
		  "RX", "RZX",
                  "beta", "theta",
		  "REML", "n_rtrms", "is_REML"))
{
    if(missing(name)) stop("'name' must not be missing")
    stopifnot(length(name <- as.character(name)) == 1,
	      is(object, "merMod"))
    name <- match.arg(name)
    rsp <- object@resp
    PR <- object@pp
    switch(name,
	   "X" = PR$X, ## ok ? - check -- use model.matrix() method instead?
	   "Z" = t(PR$Zt),
	   "Zt"= PR$Zt,
           "u" = PR$ u0 + PR$ delu(),
           ## "Gp"
	   "L"= PR$ L(),
	   "Lambda"= t(PR$ Lambdat),
	   "Lambdat"= PR$ Lambdat,
	   "RX" = PR $ RX,
	   "RZX" = PR $ RZX,

           "beta" = object@beta,
           "theta"= object@theta, ## *OR*  PR $ theta  --- which one ??

	   "REML" = rsp $ REML,
	   "is_REML" = as.logical(rsp $ REML),# correct ??

	   "n_rtrms" =, ## FIXME length(PR$flist), ##  = #{random-effect terms in the formula}
	   "Gp"= ,  # FIXME
	   "..foo.." =# placeholder!
	   stop(gettextf("'%s' is not implemented yet",
			 sprintf("getME(*, \"%s\")", name))),
	   ## otherwise
	   stop(sprintf("Mixed-Effects extraction of '%s' is not available for class \"%s\"",
			name, class(object))))
}## {getME}


isREML <- function(x) UseMethod("isREML")
isREML.merMod <- function(x) as.logical(x@devcomp$dims["REML"])

refitML <- function(x) UseMethod("refitML")
## Can't change x in the call to this function.
refitML.merMod <- function (x) {
    if (!isREML(x)) return(x)
    stopifnot(is(rr <- x@resp, "lmerResp"))
    resp <- new(class(rr), y=rr$y)
    resp$offset <- rr$offset
    resp$weights <- rr$weights
    resp$REML <- 0L
    xpp <- x@pp
    pp <- new(class(xpp), X=xpp$X, Zt=xpp$Zt, Lambdat=xpp$Lambdat,
	      Lind=xpp$Lind, theta=xpp$theta)
    opt <- bobyqa(x@theta, mkdevfun(pp, resp), x@lower)
    n <- length(rr$y)
    p <- ncol(pp$X)
    dims <- c(N=n, n=n, nmp=n-p, nth=length(pp$theta), p=p, q=nrow(pp$Zt),
	      nAGQ=NA_integer_, useSc=1L, reTrms=length(x@cnms),
	      spFe=0L, REML=0L, GLMM=0L, NLMM=0L)
    wrss <- resp$wrss()
    ussq <- pp$sqrL(1)
    pwrss <- wrss + ussq
    cmp <- c(ldL2=pp$ldL2(), ldRX2=pp$ldRX2(), wrss=wrss, ussq=ussq,
	     pwrss=pwrss, drsum=NA, dev=opt$fval, REML=NA,
	     sigmaML=sqrt(pwrss/n), sigmaREML=sqrt(pwrss/(n-p)))

### FIXME: Should modify the call slot to set REML=FALSE.  It is
### tricky to do so without causing the call to be evaluated
    new("lmerMod", call=x@call, frame=x@frame, flist=x@flist,
	cnms=x@cnms, theta=pp$theta, beta=pp$delb(), u=pp$delu(),
	lower=x@lower, devcomp=list(cmp=cmp, dims=dims), pp=pp, resp=resp)
}

getCall.merMod <- function(x, ...) x@call


##' <description>
##'
##' <details>
##' @title vcov(): Extract conditional covariance matrix of fixed effects
##' @param sigma = sigma(object)
##' @param correlation
##' @param ...
mkVcov <- function(sigma, unsc, nmsX, correlation = TRUE, ...) {
    V <- sigma^2 * unsc
    if(is.null(rr <- tryCatch(as(V, "dpoMatrix"),
			      error = function(e) NULL)))
	stop("Computed variance-covariance matrix is not positive definite")
    dimnames(rr) <- list(nmsX, nmsX)
    if(correlation)
	rr@factors$correlation <- as(rr, "corMatrix")
    rr
}

vcov.merMod <- function(object, correlation = TRUE, sigm = sigma(object), ...)
    mkVcov(sigm, unsc = object@pp$unsc(), nmsX = colnames(object@pp$X),
	   correlation=correlation, ...)

vcov.summary.mer <- function(object, correlation = TRUE, ...)
{
    if(!is.null(object$vcov))
	vcov
    else if(!is.null(PP <- object$pp))
	mkVcov(object$sigma, RX = , nmsX = colnames(FE@X),
	       correlation=correlation, ...)
    else stop("Both 'vcov' and 'fe' components are missing.  You need\n",
	      "at least one TRUE in summary(..,	 varcov = *, keep.X = *)")
}

mkVarCorr <- function(sc, cnms, nc, theta, nms) {
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
    ans <- lapply(ncseq, function(i)
	      {
		  ## Li := \Lambda_i, the i-th block diagonal of \Lambda(\theta)
		  Li <- diag(nrow = nc[i])
		  Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
		  rownames(Li) <- cnms[[i]]
		  ## val := \Sigma_i = \sigma^2 \Lambda_i \Lambda_i', the
		  val <- tcrossprod(sc * Li) # variance-covariance
		  stddev <- sqrt(diag(val))
		  correl <- t(val / stddev)/stddev
		  diag(correl) <- 1
		  attr(val, "stddev") <- stddev
		  attr(val, "correlation") <- correl
		  val
	      })
    if(is.character(nms)) names(ans) <- nms
    attr(ans, "sc") <- sc
    ans
}

VarCorr.merMod <- function(x, sigma, rdig)# <- 3 args from nlme
{
    if (is.null(cnms <- x@cnms))
	stop("VarCorr methods require reTrms, not just reModule")
    if(missing(sigma)) # "bug": fails via default 'sigma=sigma(x)'
	sigma <- lme4Eigen::sigma(x)
    nc <- sapply(cnms, length)	  # no. of columns per term
    mkVarCorr(sigma, cnms=cnms, nc=nc, theta = x@theta,
	      nms = {fl <- x@flist; names(fl)[attr(fl, "assign")]})
}

##' <description>
##' Compute standard errors of fixed effects from an lmer()
##'
##' <details>
##' @title
##' @param object "lmerenv" object,
##' @param RX the Cholesky factor (CHMfactor) L of ...
##' @return numeric vector of length length(fixef(.))
##' @author Doug Bates & Martin Maechler
##' currently *not* exported on purpose
unscaledVar <- function(object, RX = object@fe@fac)
{
    if (is(RX, "Cholesky")) return(diag(chol2inv(RX)))
    stopifnot(is(RX, "CHMfactor"))
    p <- ncol(RX)
    if (p < 1) return(numeric(0))
    p1 <- p - 1L
    ei <- as(c(1, rep.int(0, p1)), "sparseMatrix")
    DI <- function(i) {
	ei@i <- i
	sum(solve(RX, ei, sys = "L")@x^2)
    }
    as.vector(solve(RX, unlist(lapply(0:p1, DI)),
		    system = "Pt"))
}

formatVC <- function(varc, digits = max(3, getOption("digits") - 2),
		     useScale)
### "format()" the 'VarCorr' matrix of the random effects -- for show()ing
{
    sc <- unname(attr(varc, "sc"))
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- c(lapply(varc, attr, "stddev"), if(useScale) list(Residual = sc))
    reLens <- unlist(c(lapply(reStdDev, length)))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- c(unlist(lapply(varc, colnames)), if(useScale) "")
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <-
	    do.call("rBind",
		    lapply(recorr,
			   function(x) {
			       x <- as(x, "matrix")
			       cc <- format(round(x, 3), nsmall = 3)
			       cc[!lower.tri(cc)] <- ""
			       nr <- dim(cc)[1]
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }))[, -maxlen, drop = FALSE]
	if (nrow(corr) < nrow(reMat))
	    corr <- rbind(corr, matrix("", nr = nrow(reMat) - nrow(corr), nc = ncol(corr)))
	colnames(corr) <- rep.int("", ncol(corr))
	colnames(corr)[1] <- "Corr"
	cbind(reMat, corr)
    } else reMat
}

##' <description>
##'
##' @title Summary Method for *mer() fits, i.e., "merMod" objects
##' @param object
##' @param varcov logical indicating if vcov(.) should be computed and stored.
##' @param keep.X logical indicating if the 'fe' component of object should be stored;
##'   the default is true when 'varcov' is false, as we then need fe for vcov()
##' @param ...
##' @return S3 class "summary.mer", basically a list .....
summary.merMod <- function(object, ...)
{
    resp <- object@resp
    devC <- object@devcomp
    dd <- devC$dims
    cmp <- devC$cmp
    useSc <- as.logical(dd["useSc"])
    sig <- sigma(object)
    REML <- isREML(object)

    link <- fam <- NULL
    if(is(resp, "glmerResp")) {
        fam <- resp$fam()
        link <- resp$link()
    }
    coefs <- cbind("Estimate" = fixef(object),
		   "Std. Error" = sig * sqrt(diag(object@pp$unsc())))
    if (nrow(coefs) > 0) {
	coefs <- cbind(coefs, coefs[,1]/coefs[,2], deparse.level=0)
	colnames(coefs)[3] <- paste(if(useSc) "t" else "z", "value")
    }
    mName <- paste(switch(1L + dd["GLMM"] * 2L + dd["NLMM"],
			  "Linear", "Nonlinear",
			  "Generalized linear", "Generalized nonlinear"),
		   "mixed model fit by",
		   ifelse(REML, "REML", "maximum likelihood"))
    llik <- logLik(object)   # returns NA for a REML fit - maybe change?
    AICstats <- {
	if (REML) cmp["REML"] # do *not* show likelihood stats here
	else {
	    c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
	      deviance = deviance(object))
	}
    }
    ## FIXME: You can't count on object@re@flist,
    ##	      nor compute VarCorr() unless is(re, "reTrms"):
    varcor <- VarCorr(object)
					# use S3 class for now
    structure(list(methTitle=mName, devcomp=devC,
                   isLmer=is(resp, "lmerResp"), useScale=useSc,
		   logLik=llik, family=fam, link=link,
		   ngrps=sapply(object@flist, function(x) length(levels(x))),
		   coefficients=coefs, sigma=sig,
		   vcov=vcov(object, correlation=TRUE, sigm=sig),
		   varcor=varcor, # and use formatVC(.) for printing.
		   AICtab=AICstats, call=object@call
		   ), class = "summary.mer")
}

summary.summary.mer <- function(object, varcov = FALSE, ...)
{
    if(varcov && is.null(object$vcov))
	object$vcov <- vcov(object, correlation=TRUE, sigm = object$sigma)
    object
}



### Plots for the ranef.mer class ----------------------------------------

dotplot.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
	if (is.null(se)) return(list())
	x <- as.numeric(x)
	hw <- 1.96 * as.numeric(se[subscripts])
	list(xlim = range(x - hw, x + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16,
			 horizontal = TRUE, col = dot.symbol$col,
			 lty = dot.line$lty, lwd = dot.line$lwd,
			 col.line = dot.line$col, levels.fos = unique(y),
			 groups = NULL, ...)
    {
	x <- as.numeric(x)
	y <- as.numeric(y)
	dot.line <- trellis.par.get("dot.line")
	dot.symbol <- trellis.par.get("dot.symbol")
	sup.symbol <- trellis.par.get("superpose.symbol")
	panel.abline(h = levels.fos, col = col.line, lty = lty, lwd = lwd)
	panel.abline(v = 0, col = col.line, lty = lty, lwd = lwd)
	if (!is.null(se)) {
	    se <- as.numeric(se[subscripts])
	    panel.segments( x - 1.96 * se, y, x + 1.96 * se, y, col = 'black')
	}
	panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(x, ...) {
	ss <- stack(x)
	ss$ind <- factor(as.character(ss$ind), levels = colnames(x))
	ss$.nn <- rep.int(reorder(factor(rownames(x)), x[[1]]), ncol(x))
	se <- NULL
	if (!is.null(pv <- attr(x, "postVar")))
	    se <- unlist(lapply(1:(dim(pv)[1]), function(i) sqrt(pv[i, i, ])))
	dotplot(.nn ~ values | ind, ss, se = se,
		prepanel = prepanel.ci, panel = panel.ci,
		xlab = NULL, ...)
    }
    lapply(x, f, ...)
}

plot.ranef.mer <- function(x, y, ...)
{
    lapply(x, function(x) {
	cn <- lapply(colnames(x), as.name)
	switch(min(ncol(x), 3),
	       qqmath(eval(substitute(~ x, list(x = cn[[1]]))), x, ...),
	       xyplot(eval(substitute(y ~ x,
				      list(y = cn[[1]],
					   x = cn[[2]]))), x, ...),
	       splom(~ x, ...))
    })
}

qqmath.ranef.mer <- function(x, data, ...)
{
    prepanel.ci <- function(x, y, se, subscripts, ...) {
	x <- as.numeric(x)
	se <- as.numeric(se[subscripts])
	hw <- 1.96 * se
	list(xlim = range(x - hw, x + hw, finite = TRUE))
    }
    panel.ci <- function(x, y, se, subscripts, pch = 16, ...)  {
	panel.grid(h = -1,v = -1)
	panel.abline(v = 0)
	x <- as.numeric(x)
	y <- as.numeric(y)
	se <- as.numeric(se[subscripts])
	panel.segments(x - 1.96 * se, y, x + 1.96 * se, y, col = 'black')
	panel.xyplot(x, y, pch = pch, ...)
    }
    f <- function(x) {
	if (!is.null(pv <- attr(x, "postVar"))) {
	    cols <- 1:(dim(pv)[1])
	    se <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
	    nr <- nrow(x)
	    nc <- ncol(x)
	    ord <- unlist(lapply(x, order)) +
		rep((0:(nc - 1)) * nr, each = nr)
	    rr <- 1:nr
	    ind <- gl(ncol(x), nrow(x), labels = names(x))
	    xyplot(rep(qnorm((rr - 0.5)/nr), ncol(x)) ~ unlist(x)[ord] | ind[ord],
		   se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
		   scales = list(x = list(relation = "free")),
		   ylab = "Standard normal quantiles",
		   xlab = NULL, ...)
	} else {
	    qqmath(~values|ind, stack(x),
		   scales = list(y = list(relation = "free")),
		   xlab = "Standard normal quantiles",
		   ylab = NULL, ...)
	}
    }
    lapply(x, f)
}

plot.coef.mer <- function(x, y, ...)
{
    ## remove non-varying columns from frames
    reduced <- lapply(x, function(el)
		      el[, !sapply(el, function(cc) all(cc == cc[1]))])
    plot.ranef.mer(reduced, ...)
}

dotplot.coef.mer <- function(x, data, ...) {
    mc <- match.call()
    mc[[1]] <- as.name("dotplot.ranef.mer")
    eval(mc)
}


## A more effective method for merMod objects is defined in lmer.R ??
## (Not sure this comment is still operative, as in Ron Ziegler's
## famous reply about "That statement is no longer operative.")

##' <description>
##'
##' <details>
##' @title anova() for both  "lmerenv" and "lmer" fitted models
##' @param object an "lmerenv", "lmer" or "lmerMod" - fitted model
##' @param ...	further such objects
##' @return an "anova" data frame; the traditional (S3) result of anova()
anovaLmer <- function(object, ...) {
    mCall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
    modp <- {
#	as.logical(.sapply(dots, is, "lmerenv")) |
	as.logical(.sapply(dots, is, "merMod")) |
#	as.logical(.sapply(dots, is, "lmer")) |
	as.logical(.sapply(dots, is, "lm")) }
    if (any(modp)) {			# multiple models - form table
	opts <- dots[!modp]
	mods <- c(list(object), dots[modp])
	## model names
	mNms <- .sapply(as.list(mCall)[c(FALSE, TRUE, modp)], deparse)
	names(mods) <- sub("@env$", '', mNms) # <- hack
	mods <- lapply(mods, refitML)

	devs <- sapply(mods, deviance)
	llks <- lapply(mods, logLik)
	ii <- order(Df <- .sapply(llks, attr, "df"))
	mods <- mods[ii]
	llks <- llks[ii]
	Df   <- Df  [ii]
	calls <- lapply(mods, getCall)
	data <- lapply(calls, "[[", "data")
	if (any(data != data[[1]]))
	    stop("all models must be fit to the same data object")
	header <- paste("Data:", data[[1]])
	subset <- lapply(calls, "[[", "subset")
	if (any(subset != subset[[1]]))
	    stop("all models must use the same subset")
	if (!is.null(subset[[1]]))
	    header <-
		c(header, paste("Subset", deparse(subset[[1]]),
				sep = ": "))
	llk <- unlist(llks)
	chisq <- 2 * pmax(0, c(NA, diff(llk)))
	dfChisq <- c(NA, diff(Df))
	val <- data.frame(Df = Df,
			  AIC = .sapply(llks, AIC),
			  BIC = .sapply(llks, BIC),
			  logLik = llk,
			  deviance = -2*llk,
			  Chisq = chisq,
			  "Chi Df" = dfChisq,
			  "Pr(>Chisq)" = pchisq(chisq, dfChisq, lower = FALSE),
			  row.names = names(mods), check.names = FALSE)
	class(val) <- c("anova", class(val))
	attr(val, "heading") <-
	    c(header, "Models:",
	      paste(rep(names(mods), times = unlist(lapply(lapply(lapply(calls,
				     "[[", "formula"), deparse), length))),
		    unlist(lapply(lapply(calls, "[[", "formula"), deparse)),
		    sep = ": "))
	return(val)
    }
    else { ## ------ single model ---------------------
	dc <- devcomp(object)
	p <- dc$dims["p"]
	asgn <- object@ X @ assign
	stopifnot(length(asgn) == p,
		  is(object@fe, "deFeMod"), # haven't worked out sparse version
		  is(object@re, "reTrms"))  # things are really weird with no re terms
	ss <- ((object@pp@RX %*% object@beta)@x)^2
	names(ss) <- colnames(object@X)
	terms <- terms(object)
	nmeffects <- setdiff(attr(terms, "term.labels"), names(object@re@flist))
	if ("(Intercept)" %in% names(ss))
	    nmeffects <- c("(Intercept)", nmeffects)
	ss <- unlist(lapply(split(ss, asgn), sum))
	stopifnot(length(ss) == length(nmeffects))
	df <- unlist(lapply(split(asgn,	 asgn), length))
	## dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
	ms <- ss/df
	f <- ms/(sigma(object)^2)
	## P <- pf(f, df, dfr, lower.tail = FALSE)
	## table <- data.frame(df, ss, ms, dfr, f, P)
	table <- data.frame(df, ss, ms, f)
	dimnames(table) <-
	    list(nmeffects,
		 ## c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		 c("Df", "Sum Sq", "Mean Sq", "F value"))
	if ("(Intercept)" %in% nmeffects)
	    table <- table[-match("(Intercept)", nmeffects), ]
	attr(table, "heading") <- "Analysis of Variance Table"
	class(table) <- c("anova", "data.frame")
	table
    }
}## {anovaLmer}

anova.merMod <- anovaLmer

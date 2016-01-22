## NB:  doc in ../man/*.Rd  ***not*** auto generated
## FIXME: need to document S3 methods better (can we pull from r-forge version?)

##' Fit a linear mixed model (LMM)
lmer <- function(formula, data=NULL, REML = TRUE,
                 control = lmerControl(), start = NULL,
                 verbose = 0L, subset, weights, na.action, offset,
                 contrasts = NULL, devFunOnly=FALSE,
                 ...)
{
    mc <- mcout <- match.call()
    missCtrl <- missing(control)
    ## see functions in modular.R for the body ...
    if (!missCtrl && !inherits(control, "lmerControl")) {
        if(!is.list(control)) stop("'control' is not a list; use lmerControl()")
        ## back-compatibility kluge
	warning("passing control as list is deprecated: please use lmerControl() instead",
		immediate.=TRUE)
        control <- do.call(lmerControl, control)
    }
    if (!is.null(list(...)[["family"]])) {
       warning("calling lmer with 'family' is deprecated; please use glmer() instead")
       mc[[1]] <- quote(lme4::glmer)
       if(missCtrl) mc$control <- glmerControl()
       return(eval(mc, parent.frame(1L)))
    }
    mc$control <- control ## update for  back-compatibility kluge

    ## https://github.com/lme4/lme4/issues/50
    ## parse data and formula
    mc[[1]] <- quote(lme4::lFormula)
    lmod <- eval(mc, parent.frame(1L))
    mcout$formula <- lmod$formula
    lmod$formula <- NULL

    ## create deviance function for covariance parameters (theta)
    devfun <- do.call(mkLmerDevfun,
		      c(lmod,
			list(start=start, verbose=verbose, control=control)))
    if (devFunOnly) return(devfun)
    ## optimize deviance function over covariance parameters
    opt <- if (control$optimizer=="none")
        list(par=NA,fval=NA,conv=1000,message="no optimization")
    else {
        optimizeLmer(devfun, optimizer = control$optimizer,
                     restart_edge = control$restart_edge,
                     boundary.tol = control$boundary.tol,
                     control = control$optCtrl,
                     verbose=verbose,
                     start=start,
                     calc.derivs=control$calc.derivs,
                     use.last.params=control$use.last.params)
    }
    cc <- checkConv(attr(opt,"derivs"), opt$par,
                    ctrl = control$checkConv,
                    lbound=environment(devfun)$lower)
    mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
             mc = mcout, lme4conv=cc) ## prepare output
}## { lmer }


##' Fit a generalized linear mixed model (GLMM)
glmer <- function(formula, data=NULL, family = gaussian,
                  control = glmerControl(), start = NULL, verbose = 0L, nAGQ = 1L,
                  subset, weights, na.action, offset,
                  contrasts = NULL, mustart, etastart, devFunOnly = FALSE, ...)
{
    if (!inherits(control, "glmerControl")) {
	if(!is.list(control)) stop("'control' is not a list; use glmerControl()")
	## back-compatibility kluge
	msg <- "Use control=glmerControl(..) instead of passing a list"
	if(length(cl <- class(control))) msg <- paste(msg, "of class", dQuote(cl[1]))
	warning(msg, immediate.=TRUE)
	control <- do.call(glmerControl, control)
    }
    mc <- mcout <- match.call()

    ## family-checking code duplicated here and in glFormula (for now) since
    ## we really need to redirect at this point; eventually deprecate formally
    ## and clean up
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if( is.function(family)) family <- family()
    if (isTRUE(all.equal(family, gaussian()))) {
        ## redirect to lmer (with warning)
        warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;",
                " please call lmer() directly")
        mc[[1]] <- quote(lme4::lmer)
        mc["family"] <- NULL            # to avoid an infinite loop
        return(eval(mc, parent.frame()))
    }

    ## see https://github.com/lme4/lme4/issues/50
    ## parse the formula and data
    mc[[1]] <- quote(lme4::glFormula)
    glmod <- eval(mc, parent.frame(1L))
    mcout$formula <- glmod$formula
    glmod$formula <- NULL

    ## create deviance function for covariance parameters (theta)

    nAGQinit <- if(control$nAGQ0initStep) 0L else 1L
    devfun <- do.call(mkGlmerDevfun, c(glmod, list(verbose = verbose,
                                                   control = control,
                                                   nAGQ = nAGQinit)))
    if (nAGQ==0 && devFunOnly) return(devfun)
    ## optimize deviance function over covariance parameters

    ## FIXME: perhaps should be in glFormula instead??
    if (is.list(start)) {
        start.bad <- setdiff(names(start),c("theta","fixef"))
        if (length(start.bad)>0) {
            stop(sprintf("bad name(s) for start vector (%s); should be %s and/or %s",
                         paste(start.bad,collapse=", "),
                         shQuote("theta"),
                         shQuote("fixef")),call.=FALSE)
        }
        if (!is.null(start$fixef) && nAGQ==0)
            stop("should not specify both start$fixef and nAGQ==0")
    }

    ## FIX ME: allow calc.derivs, use.last.params etc. if nAGQ=0
    if(control$nAGQ0initStep) {
        opt <- optimizeGlmer(devfun,
                             optimizer = control$optimizer[[1]],
                             ## DON'T try fancy edge tricks unless nAGQ=0 explicitly set
                             restart_edge=if (nAGQ==0) control$restart_edge else FALSE,
                             boundary.tol=if (nAGQ==0) control$boundary.tol else 0,
                             control = control$optCtrl,
                             start=start,
                             nAGQ = 0,
                             verbose=verbose,
                             calc.derivs=FALSE)
    }

    if(nAGQ > 0L) {

        start <- updateStart(start,theta=opt$par)

        ## update deviance function to include fixed effects as inputs
        devfun <- updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)

        if (devFunOnly) return(devfun)
        ## reoptimize deviance function over covariance parameters and fixed effects
        opt <- optimizeGlmer(devfun,
                             optimizer = control$optimizer[[2]],
                             restart_edge=control$restart_edge,
                             boundary.tol=control$boundary.tol,
                             control = control$optCtrl,
                             start=start,
                             nAGQ=nAGQ,
                             verbose = verbose,
                             stage=2,
                             calc.derivs=control$calc.derivs,
                             use.last.params=control$use.last.params)
    }
    cc <- if (!control$calc.derivs) NULL else {
        if (verbose > 10) cat("checking convergence\n")
        checkConv(attr(opt,"derivs"),opt$par,
                  ctrl = control$checkConv,
                  lbound=environment(devfun)$lower)
    }

    ## prepare output
    mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr,
             mc = mcout, lme4conv=cc)

}## {glmer}

##' Fit a nonlinear mixed-effects model
nlmer <- function(formula, data=NULL, control = nlmerControl(), start = NULL, verbose = 0L,
                  nAGQ = 1L, subset, weights, na.action, offset,
                  contrasts = NULL, devFunOnly = FALSE, ...)
{

    vals <- nlformula(mc <- match.call())
    p <- ncol(X <- vals$X)
    if ((rankX <- rankMatrix(X)) < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))

    rho <- list2env(list(verbose=verbose,
                         tolPwrss=0.001, # this is reset to the tolPwrss argument's value later
                         resp=vals$resp,
                         lower=vals$reTrms$lower),
                    parent=parent.frame())
    rho$pp <- do.call(merPredD$new,
                      c(vals$reTrms[c("Zt","theta","Lambdat","Lind")],
                        list(X=X, n=length(vals$respMod$mu), Xwts=vals$respMod$sqrtXwt,
                             beta0=qr.coef(qr(X), unlist(lapply(vals$pnames, get,
                             envir = rho$resp$nlenv))))))
    rho$u0 <- rho$pp$u0
    rho$beta0 <- rho$pp$beta0
    ## deviance as a function of theta only :
    devfun <- mkdevfun(rho, 0L, verbose=verbose, control=control)
    if (devFunOnly && !nAGQ) return(devfun)
    devfun(rho$pp$theta) # initial coarse evaluation to get u0 and beta0
    rho$u0 <- rho$pp$u0
    rho$beta0 <- rho$pp$beta0
    rho$tolPwrss <- control$tolPwrss # Reset control parameter (the initial optimization is coarse)

    opt <- optwrap(control$optimizer[[1]], devfun, rho$pp$theta, lower=rho$lower,
                   control=control$optCtrl, adj=FALSE)
    rho$control <- attr(opt,"control")

    if (nAGQ > 0L) {
        rho$lower <- c(rho$lower, rep.int(-Inf, length(rho$beta0)))
        rho$u0    <- rho$pp$u0
        rho$beta0 <- rho$pp$beta0
        rho$dpars <- seq_along(rho$pp$theta)
        if (nAGQ > 1L) {
            if (length(vals$reTrms$flist) != 1L || length(vals$reTrms$cnms[[1]]) != 1L)
                stop("nAGQ > 1 is only available for models with a single, scalar random-effects term")
            rho$fac <- vals$reTrms$flist[[1]]
        }
        devfun <- mkdevfun(rho, nAGQ, verbose=verbose, control=control)
        if (devFunOnly) return(devfun)

        opt <- optwrap(control$optimizer[[2]], devfun,
                       par = c(rho$pp$theta, rho$beta0),
                       lower = rho$lower, control = control$optCtrl,
                       adj = TRUE, verbose=verbose)

    }
    mkMerMod(environment(devfun), opt, vals$reTrms, fr = vals$frame, mc = mc)
}## {nlmer}

## R 3.1.0 devel [2013-08-05]: This does not help yet
if(getRversion() >= "3.1.0") utils::suppressForeignCheck("nlmerAGQ")
if(getRversion() < "3.1.0") dontCheck <- identity

##' Create a deviance evaluation function from a predictor and a response module
mkdevfun <- function(rho, nAGQ=1L, maxit=100L, verbose=0, control=list()) {
    ## FIXME: should nAGQ be automatically embedded in rho?
    stopifnot(is.environment(rho), is(rho$resp, "lmResp"))

    ## silence R CMD check warnings *locally* in this function
    ## (clearly preferred to using globalVariables() !]
    fac <- pp <- resp <- lp0 <- compDev <- dpars <- baseOffset <- tolPwrss <-
	pwrssUpdate <- ## <-- even though it's a function below
	GQmat <- nlmerAGQ <- NULL

    ## The deviance function (to be returned):
    ff <-
    if (is(rho$resp, "lmerResp")) {
	rho$lmer_Deviance <- lmer_Deviance
	function(theta) .Call(lmer_Deviance, pp$ptr(), resp$ptr(), as.double(theta))
    } else if (is(rho$resp, "glmResp")) {
        ## control values will override rho values *if present*
        if (!is.null(tp <- control$tolPwrss)) rho$tolPwrss <- tp
        if (!is.null(cd <- control$ compDev)) rho$compDev <- cd
	if (nAGQ == 0L)
	    function(theta) {
		resp$updateMu(lp0)
		pp$setTheta(theta)
		p <- pwrssUpdate(pp, resp, tol=tolPwrss, GQmat=GHrule(0L),
                                 compDev=compDev, maxit=maxit, verbose=verbose)
                resp$updateWts()
                p
            }
	else  ## nAGQ > 0
	    function(pars) {
                ## pp$setDelu(rep(0, length(pp$delu)))
                resp$setOffset(baseOffset)
		resp$updateMu(lp0)
		pp$setTheta(as.double(pars[dpars])) # theta is first part of pars
                spars <- as.numeric(pars[-dpars])
                offset <- if (length(spars)==0) baseOffset else baseOffset + pp$X %*% spars
		resp$setOffset(offset)
		p <- pwrssUpdate(pp, resp, tol=tolPwrss, GQmat=GQmat,
                                 compDev=compDev, grpFac=fac, maxit=maxit, verbose=verbose)
                resp$updateWts()
                p
	    }
    } else if (is(rho$resp, "nlsResp")) {
	if (nAGQ < 2L) {
	    rho$nlmerLaplace <- nlmerLaplace
            rho$tolPwrss <- control$tolPwrss
	    switch(nAGQ + 1L,
                   function(theta)
                       .Call(nlmerLaplace, pp$ptr(), resp$ptr(), as.double(theta),
                             as.double(u0), beta0, verbose, FALSE, tolPwrss),
                   function(pars)
                       .Call(nlmerLaplace, pp$ptr(), resp$ptr(), pars[dpars],
                             u0,  pars[-dpars], verbose, TRUE, tolPwrss))
	} else {
            stop("AGQ>1 not yet implemented for nlmer models")
	    rho$nlmerAGQ <- nlmerAGQ
	    rho$GQmat	 <- GHrule(nAGQ)
	    ## function(pars) {
            ## .Call(nlmerAGQ, ## <- dontCheck(nlmerAGQ)  should work according to docs but does not
            ## pp$ptr(), resp$ptr(), fac, GQmat, pars[dpars],
            ## u0, pars[-dpars], tolPwrss)
            ##}
	}
    }
    else stop("code not yet written")
    environment(ff) <- rho
    ff
}

## Determine a step factor that will reduce the pwrss
##
## The penalized, weighted residual sum of squares (pwrss) is the sum
## of the weighted residual sum of squares from the resp module and
## the squared length of u from the predictor module.  The predictor module
## contains a base value and an increment for the coefficients.
## @title Determine a step factor
## @param pp predictor module
## @param resp response module
## @param verbose logical value determining verbose output
## @return NULL if successful
## @note Typically all this is done in the C++ code.
##     The R code is for debugging and comparisons of
##     results.
## stepFac <- function(pp, resp, verbose, maxSteps = 10) {
##     stopifnot(is.numeric(maxSteps), maxSteps >= 2)
##     pwrss0 <- resp$wrss() + pp$sqrL(0)
##     for (fac in 2^(-(0:maxSteps))) {
## 	wrss <- resp$updateMu(pp$linPred(fac))
## 	pwrss1 <- wrss + pp$sqrL(fac)
## 	if (verbose > 3L)
## 	    cat(sprintf("pwrss0=%10g, diff=%10g, fac=%6.4f\n",
## 			pwrss0, pwrss0 - pwrss1, fac))
## 	if (pwrss1 <= pwrss0) {
## 	    pp$installPars(fac)
## 	    return(NULL)
## 	}
##     }
##     stop("step factor reduced below ",signif(2^(-maxSteps),2)," without reducing pwrss")
## }

RglmerWrkIter <- function(pp, resp, uOnly=FALSE) {
    pp$updateXwts(resp$sqrtWrkWt())
    pp$updateDecomp()
    pp$updateRes(resp$wtWrkResp())
    if (uOnly) pp$solveU() else pp$solve()
    resp$updateMu(pp$linPred(1))	# full increment
    resp$resDev() + pp$sqrL(1)
}

##' @param pp pred module
##' @param resp resp module
##' @param tol numeric tolerance
##' @param GQmat matrix of Gauss-Hermite quad info
##' @param compDev compute in C++ (as opposed to doing as much as possible in R)
##' @param grpFac grouping factor (normally found in environment ...)
##' @param verbose verbosity, of course
glmerPwrssUpdate <- function(pp, resp, tol, GQmat, compDev=TRUE, grpFac=NULL, maxit = 70L, verbose=0) {
    nAGQ <- nrow(GQmat)
    if (compDev) {
        if (nAGQ < 2L)
            return(.Call(glmerLaplace, pp$ptr(), resp$ptr(),
                         nAGQ, tol, as.integer(maxit),
                         verbose))
        return(.Call(glmerAGQ, pp$ptr(), resp$ptr(),
                     tol, as.integer(maxit),
                     GQmat, grpFac, verbose))
    }
    oldpdev <- .Machine$double.xmax
    uOnly   <- nAGQ == 0L
    i <- 0
    repeat {
        ## oldu <- pp$delu
        ## olddelb <- pp$delb
        pdev <- RglmerWrkIter(pp, resp, uOnly=uOnly)
        if (verbose > 2) cat(i,": ",pdev,"\n",sep="")
        ## check convergence first so small increases don't trigger errors
        if (is.na(pdev)) stop("encountered NA in PWRSS update")
        if (abs((oldpdev - pdev) / pdev) < tol)
            break
        ## if (pdev > oldpdev) {
        ##     ## try step-halving
        ##     ## browser()
        ##     k <- 0
        ##     while (k < 10 && pdev > oldpdev) {
        ##         pp$setDelu((oldu + pp$delu)/2.)
        ##         if (!uOnly) pp$setDelb((olddelb + pp$delb)/2.)
        ##         pdev <- RglmerWrkIter(pp, resp, uOnly=uOnly)
        ##         k <- k+1
        ##     }
        ## }
        if (pdev > oldpdev) stop("PIRLS update failed")
        oldpdev <- pdev
        i <- i+1
    }
    resp$Laplace(pp$ldL2(), 0., pp$sqrL(1))  ## FIXME: should 0. be pp$ldRX2 ?
}

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

## bootMer() ---> now in ./bootMer.R


## Methods for the merMod class

## Anova for merMod objects
##
## @title anova() for merMod objects
## @param a merMod object
## @param ...	further such objects
## @param refit should objects be refitted with ML (if applicable)
## @return an "anova" data frame; the traditional (S3) result of anova()
anovaLmer <- function(object, ..., refit = TRUE, model.names=NULL) {
    mCall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
    modp <- (as.logical(vapply(dots, is, NA, "merMod")) |
             as.logical(vapply(dots, is, NA, "lm")))
    if (any(modp)) {			# multiple models - form table
	## opts <- dots[!modp]
	mods <- c(list(object), dots[modp])
        nobs.vec <- vapply(mods, nobs, 1L)
        if (var(nobs.vec) > 0)
            stop("models were not all fitted to the same size of dataset")

	## model names
        if (is.null(mNms <- model.names))
            mNms <- vapply(as.list(mCall)[c(FALSE, TRUE, modp)], safeDeparse, "")

        ## HACK to try to identify model names in situations such as
        ## 'do.call(anova,list(model1,model2))' where the model names
        ## are lost in the call stack ... this doesn't quite work but might
        ## be useful for future attempts?
        ## maxdepth <- -2
        ## depth <- -1
        ## while (depth >= maxdepth &
        ##        all(grepl("S4 object of class structure",mNms))) {
        ##     xCall <- match.call(call=sys.call(depth))
        ##     mNms <- .sapply(as.list(xCall)[c(FALSE, TRUE, modp)], deparse)
        ##     depth <- depth-1
        ## }
        ## if (depth < maxdepth) {
        if (any(duplicated(mNms))) {
            warning("failed to find unique model names, assigning generic names")
            mNms <- paste0("MODEL",seq_along(mNms))
        }
        if (length(mNms) != length(mods))
            stop("model names vector and model list have different lengths")
	names(mods) <- sub("@env$", '', mNms) # <- hack
	models.reml <- vapply(mods, function(x) is(x,"merMod") && isREML(x), NA)
        models.GHQ <- vapply(mods, function(x) is(x,"glmerMod") && getME(x,"devcomp")$dims["nAGQ"]>1 , NA)
        if (any(models.GHQ) && any(vapply(mods, function(x) is(x,"glm"), NA)))
            stop("GLMMs with nAGQ>1 have log-likelihoods incommensurate with glm() objects")
	if (refit) {
	    ## message only if at least one models is REML:
	    if (any(models.reml)) message("refitting model(s) with ML (instead of REML)")
	    mods[models.reml] <- lapply(mods[models.reml], refitML)
	} else { ## check that models are consistent (all REML or all ML)
	    if(any(models.reml) && any(!models.reml))
		warning("some models fit with REML = TRUE, some not")
	}
	## devs <- sapply(mods, deviance)
	llks <- lapply(mods, logLik)
        ## Order models by increasing degrees of freedom:
	ii <- order(Df <- vapply(llks, attr, FUN.VALUE=numeric(1), "df"))
	mods <- mods[ii]
	llks <- llks[ii]
	Df   <- Df  [ii]
	calls <- lapply(mods, getCall)
	data <- lapply(calls, `[[`, "data")
	if(!all(vapply(data, identical, NA, data[[1]])))
	    stop("all models must be fit to the same data object")
	header <- paste("Data:", abbrDeparse(data[[1]]))
	subset <- lapply(calls, `[[`, "subset")
	if(!all(vapply(subset, identical, NA, subset[[1]])))
	    stop("all models must use the same subset")
	if (!is.null(subset[[1]]))
	    header <- c(header, paste("Subset:", abbrDeparse(subset[[1]])))
	llk <- unlist(llks)
	chisq <- 2 * pmax(0, c(NA, diff(llk)))
	dfChisq <- c(NA, diff(Df))
	val <- data.frame(Df = Df,
			  AIC = .sapply(llks, AIC), # FIXME? vapply()
			  BIC = .sapply(llks, BIC), #  "       "
                          logLik = llk,
			  deviance = -2*llk,
			  Chisq = chisq,
			  "Chi Df" = dfChisq,
			  "Pr(>Chisq)" = pchisq(chisq, dfChisq, lower.tail = FALSE),
			  row.names = names(mods), check.names = FALSE)
	class(val) <- c("anova", class(val))
	forms <- lapply(lapply(calls, `[[`, "formula"), deparse)
	structure(val,
		  heading = c(header, "Models:",
			      paste(rep(names(mods), times = lengths(forms)),
				    unlist(forms), sep = ": ")))
    }
    else { ## ------ single model ---------------------
	dc <- getME(object, "devcomp")
        X <- getME(object, "X")
	stopifnot(length(asgn <- attr(X, "assign")) == dc$dims[["p"]])
	ss <- as.vector(object@pp$RX() %*% object@beta)^2
	names(ss) <- colnames(X)
	terms <- terms(object)
        nmeffects <- attr(terms, "term.labels")[unique(asgn)]
	if ("(Intercept)" %in% names(ss))
	    nmeffects <- c("(Intercept)", nmeffects)
	ss <- unlist(lapply(split(ss, asgn), sum))
	stopifnot(length(ss) == length(nmeffects))
	df <- lengths(split(asgn, asgn))
	## dfr <- unlist(lapply(split(dfr, asgn), function(x) x[1]))
	ms <- ss/df
	f <- ms/(sigma(object)^2)
        ## No longer provide p-values, but still the F statistic (may not be F distributed):
        ##
	## P <- pf(f, df, dfr, lower.tail = FALSE)
	## table <- data.frame(df, ss, ms, dfr, f, P)
	table <- data.frame(df, ss, ms, f)
	dimnames(table) <-
	    list(nmeffects,
		 ## c("Df", "Sum Sq", "Mean Sq", "Denom", "F value", "Pr(>F)"))
		 c("Df", "Sum Sq", "Mean Sq", "F value"))
	if ("(Intercept)" %in% nmeffects)
	    table <- table[-match("(Intercept)", nmeffects), ]
        structure(table, heading = "Analysis of Variance Table",
                  class = c("anova", "data.frame"))
    }
}## {anovaLmer}

##' @importFrom stats anova
##' @S3method anova merMod
anova.merMod <- anovaLmer

##' @S3method as.function merMod
as.function.merMod <- function(x, ...) {
    rho <- list2env(list(resp  = x@resp$copy(),
                         pp    = x@pp$copy(),
                         beta0 = x@beta,
                         u0   =  x@u),
                    parent=as.environment("package:lme4"))
    ## FIXME: extract verbose [, maxit] and control
    mkdevfun(rho, getME(x, "devcomp")$dims[["nAGQ"]])
}

## coef() method for all kinds of "mer", "*merMod", ... objects
## ------  should work with fixef() + ranef()  alone
coefMer <- function(object, ...)
{
    if (length(list(...)))
	warning('arguments named "', paste(names(list(...)), collapse = ", "),
                '" ignored')
    fef <- data.frame(rbind(fixef(object)), check.names = FALSE)
    ref <- ranef(object)
    ## check for variables in RE but missing from FE, fill in zeros in FE accordingly
    refnames <- unlist(lapply(ref,colnames))
    nmiss <- length(missnames <- setdiff(refnames,names(fef)))
    if (nmiss > 0) {
        fillvars <- setNames(data.frame(rbind(rep(0,nmiss))),missnames)
        fef <- cbind(fillvars,fef)
    }
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

##' @importFrom stats coef
##' @S3method coef merMod
coef.merMod <- coefMer

## FIXME: should these values (i.e. ML criterion for REML models
##  and vice versa) be computed and stored in the object in the first place?
##' @importFrom stats deviance
##' @S3method deviance merMod
deviance.merMod <- function(object, REML = NULL, ...) {
                            ## type = c("conditional", "unconditional", "penalized"),
                            ## relative = TRUE, ...) {
    if (isGLMM(object)) {
        return(sum(residuals(object,type="deviance")^2))
        ## ------------------------------------------------------------
        ## proposed change to deviance function for GLMMs
        ## ------------------------------------------------------------
        ## @param type Type of deviance (can be unconditional,
        ## penalized, conditional)
        ## @param relative Should deviance be shifted relative to a
        ## saturated model? (only available with type == penalized or
        ## conditional)
        ## ------------------------------------------------------------
        ## ans <- switch(type[1],
        ##               unconditional = {
        ##                   if (relative) {
        ##                       stop("unconditional and relative deviance is undefined")
        ##                   }
        ##                   c(-2 * logLik(object))
        ##               },
        ##               penalized = {
        ##                   sqrL <- object@pp$sqrL(1)
        ##                   if (relative) {
        ##                       object@resp$resDev() + sqrL
        ##                   } else {
        ##                       useSc <- unname(getME(gm1, "devcomp")$dims["useSc"])
        ##                       qLog2Pi <- unname(getME(object, "q")) * log(2 * pi)
        ##                       object@resp$aic() - (2 * useSc) + sqrL + qLog2Pi
        ##                   }
        ##               },
        ##               conditional = {
        ##                   if (relative) {
        ##                       object@resp$resDev()
        ##                   } else {
        ##                       useSc <- unname(getME(gm1, "devcomp")$dims["useSc"])
        ##                       object@resp$aic() - (2 * useSc)
        ##                   }
        ##               })
        ## return(ans)
    }
    if (isREML(object) && is.null(REML)) {
        warning("deviance() is deprecated for REML fits; use REMLcrit for the REML criterion or deviance(.,REML=FALSE) for deviance calculated at the REML fit")
        return(devCrit(object, REML=TRUE))
    }
    devCrit(object, REML=FALSE)
}

REMLcrit <- function(object) {
    devCrit(object, REML=TRUE)
}

## original deviance.merMod -- now wrapped by REMLcrit
## REML=NULL:
##    if REML fit return REML criterion
##    if ML fit, return deviance
## REML=TRUE:
##    if not LMM, stop.
##    if ML fit, compute and return REML criterion
##    if REML fit, return REML criterion
## REML=FALSE:
##    if ML fit, return deviance
##    if REML fit, compute and return deviance
devCrit <- function(object, REML = NULL) {
    ## cf. (1) lmerResp::Laplace in respModule.cpp
    ##     (2) section 5.6 of lMMwR, listing lines 34-42
    if (isTRUE(REML) && !isLMM(object))
        stop("can't compute REML deviance for a non-LMM")
    cmp <- object@devcomp$cmp
    if (is.null(REML) || is.na(REML[1]))
        REML <- isREML(object)
    if (REML) {
        if (isREML(object)) {
            cmp[["REML"]]
        } else {
            ## adjust ML results to REML
	    lnum <- log(2*pi*cmp[["pwrss"]])
	    n <- object@devcomp$dims[["n"]]
	    nmp <- n - length(object@beta)
            ldW <- sum(log(weights(object, method = "prior")))
            - ldW + cmp[["ldL2"]] + cmp[["ldRX2"]] + nmp*(1 + lnum - log(nmp))
        }
    } else {
        if (!isREML(object)) {
            cmp[["dev"]]
        } else {
            ## adjust REML results to ML
            n <- object@devcomp$dims[["n"]]
            lnum <- log(2*pi*cmp[["pwrss"]])
            ldW <- sum(log(weights(object, method = "prior")))
            - ldW + cmp[["ldL2"]] + n*(1 + lnum - log(n))
        }
    }
}

## copied from stats:::safe_pchisq
safe_pchisq <- function (q, df, ...) {
    df[df <= 0] <- NA
    pchisq(q = q, df = df, ...)
}

##' @importFrom stats drop1
##' @S3method drop1 merMod
drop1.merMod <- function(object, scope, scale = 0, test = c("none", "Chisq", "user"),
                         k = 2, trace = FALSE,
                         sumFun=NULL, ...) {
    evalhack <- "formulaenv"
    test <- match.arg(test)
    if ((test=="user" && is.null(sumFun)) ||
        ((test!="user" && !is.null(sumFun))))
        stop(sQuote("sumFun"),' must be specified if (and only if) test=="user"')
    tl <- attr(terms(object), "term.labels")
    if(missing(scope)) scope <- drop.scope(object)
    else {
	if(!is.character(scope)) {
	    scope <- attr(terms(getFixedFormula(update.formula(object, scope))),
                                "term.labels")
        }
	if(!all(match(scope, tl, 0L) > 0L))
	    stop("scope is not a subset of term labels")
    }
    ns <- length(scope)
    if (is.null(sumFun)) {
        sumFun <- function(x,scale,k,...)
            setNames(extractAIC(x,scale,k,...),c("df","AIC"))
    }
    ss <- sumFun(object, scale=scale, k=k, ...)
    ans <- matrix(nrow = ns + 1L, ncol = length(ss),
                  dimnames =  list(c("<none>", scope), names(ss)))
    ans[1, ] <- ss
    n0 <- nobs(object, use.fallback = TRUE)
    env <- environment(formula(object)) # perhaps here is where trouble begins??
    for(i in seq_along(scope)) {  ## was seq(ns), failed on empty scope
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying -", tt, "\n", sep='')
	    flush.console()
        }
        ## FIXME: make this more robust, somehow?
        ## three choices explored so far:
        ##  (1) evaluate nfit in parent frame: tests in inst/tests/test-formulaEval.R
        ##      will fail on lapply(m_data_List,drop1)
        ##      (formula environment contains r,x,y,z but not d)
        ##  (2) evaluate nfit in frame of formula: tests will fail when data specified and formula is character
        ##  (3) update with data=NULL: fails when ...
        ##
        if (evalhack %in% c("parent","formulaenv")) {
            nfit <- update(object,
                           as.formula(paste("~ . -", tt)),
                           evaluate = FALSE)
            ## nfit <- eval(nfit, envir = env) # was  eval.parent(nfit)
            if (evalhack=="parent") {
                nfit <- eval.parent(nfit)
            } else if (evalhack=="formulaenv") {
                nfit <- eval(nfit,envir=env)
            }
        } else {
            nfit <- update(object,
                           as.formula(paste("~ . -", tt)),data=NULL,
                           evaluate = FALSE)
            nfit <- eval(nfit,envir=env)
        }
	if (test=="user") {
            ans[i+1, ] <- sumFun(object, nfit, scale=scale, k=k, ...)
        } else {
            ans[i+1, ] <- sumFun(nfit, scale, k = k, ...)
        }
        nnew <- nobs(nfit, use.fallback = TRUE)
        if(all(is.finite(c(n0, nnew))) && nnew != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    if (test=="user") {
        aod <- as.data.frame(ans)
    } else {
        dfs <- ans[1L, 1L] - ans[, 1L]
        dfs[1L] <- NA
        aod <- data.frame(Df = dfs, AIC = ans[,2])
        if(test == "Chisq") {
            ## reconstruct deviance from AIC (ugh)
            dev <- ans[, 2L] - k*ans[, 1L]
            dev <- dev - dev[1L] ; dev[1L] <- NA
            nas <- !is.na(dev)
            P <- dev
            P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
            aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
        } else if (test == "F") {
            ## FIXME: allow this if denominator df are specified externally?
            stop("F test STUB -- unfinished maybe forever")
            dev <- ans[, 2L] - k*ans[, 1L]
            dev <- dev - dev[1L] ; dev[1L] <- NA
            nas <- !is.na(dev)
            P <- dev
            P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
            aod[, c("LRT", "Pr(F)")] <- list(dev, P)
        }
    }
    head <- c("Single term deletions", "\nModel:", deparse(formula(object)),
	      if(scale > 0) paste("\nscale: ", format(scale), "\n"))
    if (!is.null(method <- attr(ss,"method"))) {
        head <- c(head,"Method: ",method,"\n")
    }
    structure(aod, heading = head, class = c("anova", "data.frame"))
}

##' @importFrom stats extractAIC
##' @S3method extractAIC merMod
extractAIC.merMod <- function(fit, scale = 0, k = 2, ...) {
    L <- logLik(refitML(fit))
    edf <- attr(L,"df")
    c(edf,-2*L + k*edf)
}

##' @importFrom stats family
##' @S3method family merMod
family.merMod <- function(object, ...) family(object@resp, ...)

##' @S3method family glmResp
family.glmResp <- function(object, ...) {
                                        # regenerate initialize
                                        # expression if necessary

    ## FIXME: may fail with user-specified/custom family?
    ## should be obsolete
    if(is.null(object$family$initialize))
        return(do.call(object$family$family,
                       list(link=object$family$link)))

    object$family
}

##' @S3method family lmResp
family.lmResp <- function(object, ...) gaussian()

##' @S3method family nlsResp
family.nlsResp <- function(object, ...) gaussian()

##' @importFrom stats fitted
##' @S3method fitted merMod
fitted.merMod <- function(object, ...) {
    xx <- object@resp$mu
    if (length(xx)==0) {
        ## handle 'fake' objects created by simulate()
        xx <- rep(NA,nrow(model.frame(object)))
    }
    if (is.null(nm <- rownames(model.frame(object)))) nm <- seq_along(xx)
    names(xx) <- nm
    if (!is.null(fit.na.action <- attr(model.frame(object),"na.action")))
        napredict(fit.na.action, xx)
    else
        xx
}

##' Extract the fixed-effects estimates
##'
##' Extract the estimates of the fixed-effects parameters from a fitted model.
##' @name fixef
##' @title Extract fixed-effects estimates
##' @aliases fixef fixed.effects fixef.merMod
##' @docType methods
##' @param object any fitted model object from which fixed effects estimates can
##' be extracted.
##' @param \dots optional additional arguments. Currently none are used in any
##' methods.
##' @return a named, numeric vector of fixed-effects estimates.
##' @keywords models
##' @examples
##' fixef(lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
##' @importFrom nlme fixef
##' @export fixef
##' @method fixef merMod
##' @export
fixef.merMod <- function(object, add.dropped=FALSE, ...) {
    X <- getME(object,"X")
    ff <- structure(object@beta, names = dimnames(X)[[2]])
    if (add.dropped) {
        if (!is.null(dd <- attr(X,"col.dropped"))) {
            ## restore positions dropped for rank deficiency
            vv <- numeric(length(ff)+length(dd))
            all.pos <- seq_along(vv)
            kept.pos <- all.pos[-dd]
            vv[kept.pos] <- ff
            names(vv)[kept.pos] <- names(ff)
            vv[dd] <- NA
            names(vv)[dd] <- names(dd)
            ff <- vv
        }
    }
    return(ff)
}


getFixedFormula <- function(form) {
    RHSForm(form) <- nobars(RHSForm(form))
    form
}

##' @importFrom stats formula
##' @S3method formula merMod
formula.merMod <- function(x, fixed.only=FALSE, random.only=FALSE, ...) {
    if (missing(fixed.only) && random.only) fixed.only <- FALSE
    if (fixed.only && random.only) stop("can't specify 'only fixed' and 'only random' terms")
    if (is.null(form <- attr(x@frame,"formula"))) {
        if (!grepl("lmer$",deparse(getCall(x)[[1]])))
            stop("can't find formula stored in model frame or call")
        form <- as.formula(formula(getCall(x),...))
    }
    if (fixed.only) {
        form <- getFixedFormula(form)
    }
    if (random.only) {
        ## from predict.R
        form <- reOnly(form,response=TRUE)
    }
    form
}

##' @S3method isREML merMod
isREML.merMod <- function(x, ...) as.logical(x@devcomp$dims[["REML"]])

##' @S3method isGLMM merMod
isGLMM.merMod <- function(x,...) {
  as.logical(x@devcomp$dims[["GLMM"]])
  ## or: is(x@resp,"glmResp")
}

##' @S3method isNLMM merMod
isNLMM.merMod <- function(x,...) {
  as.logical(x@devcomp$dims[["NLMM"]])
  ## or: is(x@resp,"nlsResp")
}

##' @S3method isLMM merMod
isLMM.merMod <- function(x,...) {
  !isGLMM(x) && !isNLMM(x)
  ## or: is(x@resp,"lmerResp") ?
}

npar.merMod <- function(object) {
    n <- length(object@beta) + length(object@theta) +
        object@devcomp[["dims"]][["useSc"]]
    ## FIXME: this is a bit of a hack: a user *might* have specified
    ## negative binomial family with a known theta, in which case we
    ## shouldn't count it as extra.  Either glmer.nb needs to set a
    ## flag somewhere, or we need class 'nbglmerMod' to extend 'glmerMod' ...
    ## We do *not* want to use the 'useSc' slot (as above), because
    ## although theta is in some sense a scale parameter, it's not
    ## one in the formal sense (and isn't stored in the 'sigma' slot)
    if (grepl("Negative Binomial",family(object)$family)) {
        n <- n+1
    }
    return(n)
    ## TODO: how do we feel about counting the scale parameter ???
}

##' @importFrom stats logLik
##' @S3method logLik merMod
logLik.merMod <- function(object, REML = NULL, ...) {
    if (is.null(REML) || is.na(REML[1]))
        REML <- isREML(object)
    val <- -devCrit(object, REML = REML)/2
    ## dc <- object@devcomp
    nobs <- nobs.merMod(object)
    structure(val,
	      nobs = nobs,
	      nall = nobs,
	      df = npar.merMod(object),
              ## length(object@beta) + length(object@theta) + dc$dims[["useSc"]],
	      class = "logLik")
}

##' @importFrom stats df.residual
##' @S3method df.residual merMod
##  TODO: not clear whether the residual df should be based
##  on p=length(beta) or p=length(c(theta,beta)) ... but
##  this is just to allow things like aods3::gof to work ...
##
df.residual.merMod <- function(object, ...) {
    nobs(object)-npar.merMod(object)
}

stripwhite <- function(x) gsub("(^ +| +$)","",x) # FIXME: never used ?
##' @importFrom stats logLik
##' @S3method model.frame merMod
model.frame.merMod <- function(formula, fixed.only=FALSE, ...) {
    fr <- formula@frame
    if (fixed.only) {
        ff <- formula(formula,fixed.only=TRUE)
        ## thanks to Thomas Leeper and Roman Lustrik, Stack Overflow
        vars <- rownames(attr(terms.formula(ff), "factors"))
        fr <- fr[vars]
    }
    fr
}

##' @importFrom stats model.matrix
##' @S3method model.matrix merMod
model.matrix.merMod <-
    function(object,
             type = c("fixed", "random", "randomListRaw"), ...) {
    switch(type[1],
           "fixed" = object@pp$X,
           "random" = getME(object, "Z"),
           "randomListRaw" = mmList(object))
}

##' Dummy variables (experimental)
##'
##' Largely a wrapper for \code{model.matrix} that
##' accepts a factor, \code{f}, and returns a dummy
##' matrix with \code{nlevels(f)-1} columns.
dummy <- function(f, levelsToKeep){
  f <- as.factor(f)
  mm <- model.matrix(~ 0 + f)
  colnames(mm) <- levels(f)

                                        # sort out levels to keep
  missingLevels <- missing(levelsToKeep)
  if(missingLevels) levelsToKeep <- levels(f)[-1]
  if(!any(levels(f) %in% levelsToKeep))
      stop("at least some of the levels in f ",
           "must also be present in levelsToKeep")
  if(!all(levelsToKeep %in% levels(f)))
      stop("all of the levelsToKeep must be levels of f")
  mm <- mm[, levelsToKeep, drop=FALSE]

  ##                                       # communicate that some usages are unlikely
  ##                                       # to help with readibility, which is the
  ##                                       # whole purpose of dummy()
  ## if((!missingLevels)&&(ncol(mm) > 1))
  ##     message("note from dummy:  explicitly specifying more than one ",
  ##             "level to keep may do little to improve readibility")
  return(mm)
}


##' @importFrom stats nobs
##' @S3method nobs merMod
nobs.merMod <- function(object, ...) nrow(object@frame)

## used in  summary.merMod():
ngrps <- function(object, ...) UseMethod("ngrps")

ngrps.default <- function(object, ...) stop("Cannot extract the number of groups from this object")

ngrps.merMod <- function(object, ...) vapply(object@flist, nlevels, 1)

ngrps.factor <- function(object, ...) nlevels(object)

##' @importFrom nlme ranef
##' @export ranef
NULL

##' Extract the modes of the random effects
##'
##' A generic function to extract the conditional modes of the random effects
##' from a fitted model object.  For linear mixed models the conditional modes
##' of the random effects are also the conditional means.
##'
##' If grouping factor i has k levels and j random effects per level the ith
##' component of the list returned by \code{ranef} is a data frame with k rows
##' and j columns.  If \code{condVar} is \code{TRUE} the \code{"postVar"}
##' attribute is an array of dimension j by j by k.  The kth face of this array
##' is a positive definite symmetric j by j matrix.  If there is only one
##' grouping factor in the model the variance-covariance matrix for the entire
##' random effects vector, conditional on the estimates of the model parameters
##' and on the data will be block diagonal and this j by j matrix is the kth
##' diagonal block.  With multiple grouping factors the faces of the
##' \code{"postVar"} attributes are still the diagonal blocks of this
##' conditional variance-covariance matrix but the matrix itself is no longer
##' block diagonal.
##' @name ranef
##' @aliases ranef ranef.merMod
##' @param object an object of a class of fitted models with random effects,
##' typically an \code{"\linkS4class{merMod}"} object.
##' @param condVar an optional logical argument indicating if the conditional
##' variance-covariance matrices of the random effects should be added as an attribute.
##' @param postVar a (deprecated) synonym for \code{condVar}
##' @param drop an optional logical argument indicating components of the return
##' value that would be data frames with a single column, usually a column
##' called \sQuote{\code{(Intercept)}}, should be returned as named vectors.
##' @param whichel an optional character vector of names of grouping factors for
##' which the random effects should be returned.  Defaults to all the grouping
##' factors.
##' @param \dots some methods for this generic function require additional
##' arguments.
##' @return A list of data frames, one for each grouping factor for the random
##' effects.  The number of rows in the data frame is the number of levels of
##' the grouping factor.  The number of columns is the dimension of the random
##' effect associated with each level of the factor.
##'
##' If \code{condVar} is \code{TRUE} each of the data frames has an attribute
##' called \code{"postVar"} which is a three-dimensional array with symmetric
##' faces.
##'
##' When \code{drop} is \code{TRUE} any components that would be data frames of
##' a single column are converted to named numeric vectors.
##' @note To produce a \dQuote{caterpillar plot} of the random effects apply
##' \code{\link[lattice:xyplot]{dotplot}} to the result of a call to
##' \code{ranef} with \code{condVar = TRUE}.
##' @examples
##' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##' fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
##' fm3 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
##' ranef(fm1)
##' str(rr1 <- ranef(fm1, condVar = TRUE))
##' dotplot(rr1)  ## default
##' ## specify free scales in order to make Day effects more visible
##' dotplot(rr1,scales = list(x = list(relation = 'free')))[["Subject"]]
##' if(FALSE) { ##-- condVar=TRUE is not yet implemented for multiple terms -- FIXME
##' str(ranef(fm2, condVar = TRUE))
##' }
##' op <- options(digits = 4)
##' ranef(fm3, drop = TRUE)
##' options(op)
##' @keywords models methods
##' @method ranef merMod
##' @export
ranef.merMod <- function(object, condVar = FALSE, drop = FALSE,
			 whichel = names(ans), postVar = FALSE, ...)
{
    if (!missing(postVar) && missing(condVar)) {
        warning(sQuote("postVar")," is deprecated: please use ",
                sQuote("condVar")," instead")
        condVar <- postVar
    }
    ans <- object@pp$b(1) ## not always == c(matrix(unlist(getME(object,"b"))))
    if (!is.null(object@flist)) {
	## evaluate the list of matrices
	levs <- lapply(fl <- object@flist, levels)
	asgn <- attr(fl, "assign")
	cnms <- object@cnms
	nc <- lengths(cnms)
	nb <- nc * lengths(levs)[asgn]
	nbseq <- rep.int(seq_along(nb), nb)
	ml <- split(ans, nbseq)
	for (i in seq_along(ml))
	    ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE,
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

	if (condVar) {
            sigsqr <- sigma(object)^2
            rp <- rePos$new(object)
            if(any(lengths(rp$terms) > 1L)) {
                # TODO: actually use condVar and then convert back to array format
                warning("conditional variances not currently available via ",
                        "ranef when there are multiple terms per factor")
            } else{
                vv <- .Call(merPredDcondVar, object@pp$ptr(), as.environment(rp))
                for (i in names(ans)) ## seq_along(ans))
                    attr(ans[[i]], "postVar") <- vv[[i]] * sigsqr
            }
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
}## ranef.merMod

print.ranef.mer <- function(x, ...) {
    print(unclass(x), ...)
    if(any(has.pv <- vapply(x, function(el)
			    !is.null(attr(el, "postVar")), NA)))
	cat('with conditional variances for',
	    paste(dQuote(names(x)[has.pv]), sep=", "), "\n")
    invisible(x)
}

## try to redo refit by calling modular structure ...
refit2.merMod <- function(object,
                          newresp=NULL) {
    ## the idea is to steal as much structure as we can from the
    ## previous fit, including
    ##  * starting parameter values
    ##  * random-effects structure
    ##  * fixed-effects structure
    ##  * model frame
    ## and jump into the modular structure at an appropriate place;
    ##  essentially, this should merge with a smart-as-possible
    ##  version of 'update' ...

}

## FIXME DRY: much of copy'n'paste from lmer() etc .. ==> become more modular (?)
refit.merMod <- function(object,
                         newresp=NULL,
                         ## formula=NULL, weights=NULL,
                         rename.response=FALSE,
                         maxit = 100L, ...)
{

    l... <- list(...)

    ctrl.arg <- NULL
    if("control" %in% names(l...)) ctrl.arg <- l...$control

    if(!all(names(l...) %in% c("control", "verbose"))) {
        warning("additional arguments to refit.merMod ignored")
    }
    ## TODO: not clear whether we should reset the names
    ##       to the new response variable.  Maybe not.

    ## retrieve name before it gets mangled by operations on newresp
    newrespSub <- substitute(newresp)

    ## for backward compatibility/functioning of refit(fit,simulate(fit))
    if (is.list(newresp)) {
        if (length(newresp)==1) {
            na.action <- attr(newresp,"na.action")
            newresp <- newresp[[1]]
            attr(newresp,"na.action") <- na.action
        } else {
            stop("refit not implemented for 'newresp' lists of length > 1: ",
                 "consider ", sQuote("lapply(object,refit)"))
        }
    }

    ## oldresp <- object@resp$y # need to set this before deep copy,
    ##                          # otherwise it gets reset with the call
    ##                          # to setResp below

    ## somewhat repeated from profile.merMod, but sufficiently
    ##  different that refactoring is slightly non-trivial
    ## "three minutes' thought would suffice ..."
    control <-
	if (!is.null(ctrl.arg)) {
	    if (length(ctrl.arg$optCtrl) == 0) { ## use object's version:
		obj.control <- object@optinfo$control
		ignore.pars <- c("xst", "xt")
		if (any(ign <- names(obj.control) %in% ignore.pars))
		    obj.control <- obj.control[!ign]
		ctrl.arg$optCtrl <- obj.control
	    }
	    ctrl.arg
	}
	else if (isGLMM(object))
	    glmerControl()
	else
	    lmerControl()

    ## we need this stuff defined before we call .glmerLaplace below ...
    pp      <- object@pp$copy()
    dc      <- object@devcomp
    nAGQ    <- dc$dims["nAGQ"] # possibly NA
    nth     <- dc$dims[["nth"]]
    verbose <- l...$verbose; if (is.null(verbose)) verbose <- 0L
    if (!is.null(newresp)) {
        ## update call and model frame with new response
	rcol <- attr(attr(model.frame(object), "terms"), "response")
        if (rename.response) {
            attr(object@frame,"formula")[[2]] <- object@call$formula[[2]] <-
                newrespSub
            names(object@frame)[rcol] <- deparse(newrespSub)
        }
        if (!is.null(na.act <- attr(object@frame,"na.action")) &&
            is.null(attr(newresp,"na.action"))) {
            ## will only get here if na.action is 'na.omit' or 'na.exclude'
            ## *and* newresp does not have an 'na.action' attribute
            ## indicating that NAs have already been filtered
            newresp <- if (is.matrix(newresp))
                newresp[-na.act, ]
            else newresp[-na.act]
        }
        object@frame[,rcol] <- newresp

        ## modFrame <- model.frame(object)
        ## modFrame[, attr(terms(modFrame), "response")] <- newresp
    }

    rr <- if(isLMM(object))
        mkRespMod(model.frame(object), REML = isREML(object))
    else if(isGLMM(object)) {
        mkRespMod(model.frame(object), family = family(object))
    } else
        stop("refit.merMod not working for nonlinear mixed models.\n",
             "try update.merMod instead.")

    if(!is.null(newresp)) {
        if(family(object)$family == "binomial") {
            ## re-do conversion of two-column matrix and factor
            ##  responses to proportion/weights format
            if (is.matrix(newresp) && ncol(newresp) == 2) {
                ntot <- rowSums(newresp)
                ## FIXME: test what happens for (0,0) rows
                newresp <- newresp[,1]/ntot
                rr$setWeights(ntot)
            }
            if (is.factor(newresp)) {
                ## FIXME: would be better to do this consistently with
                ## whatever machinery is used in glm/glm.fit/glmer ... ??
                newresp <- as.numeric(newresp)-1
            }
        }
        ## if (isGLMM(object) && rr$family$family=="binomial") {

        ## }
        stopifnot(length(newresp <- as.numeric(as.vector(newresp))) ==
                  length(rr$y))

    }


    if (isGLMM(object)) {
        GQmat <- GHrule(nAGQ)

        if (nAGQ <= 1) {
            glmerPwrssUpdate(pp,rr, control$tolPwrss, GQmat, maxit=maxit)
        } else {
            glmerPwrssUpdate(pp,rr, control$tolPwrss, GQmat, maxit=maxit,
                             grpFac = object@flist[[1]])
        }
    }

    devlist <-
	if (isGLMM(object)) {
	    baseOffset <- forceCopy(object@resp$offset)

	    list(tolPwrss= dc$cmp [["tolPwrss"]],
		 compDev = dc$dims[["compDev"]],
		 nAGQ = unname(nAGQ),
		 lp0 = pp$linPred(1), ## object@resp$eta - baseOffset,
		 baseOffset = baseOffset,
		 pwrssUpdate = glmerPwrssUpdate,
		 ## save GQmat in the object and use that instead of nAGQ
		 GQmat = GHrule(nAGQ),
		 fac = object@flist[[1]],
		 pp=pp, resp=rr, u0=pp$u0, verbose=verbose, dpars=seq_len(nth))
	} else
	    list(pp=pp, resp=rr, u0=pp$u0, verbose=verbose, dpars=seq_len(nth))
    ff <- mkdevfun(list2env(devlist), nAGQ=nAGQ, maxit=maxit, verbose=verbose)
    ## rho <- environment(ff)
    xst       <- rep.int(0.1, nth)
    x0        <- pp$theta
    lower     <- object@lower
    if (!is.na(nAGQ) && nAGQ > 0L) {
        xst   <- c(xst, sqrt(diag(pp$unsc())))
        x0    <- c(x0, unname(fixef(object)))
        lower <- c(lower, rep(-Inf,length(x0)-length(lower)))
    }
    ## control <- c(control,list(xst=0.2*xst, xt=xst*0.0001))
    ## FIX ME: allow use.last.params to be passed through
    calc.derivs <- !is.null(object@optinfo$derivs)
    ## if(isGLMM(object)) {
    ##     rho$resp$updateWts()
    ##     rho$pp$updateDecomp()
    ##     rho$lp0 <- rho$pp$linPred(1)
    ## }
    opt <- optwrap(object@optinfo$optimizer,
                   ff, x0, lower=lower, control=control$optCtrl,
                   calc.derivs=calc.derivs)
    cc <- checkConv(attr(opt,"derivs"),opt$par,
                    ## FIXME: was there a reason that ctrl was passed
                    ## via the call slot?  it was causing problems
                    ## when optTheta called refit (github issue #173)
		    # ctrl = eval(object@call$control)$checkConv,
                    ctrl = control$checkConv,
                    # ctrl = cntrl$checkConv,
                    lbound=lower)
    if (isGLMM(object)) rr$setOffset(baseOffset)
    mkMerMod(environment(ff), opt,
             list(flist=object@flist, cnms=object@cnms,
                  Gp=object@Gp, lower=object@lower),
             object@frame, getCall(object), cc)
}

refitML.merMod <- function (x, optimizer="bobyqa", ...) {
    ## FIXME: optimizer is set to 'bobyqa' for back-compatibility, but that's not
    ##  consistent with lmer (default NM).  Should be based on internally stored 'optimizer' value
    if (!isREML(x)) return(x)
    stopifnot(is(rr <- x@resp, "lmerResp"))
    rho <- new.env(parent=parent.env(environment()))
    rho$resp <- new(class(rr), y=rr$y, offset=rr$offset, weights=rr$weights, REML=0L)
    xpp <- x@pp$copy()
    rho$pp <- new(class(xpp), X=xpp$X, Zt=xpp$Zt, Lambdat=xpp$Lambdat,
                  Lind=xpp$Lind, theta=xpp$theta, n=nrow(xpp$X))
    devfun <- mkdevfun(rho, 0L) # also pass ?? (verbose, maxit, control)
    opt <- ## "smart" calc.derivs rules
	if(optimizer == "bobyqa" && !any("calc.derivs" == names(list(...))))
	    optwrap(optimizer, devfun, x@theta, lower=x@lower, calc.derivs=TRUE, ...)
	else
	    optwrap(optimizer, devfun, x@theta, lower=x@lower, ...)
    ## FIXME: Should be able to call mkMerMod() here, and be done
    n <- length(rr$y)
    pp <- rho$pp
    p <- ncol(pp$X)
    dims <- c(N=n, n=n, p=p, nmp=n-p, q=nrow(pp$Zt), nth=length(pp$theta),
	      useSc=1L, reTrms=length(x@cnms),
	      spFe=0L, REML=0L, GLMM=0L, NLMM=0L)#, nAGQ=NA_integer_)
    wrss <- rho$resp$wrss()
    ussq <- pp$sqrL(1)
    pwrss <- wrss + ussq
    cmp <- c(ldL2=pp$ldL2(), ldRX2=pp$ldRX2(), wrss=wrss, ussq=ussq,
	     pwrss=pwrss, drsum=NA, dev=opt$fval, REML=NA,
	     sigmaML=sqrt(pwrss/n), sigmaREML=sqrt(pwrss/(n-p)))
    ## modify the call  to have REML=FALSE. (without evaluating the call!)
    cl <- x@call
    cl[["REML"]] <- FALSE
    new("lmerMod", call = cl, frame=x@frame, flist=x@flist,
	cnms=x@cnms, theta=pp$theta, beta=pp$delb, u=pp$delu,
	optinfo = .optinfo(opt),
	lower=x@lower, devcomp=list(cmp=cmp, dims=dims), pp=pp, resp=rho$resp)
}

##' residuals of merMod objects 		--> ../man/residuals.merMod.Rd
##' @param object a fitted [g]lmer (\code{merMod}) object
##' @param type type of residuals
##' @param scaled scale residuals by residual standard deviation (=scale parameter)?
##' @param \dots additional arguments (ignored: for method compatibility)
residuals.merMod <- function(object,
                             type = if(isGLMM(object)) "deviance" else "response",
                             scaled = FALSE, ...) {
    r <- residuals(object@resp, type,...)
    fr <- model.frame(object)
    if (is.null(nm <- rownames(fr))) nm <- seq_along(r)
    names(r) <- nm
    if (scaled) r <- r/sigma(object)
    if (!is.null(na.action <- attr(fr, "na.action")))
        naresid(na.action, r)
    else
        r
}

##' @rdname residuals.merMod
##' @S3method residuals lmResp
##' @method residuals lmResp
residuals.lmResp <- function(object,
                             type = c("working", "response", "deviance",
                                      "pearson", "partial"), ...) {
    y <- object$y
    r <- object$wtres
    mu <- object$mu
    switch(match.arg(type),
           working =,
           response = y-mu,
           deviance =,
           pearson = r,
           partial = stop(gettextf("partial residuals are not implemented yet"),
                          call. = FALSE)
           )
}

##' @rdname residuals.merMod
##' @S3method residuals glmResp
##' @method residuals glmResp
residuals.glmResp <- function(object, type = c("deviance", "pearson",
                                      "working", "response", "partial"),
                              ...) {
    type <- match.arg(type)
    y <- object$y
    mu <- object$mu
    switch(type,
           deviance = {
               d.res <- sqrt(object$devResid())
               ifelse(y > mu, d.res, -d.res)
           },
           pearson = object$wtres,
           working = object$wrkResids(),
           response = y - mu,
           partial = stop(gettextf("partial residuals are not implemented yet"),
                    call. = FALSE)
       )

}

hatvalues.merMod <- function(model, fullHatMatrix = FALSE, ...) {
    if(isGLMM(model)) warning("the hat matrix may not make sense for GLMMs")
    ## FIXME:  add restriction for NLMMs?

    ## prior weights, W ^ {1/2} :
    sqrtW <- Diagonal(x = sqrt(weights(model, type = "prior")))
    with(getME(model, c("L", "Lambdat", "Zt", "RX", "X", "RZX")), {
        ## CL:= right factor of the random-effects component of the hat matrix	(64)
        CL <- solve(L, solve(L, Lambdat %*% Zt %*% sqrtW,
                             system = "P"), system = "L")
	## CR:= right factor of the fixed-effects component of the hat matrix	(65)
	##	{MM (FIXME Matrix):  t(.) %*% here faster than crossprod()}
	CR <- solve(t(RX), t(X) %*% sqrtW - crossprod(RZX, CL))
	if(fullHatMatrix) ## H = (C_L^T C_L + C_R^T C_R)			(63)
	    crossprod(CL) + crossprod(CR)
	else ## diagonal of the hat matrix, diag(H) :
	    colSums(CR^2) + colSums(CL^2)
    })
}


##' @S3method sigma merMod
sigma.merMod <- function(object, ...) {
    dc <- object@devcomp
    dd <- dc$dims
    if(dd[["useSc"]])
        dc$cmp[[if(dd[["REML"]]) "sigmaREML" else "sigmaML"]] else 1.
}

##' @importFrom stats terms
##' @S3method terms merMod
terms.merMod <- function(x, fixed.only=TRUE, random.only=FALSE, ...) {
    if (missing(fixed.only) && random.only) fixed.only <- FALSE
    if (fixed.only && random.only) stop("can't specify 'only fixed' and 'only random' terms")
    tt <- attr(x@frame,"terms")
    if (fixed.only) {
        tt <- terms.formula(formula(x,fixed.only=TRUE))
        attr(tt,"predvars") <- attr(terms(x@frame),"predvars.fixed")
    }
    if (random.only) {
        tt <- terms.formula(subbars(formula(x,random.only=TRUE)))
        ## FIXME: predvars should be random-only
        attr(tt,"predvars") <- attr(terms(x@frame),"predvars.random")
    }
    tt
}

##' @importFrom stats update
##' @S3method update merMod
update.merMod <- function(object, formula., ..., evaluate = TRUE) {
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
    if (evaluate) {
        ff <- environment(formula(object))
        pf <- parent.frame()  ## save parent frame in case we need it
        sf <- sys.frames()[[1]]
        tryCatch(eval(call, envir=ff),
                 error=function(e) {
                     tryCatch(eval(call, envir=sf),
                              error=function(e) {
                                  eval(call, pf)
                              })
                 })
    } else call
}

###----- Printing etc ----------------------------

## lme4.0, for GLMM had
## 'Generalized linear mixed model fit by the Laplace approximation'
## 'Generalized linear mixed model fit by the adaptive Gaussian Hermite approximation'
## so did *not* mention  "maximum likelihood" at all in the GLMM case


methTitle <- function(dims) { # dims == object@devcomp$dims
    GLMM <- dims[["GLMM"]]
    kind <- switch(1L + GLMM * 2L + dims[["NLMM"]],
		   "Linear", "Nonlinear",
		   "Generalized linear", "Generalized nonlinear")
    paste(kind, "mixed model fit by",
	  if(dims[["REML"]]) "REML"
	  else paste("maximum likelihood",
		     if(GLMM) {
			 ## TODO? Use shorter wording here, for (new) 'long = FALSE' argument
			 if((nAGQ <- dims[["nAGQ"]]) == 1)
			     "(Laplace Approximation)"
			 else
			     sprintf("(Adaptive Gauss-Hermite Quadrature, nAGQ = %d)",
				     nAGQ)
			 }))
}


cat.f <- function(...) cat(..., fill = TRUE)

famlink <- function(object, resp = object@resp) {
    if(is(resp, "glmResp"))
	resp$family[c("family", "link")]
    else list(family = NULL, link = NULL)
}

##' @title print method title
##' @param mtit the result of methTitle(obj)
##' @param class typically class(obj)
.prt.methTit <- function(mtit, class) {
    if(nchar(mtit) + 5 + nchar(class) > (w <- getOption("width"))) {
	## wrap around
	mtit <- strwrap(mtit, width = w - 2, exdent = 2)
	cat(mtit, " [",class,"]", sep = "", fill = TRUE)
    } else ## previous: simple one-liner
	cat(sprintf("%s ['%s']\n", mtit, class))
}

.prt.family <- function(famL) {
    if (!is.null(f <- famL$family)) {
	cat.f(" Family:", f,
	      if(!is.null(ll <- famL$link)) paste(" (", ll, ")"))
    }
}

.prt.resids <- function(resids, digits, title = "Scaled residuals:", ...) {
    cat(title,"\n")
    ## FIXME: need testing code
    rq <- setNames(zapsmall(quantile(resids, na.rm=TRUE), digits + 1L),
                   c("Min", "1Q", "Median", "3Q", "Max"))
    print(rq, digits = digits, ...)
    cat("\n")
}

.prt.call <- function(call, long = TRUE) {
    if (!is.null(cc <- call$formula))
	cat.f("Formula:", deparse(cc))
    if (!is.null(cc <- call$data))
	cat.f("   Data:", deparse(cc))
    if (!is.null(cc <- call$weights))
        cat.f("Weights:", deparse(cc))
    if (!is.null(cc <- call$offset))
        cat.f(" Offset:", deparse(cc))
    if (long && length(cc <- call$control) &&
	!identical((dc <- deparse(cc)), "lmerControl()"))
	## && !identical(eval(cc), lmerControl()))
	cat.f("Control:", dc)
    if (!is.null(cc <- call$subset))
	cat.f(" Subset:", deparse(cc))
}

##' @title Extract Log Likelihood, AIC, and related statics from a Fitted LMM
##' @param object a LMM model fit
##' @param devianceFUN the function to be used for computing the deviance; should not be changed for  \pkg{lme4} created objects.
##' @param chkREML optional logical indicating \code{object} maybe a REML fit.
##' @param devcomp
##' @return a list with components \item{logLik} and \code{AICtab} where the first is \code{\link{logLik(object)}} and \code{AICtab} is a "table" of AIC, BIC, logLik, deviance, df.residual() values.
llikAIC <- function(object, devianceFUN = devCrit, chkREML = TRUE, devcomp = object@devcomp) {
    llik <- logLik(object)   # returns NA for a REML fit - maybe change?
    AICstats <- {
	if(chkREML && devcomp$dims[["REML"]])
            devcomp$cmp["REML"] # *no* likelihood stats here
	else {
	    c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
	      deviance = devianceFUN(object),
              df.resid = df.residual(object))
	}
    }
    list(logLik = llik, AICtab = AICstats)
}

.prt.aictab <- function(aictab, digits = 1) {
    t.4 <- round(aictab, digits)
    if (length(aictab) == 1 && names(aictab) == "REML")
	cat.f("REML criterion at convergence:", t.4)
    else {
        ## slight hack to get residual df formatted as an integer
        t.4F <- format(t.4)
        t.4F["df.resid"] <- format(t.4["df.resid"])
        print(t.4F, quote = FALSE)
    }
}

.prt.VC <- function(varcor, digits, comp, formatter = format, ...) {
    cat("Random effects:\n")
    fVC <- if(missing(comp))
	formatVC(varcor, digits = digits, formatter = formatter)
    else
	formatVC(varcor, digits = digits, formatter = formatter, comp = comp)
    print(fVC, quote = FALSE, digits = digits, ...)
}

.prt.grps <- function(ngrps, nobs) {
    cat(sprintf("Number of obs: %d, groups: ", nobs),
        paste(paste(names(ngrps), ngrps, sep = ", "), collapse = "; "),
        fill = TRUE)
}

## FIXME: print header ("Warnings:\n") ?
##  change position in output? comes at the very end, could get lost ...
.prt.warn <- function(optinfo, summary=FALSE, ...) {
    if(length(optinfo) == 0) return() # currently, e.g., from refitML()
    ## check all warning slots: print numbers of warnings (if any)
    cc <- optinfo$conv$opt
    msgs <- unlist(optinfo$conv$lme4$messages)
    ## can't put nmsgs/nwarnings compactly into || expression
    ##   because of short-circuiting
    nmsgs <- length(msgs)
    warnings <- optinfo$warnings
    nwarnings <- length(warnings)
    if (cc > 0 || nmsgs > 0 || nwarnings > 0) {
        if (summary) {
            cat(sprintf("convergence code %d; %d optimizer warnings; %d lme4 warnings",
                cc,nmsgs,nwarnings),"\n")
        } else {
            cat(sprintf("convergence code: %d", cc),
                msgs,
                unlist(warnings),
                sep="\n")
            cat("\n")
        }
    }
}

##  options(lme4.summary.cor.max = 20)  --> ./hooks.R
##                                           ~~~~~~~~
## was      .summary.cor.max <- 20    a lme4-namespace hidden global variable

## This is modeled a bit after	print.summary.lm :
## Prints *both*  'mer' and 'merenv' - as it uses summary(x) mainly
##' @S3method print summary.merMod
print.summary.merMod <- function(x, digits = max(3, getOption("digits") - 3),
                                 correlation = NULL, symbolic.cor = FALSE,
                                 signif.stars = getOption("show.signif.stars"),
                                 ranef.comp = c("Variance", "Std.Dev."),
                                 show.resids = TRUE, ...)
{
    .prt.methTit(x$methTitle, x$objClass)
    .prt.family(x)
    .prt.call(x$call); cat("\n")
    .prt.aictab(x$AICtab); cat("\n")
    if (show.resids)
        ## need residuals.merMod() rather than residuals():
        ##  summary.merMod has no residuals method
        .prt.resids(x$residuals, digits = digits)
    .prt.VC(x$varcor, digits = digits, useScale = x$useScale,
	    comp = ranef.comp, ...)
    .prt.grps(x$ngrps, nobs = x$devcomp$dims[["n"]])

    p <- nrow(x$coefficients)
    if (p > 0) {
	cat("\nFixed effects:\n")
	printCoefmat(x$coefficients, zap.ind = 3, #, tst.ind = 4
		     digits = digits, signif.stars = signif.stars)
        ## do not show correlation when   summary(*, correlation=FALSE)  was used:
        hasCor <- !is.null(VC <- x$vcov) && !is.null(VC@factors$correlation)
	if(is.null(correlation)) { # default
	    cor.max <- getOption("lme4.summary.cor.max")
	    correlation <- hasCor && p <= cor.max
	    if(!correlation && p > cor.max) {
		nam <- deparse(substitute(x))
		if(length(nam) > 1 || nchar(nam) >= 32) nam <- "...."
		message(sprintf(paste(
		    "\nCorrelation matrix not shown by default, as p = %d > %d.",
		    "Use print(%s, correlation=TRUE)  or",
		    "	 vcov(%s)	 if you need it\n", sep = "\n"),
				p, cor.max, nam, nam))
	    }
	}
	else if(!is.logical(correlation)) stop("'correlation' must be NULL or logical")
	if(correlation) {
	    if(is.null(VC)) VC <- vcov(x, correlation = TRUE)
	    corF <- VC@factors$correlation
	    if (is.null(corF)) { # can this still happen?
		message("\nCorrelation of fixed effects could have been required in summary()")
		corF <- cov2cor(VC)
	    }
	    p <- ncol(corF)
	    if (p > 1) {
		rn <- rownames(x$coefficients)
		rns <- abbreviate(rn, minlength = 11)
		cat("\nCorrelation of Fixed Effects:\n")
		if (is.logical(symbolic.cor) && symbolic.cor) {
		    corf <- as(corF, "matrix")
		    dimnames(corf) <- list(rns,
					   abbreviate(rn, minlength = 1, strict = TRUE))
		    print(symnum(corf))
		} else {
		    corf <- matrix(format(round(corF@x, 3), nsmall = 3),
				   ncol = p,
				   dimnames = list(rns, abbreviate(rn, minlength = 6)))
		    corf[!lower.tri(corf)] <- ""
		    print(corf[-1, -p, drop = FALSE], quote = FALSE)
		} ## !symbolic.cor
	    }  ## if (p > 1)
        } ## if (correlation)
    } ## if (p>0)

    if(length(x$fitMsgs) && any(nchar(x$fitMsgs) > 0)) {
        cat("fit warnings:\n"); writeLines(x$fitMsgs)
    }
    .prt.warn(x$optinfo,summary=FALSE)
    invisible(x)
}## print.summary.merMod


##' @S3method print merMod
print.merMod <- function(x, digits = max(3, getOption("digits") - 3),
                         correlation = NULL, symbolic.cor = FALSE,
                         signif.stars = getOption("show.signif.stars"),
			 ranef.comp = "Std.Dev.", ...)
{
    dims <- x@devcomp$dims
    .prt.methTit(methTitle(dims), class(x))
    .prt.family(famlink(x, resp = x@resp))
    .prt.call(x@call, long = FALSE)
    ## useScale <- as.logical(dims[["useSc"]])

    llAIC <- llikAIC(x)
    .prt.aictab(llAIC$AICtab, 4)
    varcor <- VarCorr(x)
    .prt.VC(varcor, digits = digits, comp = ranef.comp, ...)
    ngrps <- vapply(x@flist, nlevels, 0L)
    .prt.grps(ngrps, nobs = dims[["n"]])
    if(length(cf <- fixef(x)) > 0) {
	cat("Fixed Effects:\n")
	print.default(format(cf, digits = digits),
		      print.gap = 2L, quote = FALSE, ...)
    } else cat("No fixed effect coefficients\n")

    fitMsgs <- .merMod.msgs(x)
    if(any(nchar(fitMsgs) > 0)) {
        cat("fit warnings:\n"); writeLines(fitMsgs)
    }
    .prt.warn(x@optinfo,summary=TRUE)
    invisible(x)
}

##' @exportMethod show
setMethod("show",  "merMod", function(object) print.merMod(object))

##' Return the deviance component list
devcomp <- function(x) {
    .Deprecated("getME(., \"devcomp\")")
    stopifnot(is(x, "merMod"))
    x@devcomp
}

##' @exportMethod getL
setMethod("getL", "merMod", function(x) {
    .Deprecated("getME(., \"L\")")
    getME(x, "L")
})

##' used by tnames
mkPfun <- function(diag.only = FALSE, old = TRUE, prefix = NULL){
    local({
        function(g,e) {
            mm <- outer(e,e,paste,sep = ".")
            if(old) {
                diag(mm) <- e
            } else {
                mm[] <- paste(mm,g,sep = "|")
                if (!is.null(prefix)) mm[] <- paste(prefix[2],mm,sep = "_")
                diag(mm) <- paste(e,g,sep = "|")
                if (!is.null(prefix))  diag(mm) <- paste(prefix[1],diag(mm),sep = "_")
            }
            mm <- if (diag.only) diag(mm) else mm[lower.tri(mm,diag = TRUE)]
            if(old) paste(g,mm,sep = ".") else mm
        }
    })
}

##' Construct names of individual theta/sd:cor components
##'
##' @param object a fitted model
##' @param diag.only include only diagonal elements?
##' @param old (logical) give backward-compatible results?
##' @param prefix a character vector with two elements giving the prefix
##' for diagonal (e.g. "sd") and off-diagonal (e.g. "cor") elements
tnames <- function(object,diag.only = FALSE,old = TRUE,prefix = NULL) {
    pfun <- mkPfun(diag.only = diag.only, old = old, prefix = prefix)
    c(unlist(mapply(pfun, names(object@cnms), object@cnms)))
}

## -> ../man/getME.Rd
getME <- function(object, name, ...) UseMethod("getME")


##' Extract or Get Generalize Components from a Fitted Mixed Effects Model
getME.merMod <- function(object,
		  name = c("X", "Z","Zt", "Ztlist", "mmList",
                  "y", "mu", "u", "b",
		  "Gp", "Tp",
		  "L", "Lambda", "Lambdat", "Lind", "Tlist", "A",
		  "RX", "RZX", "sigma",
                  "flist",
                  "fixef", "beta", "theta", "ST",
		  "REML", "is_REML",
                  "n_rtrms", "n_rfacs",
                  "N", "n", "p", "q",
                  "p_i", "l_i", "q_i", "k", "m_i", "m",
                  "cnms",
                  "devcomp", "offset", "lower",
                  "devfun",
                  "glmer.nb.theta"
                  ), ...)
{
    if(missing(name)) stop("'name' must not be missing")
    ## Deal with multiple names -- "FIXME" is inefficiently redoing things
    if (length(name <- as.character(name)) > 1) {
        names(name) <- name
        return(lapply(name, getME, object = object))
    }
    if(name == "ALL") ## recursively get all provided components
        return(sapply(eval(formals()$name),
                      getME.merMod, object=object, simplify=FALSE))

    stopifnot(is(object,"merMod"))
    name <- match.arg(name)
    rsp  <- object@resp
    PR   <- object@pp
    dc   <- object@devcomp
    th   <- object@theta
    ## cmp  <- dc $ cmp
    cnms <- object@cnms
    dims <- dc $ dims
    Tpfun <- function(cnms) {
	ltsize <- function(n) n*(n+1)/2 # lower triangle size
	cLen <- cumsum(ltsize(lengths(cnms)))
	setNames(c(0, cLen),
		 c("beg__", names(cnms))) ## such that diff( Tp ) is well-named
    }
    ## mmList. <- mmList(object)
    if(any(name == c("p_i", "q_i", "m_i")))
	p_i <- vapply(mmList(object), ncol, 1L)
    ## l_i <- vapply(object@flist, nlevels, 1L)
    switch(name,
	   "X" = PR$X, ## ok ? - check -- use model.matrix() method instead?
	   "Z" = t(PR$Zt),
	   "Zt" = PR$Zt,
           "Ztlist" =
           {
               getInds <- function(i) {
                   n <- diff(object@Gp)[i] ## number of elements in this block
                   nt <- length(cnms[[i]]) ## number of REs
                   inds <- lapply(seq(nt), seq, to = n, by = nt)  ## pull out individual RE indices
                   inds <- lapply(inds,function(x) x + object@Gp[i])  ## add group offset
               }
               inds <- do.call(c,lapply(seq_along(cnms),getInds))
               setNames(lapply(inds,function(i) PR$Zt[i,]),
                        tnames(object,diag.only = TRUE))
           },
	   "mmList" = mmList.merMod(object),
           "y" = rsp$y,
           "mu" = rsp$mu,
           "u" = object@u,
	   "b" = crossprod(PR$Lambdat, object@u), # == Lambda %*% u
	   "L" = PR$ L(),
	   "Lambda" = t(PR$ Lambdat),
	   "Lambdat" = PR$ Lambdat,
           "A" = PR$Lambdat %*% PR$Zt,
           "Lind" = PR$ Lind,
	   "RX" = structure(PR$RX(), dimnames = list(colnames(PR$X), colnames(PR$X))), ## maybe add names elsewhere?
	   "RZX" = structure(PR$RZX, dimnames = list(NULL, colnames(PR$X))), ## maybe add names elsewhere?
           "sigma" = sigma(object),
           "Gp" = object@Gp,
           "Tp" = Tpfun(cnms), # "term-wise theta pointer"
           "flist" = object@flist,
           "fixef" = fixef(object),
	   "beta"  = object@beta,
	   "theta" = setNames(th, tnames(object)),
	   "ST" = setNames(vec2STlist(object@theta, n = lengths(cnms)),
			   names(cnms)),
           "Tlist" = {
               nc <- lengths(cnms) # number of columns per term
               nt <- length(nc)    # number of random-effects terms
               ans <- vector("list",nt)
               names(ans) <- names(nc)
               pos <- 0L
               for (i in 1:nt) { # will need modification for more general Lambda
                   nci <- nc[i]
                   tt <- matrix(0., nci, nci)
                   inds <- lower.tri(tt,diag=TRUE)
                   nthi <- sum(inds)
                   tt[lower.tri(tt, diag=TRUE)] <- th[pos + (1L:nthi)]
                   pos <- pos + nthi
                   ans[[i]] <- tt
               }
               ans
           },
	   "REML" = dims[["REML"]],
	   "is_REML" = isREML(object),
           ## number of random-effects terms
	   "n_rtrms" = length(cnms),
           ## number of random-effects grouping factors
           "n_rfacs" = length(object@flist),
           "N" = dims["N"],
           "n" = dims["n"],
           "p" = dims["p"],
           "q" = dims["q"],
           "p_i" = p_i,
	   "l_i" = vapply(object@flist, nlevels, 1L),
	   "q_i" = ## == p_i * l_i
	       p_i*vapply(object@flist, nlevels, 1L),
           "k" = length(cnms),
           "m_i" = choose(p_i + 1, 2),
           "m" = dc$dims[["nth"]],
           "cnms" = cnms,
           "devcomp" = dc,
           "offset" = rsp$offset,
           "lower" = object@lower,
           "devfun" = {
               ## copied from refit ... DRY ...
               verbose <- getCall(object)$verbose; if (is.null(verbose)) verbose <- 0L
               devlist <-
                   if (isGLMM(object)) {
                       stop("getME('devfun') not yet working for GLMMs")## FIXME
                       baseOffset <- rsp$offset
                       nAGQ <- unname(dc$dims["nAGQ"])
                       list(tolPwrss= dc$cmp ["tolPwrss"],
                            compDev = dc$dims["compDev"],
                            nAGQ = nAGQ,
                            lp0  = rsp$eta - baseOffset,
                            baseOffset  = baseOffset,
                            pwrssUpdate = glmerPwrssUpdate,
                            GQmat = GHrule(nAGQ),
                            fac = object@flist[[1]],
                            pp=PR, resp=rsp, u0=PR$u0, dpars=seq_along(PR$theta),
                            verbose=verbose)
                   }
                   else
                       list(pp=PR, resp=rsp, u0=PR$u0, dpars=seq_along(PR$theta), verbose=verbose)
               mkdevfun(rho=list2env(devlist),
                        ## FIXME: fragile ...
                        verbose=verbose, control=object@optinfo$control)
           },
           ## FIXME: current version gives lower bounds for theta parameters only:
           ## -- these must be extended for [GN]LMMs -- give extended value including -Inf values for beta values?

	   "glmer.nb.theta" = if(isGLMM(object) && isNBfamily(rsp$family$family))
	       getNBdisp(object) else NA,

	   "..foo.." = # placeholder!
	   stop(gettextf("'%s' is not implemented yet",
			 sprintf("getME(*, \"%s\")", name))),
	   ## otherwise
	   stop(sprintf("Mixed-Effects extraction of '%s' is not available for class \"%s\"",
			name, class(object))))
}## {getME}

##' @importMethodsFrom Matrix t %*% crossprod diag tcrossprod
##' @importClassesFrom Matrix dgCMatrix dpoMatrix corMatrix
NULL

## Extract the conditional variance-covariance matrix of the fixed-effects
## parameters
vcov.merMod <- function(object, correlation = TRUE, sigm = sigma(object),
                        use.hessian = NULL, ...)
{
    hess.avail <- (!is.null(h <- object@optinfo$derivs$Hessian) &&
                   nrow(h) > (ntheta <- length(getME(object,"theta"))))
    if (is.null(use.hessian)) use.hessian <- hess.avail
    if (use.hessian && !hess.avail) stop(shQuote("use.hessian"),
                                         "=TRUE specified, ",
                                         "but Hessian is unavailable")
    calc.vcov.hess <- function(h) {
	## ~= forceSymmetric(solve(h/2)[i,i]) : solve(h/2) = 2*solve(h)
        h <- tryCatch(solve(h),
                      error=function(e) matrix(NA,nrow=nrow(h),ncol=ncol(h)))
        i <- -seq_len(ntheta)
	h <- h[i,i]
	forceSymmetric(h + t(h))
    }

    ## OBSOLETE??  checks for symmetry, but symmetry is now forced
    ## within calc.vcov.hess anyway
    ## symmetrize <-  function(v,warnTol = 0,stopTol = sqrt(.Machine$double.eps)) {
    ##     if(nrow(v) == 1L) return(v)       # 1-by-1 matrices are always symmetrical
    ##     nonSymm <- max(abs(v[lower.tri(v)]-t(v)[lower.tri(v)]))
    ##     warnTol <- max(warnTol,stopTol)
    ##     if (nonSymm > stopTol) stop(sprintf(
    ##         "calculated variance-covariance matrix is non-symmetric (tol=%f)",
    ##                                       stopTol))
    ##     if (nonSymm > warnTol) stop(sprintf(
    ##         "calculated variance-covariance matrix is non-symmetric (tol=%f)",
    ##                                       warnTol))
    ##     v[upper.tri(v)] <- t(v)[upper.tri(v)]
    ##     v
    ## }

    V <- sigm^2 * object@pp$unsc()

    if (hess.avail) {
        V.hess <- calc.vcov.hess(h)
        bad.V.hess <- any(is.na(V.hess))
        if (!bad.V.hess) {
            e.hess <- eigen(V.hess,symmetric = TRUE,only.values = TRUE)$values
            if (min(e.hess) <= 0) bad.V.hess <- TRUE
        }
    }
    if (!use.hessian && hess.avail) {
        ## if hessian is available, go ahead and check
        ## for similarity with the RX-based estimate
        ## (inverting the hessian isn't *too* expensive)
        var.hess.tol <- 1e-4 # FIXME: should var.hess.tol be user controlled?
        if (!bad.V.hess && any(abs(V-V.hess) > var.hess.tol * V.hess))
            warning("variance-covariance matrix computed ",
                    "from finite-difference Hessian\nand ",
                    "from RX differ by >",var.hess.tol,": ",
                    "consider ",shQuote("use.hessian=TRUE"))
    }
    if (use.hessian) {
        if (!bad.V.hess) {
            V <- V.hess
        } else {
            warning("variance-covariance matrix computed ",
                    "from finite-difference Hessian is\n",
                    "not positive definite or contains NA values: falling back to ",
                    "var-cov estimated from RX")
        }
    }

    ## FIXME: try to catch non-PD matrices
    rr <- tryCatch(as(V, "dpoMatrix"), error = function(e)e)
    if (inherits(rr, "error")) {
	warning(gettextf("Computed variance-covariance matrix problem: %s;\nreturning NA matrix",
                         rr$message), domain = NA)
        rr <- matrix(NA,nrow(V),ncol(V))
    }

    nmsX <- colnames(object@pp$X)
    dimnames(rr) <- list(nmsX,nmsX)

    if(correlation)
	rr@factors$correlation <-
	    if(!is.na(sigm)) as(rr, "corMatrix") else rr # (is NA anyway)
    rr
}

##' @importFrom stats vcov
##' @S3method vcov summary.merMod
vcov.summary.merMod <- function(object, ...) {
    if(is.null(object$vcov)) stop("logic error in summary of merMod object")
    object$vcov
}

##' Make variance and correlation matrices from \code{theta}
##'
##' @param sc scale factor (residual standard deviation)
##' @param cnms component names
##' @param nc numeric vector: number of terms in each RE component
##' @param theta theta vector (lower-triangle of Cholesky factors)
##' @param nms component names (FIXME: nms/cnms redundant: nms=names(cnms)?)
##' @seealso \code{\link{VarCorr}}
##' @return A matrix
##' @export
mkVarCorr <- function(sc, cnms, nc, theta, nms) {
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
    if(!all(nms == names(cnms))) ## the above FIXME
	warning("nms != names(cnms)  -- whereas lme4-authors thought they were --\n",
		"Please report!", immediate. = TRUE)
    ans <- lapply(ncseq, function(i)
	      {
		  ## Li := \Lambda_i, the i-th block diagonal of \Lambda(\theta)
		  Li <- diag(nrow = nc[i])
		  Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
		  rownames(Li) <- cnms[[i]]
		  ## val := \Sigma_i = \sigma^2 \Lambda_i \Lambda_i', the
		  val <- tcrossprod(sc * Li) # variance-covariance
		  stddev <- sqrt(diag(val))
		  corr <- t(val / stddev)/stddev
		  diag(corr) <- 1
		  structure(val, stddev = stddev, correlation = corr)
	      })
    if(is.character(nms)) {
	## FIXME: do we want this?  Maybe not.
	## Potential problem: the names of the elements of the VarCorr() list
	##  are not necessarily unique (e.g. fm2 from example("lmer") has *two*
	##  Subject terms, so the names are "Subject", "Subject".  The print method
	##  for VarCorrs handles this just fine, but it's a little awkward if we
	##  want to dig out elements of the VarCorr list ... ???
	if (anyDuplicated(nms))
	    nms <- make.names(nms, unique = TRUE)
	names(ans) <- nms
    }
    structure(ans, sc = sc)
}

##' Extract variance and correlation components
##'
VarCorr.merMod <- function(x, sigma = 1, ...)
{
  ## TODO: now that we have '...', add  type=c("varcov","sdcorr","logs" ?
    if (is.null(cnms <- x@cnms))
	stop("VarCorr methods require reTrms, not just reModule")
    if(missing(sigma))
	sigma <- sigma(x)
    nc <- lengths(cnms) # no. of columns per term
    structure(mkVarCorr(sigma, cnms = cnms, nc = nc, theta = x@theta,
			nms = { fl <- x@flist; names(fl)[attr(fl, "assign")]}),
	      useSc = as.logical(x@devcomp$dims[["useSc"]]),
	      class = "VarCorr.merMod")
}

if(FALSE)## *NOWHERE* used _FIXME_ ??
## Compute standard errors of fixed effects from an merMod object
##
## @title Standard errors of fixed effects
## @param object "merMod" object,
## @param ... additional, optional arguments.  None are used at present.
## @return numeric vector of length length(fixef(.))
unscaledVar <- function(object, ...) {
    stopifnot(is(object, "merMod"))
    sigma(object) * diag(object@pp$unsc())
}

##' @S3method print VarCorr.merMod
print.VarCorr.merMod <- function(x, digits = max(3, getOption("digits") - 2),
		   comp = "Std.Dev.", formatter = format, ...) {
    print(formatVC(x, digits = digits, comp = comp, formatter = formatter), quote = FALSE, ...)
    invisible(x)
}

##' __NOT YET EXPORTED__
##' "format()" the 'VarCorr' matrix of the random effects -- for
##' print()ing and show()ing
##'
##' @title Format the 'VarCorr' Matrix of Random Effects
##' @param varc a \code{\link{VarCorr}} (-like) matrix with attributes.
##' @param digits the number of significant digits.
##' @param comp character vector of length one or two indicating which
##' columns out of "Variance" and "Std.Dev." should be shown in the
##' formatted output.
##' @param formatter the \code{\link{function}} to be used for
##' formatting the standard deviations and or variances (but
##' \emph{not} the correlations which (currently) are always formatted
##' as "0.nnn"
##' @param ... optional arguments for \code{formatter(*)} in addition
##' to the first (numeric vector) and \code{digits}.
##' @return a character matrix of formatted VarCorr entries from \code{varc}.
formatVC <- function(varcor, digits = max(3, getOption("digits") - 2),
		     comp = "Std.Dev.", formatter = format,
                     useScale = attr(varcor, "useSc"),
                     ...)
{
    c.nms <- c("Groups", "Name", "Variance", "Std.Dev.")
    avail.c <- c.nms[-(1:2)]
    if(anyNA(mcc <- pmatch(comp, avail.c)))
	stop("Illegal 'comp': ", comp[is.na(mcc)])
    nc <- length(colnms <- c(c.nms[1:2], (use.c <- avail.c[mcc])))
    if(length(use.c) == 0)
	stop("Must *either* show variances or standard deviations")
    reStdDev <- c(lapply(varcor, attr, "stddev"),
		  if(useScale) list(Residual = unname(attr(varcor, "sc"))))
    reLens <- lengths(reStdDev)
    nr <- sum(reLens)
    reMat <- array('', c(nr, nc), list(rep.int('', nr), colnms))
    reMat[1+cumsum(reLens)-reLens, "Groups"] <- names(reLens)
    reMat[,"Name"] <- c(unlist(lapply(varcor, colnames)), if(useScale) "")
    if(any("Variance" == use.c))
	reMat[,"Variance"] <- formatter(unlist(reStdDev)^2, digits = digits, ...)
    if(any("Std.Dev." == use.c))
	reMat[,"Std.Dev."] <- formatter(unlist(reStdDev),   digits = digits, ...)
    if (any(reLens > 1)) {
	maxlen <- max(reLens)
	recorr <- lapply(varcor, attr, "correlation")
	corr <-
	    do.call(rBind,
		    lapply(recorr,
			   function(x) {
			       x <- as(x, "matrix")
			       dig <- max(2, digits - 2) # use 'digits' !
                               ## not using formatter() for correlations
			       cc <- format(round(x, dig), nsmall = dig)
			       cc[!lower.tri(cc)] <- ""
			       nr <- nrow(cc)
			       if (nr >= maxlen) return(cc)
			       cbind(cc, matrix("", nr, maxlen-nr))
			   }))[, -maxlen, drop = FALSE]
	if (nrow(corr) < nrow(reMat))
	    corr <- rbind(corr, matrix("", nrow(reMat) - nrow(corr), ncol(corr)))
	colnames(corr) <- c("Corr", rep.int("", max(0L, ncol(corr)-1L)))
	cbind(reMat, corr)
    } else reMat
}

##' @S3method summary merMod
summary.merMod <- function(object,
                           correlation = (p <= getOption("lme4.summary.cor.max")),
                           use.hessian = NULL,
                           ...)
{
    if (length(list(...)) > 0) {
        ## FIXME: need testing code
        warning("additional arguments ignored")
    }
    ## se.calc:
    hess.avail <- (!is.null(h <- object@optinfo$derivs$Hessian) &&
                   nrow(h) > length(getME(object,"theta")))
    if (is.null(use.hessian)) use.hessian <- hess.avail
    if (use.hessian && !hess.avail)
	stop("'use.hessian=TRUE' specified, but Hessian is unavailable")
    resp <- object@resp
    devC <- object@devcomp
    dd <- devC$dims
    ## cmp <- devC$cmp
    useSc <- as.logical(dd[["useSc"]])
    sig <- sigma(object)
    ## REML <- isREML(object)

    famL <- famlink(resp = resp)
    p <- length(coefs <- fixef(object))

    vc <- vcov(object, use.hessian = use.hessian)
    stdError <- sqrt(diag(vc))
    coefs <- cbind("Estimate" = coefs,
                   "Std. Error" = stdError)
    if (p > 0) {
	coefs <- cbind(coefs, (cf3 <- coefs[,1]/coefs[,2]), deparse.level = 0)
	colnames(coefs)[3] <- paste(if(useSc) "t" else "z", "value")
        if (isGLMM(object)) # FIXME: if "t" above, cannot have "z" here
            coefs <- cbind(coefs, "Pr(>|z|)" =
                           2*pnorm(abs(cf3), lower.tail = FALSE))
    }

    llAIC <- llikAIC(object)
    ## FIXME: You can't count on object@re@flist,
    ##	      nor compute VarCorr() unless is(re, "reTrms"):
    varcor <- VarCorr(object)
					# use S3 class for now
    structure(list(methTitle = methTitle(dd),
                   objClass = class(object),
                   devcomp = devC,
                   isLmer = is(resp, "lmerResp"), useScale = useSc,
                   logLik = llAIC[["logLik"]],
                   family = famL$fami, link = famL$link,
		   ngrps = ngrps(object),
		   coefficients = coefs, sigma = sig,
		   vcov = vcov(object, correlation = correlation, sigm = sig),
		   varcor = varcor, # and use formatVC(.) for printing.
		   AICtab = llAIC[["AICtab"]], call = object@call,
                   residuals = residuals(object,"pearson",scaled = TRUE),
		   fitMsgs = .merMod.msgs(object),
                   optinfo = object@optinfo
		   ), class = "summary.merMod")
}

## TODO: refactor?
##' @S3method summary summary.merMod
summary.summary.merMod <- function(object, varcov = TRUE, ...) {
    if(varcov && is.null(object$vcov))
	object$vcov <- vcov(object, correlation = TRUE, sigm = object$sigma)
    object
}

### Plots for the ranef.mer class ----------------------------------------

##' @importFrom lattice dotplot
##' @S3method  dotplot ranef.mer
dotplot.ranef.mer <- function(x, data, main = TRUE, ...)
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
    f <- function(nx, ...) {
        xt <- x[[nx]]
	ss <- stack(xt)
        mtit <- if(main) nx # else NULL
	ss$ind <- factor(as.character(ss$ind), levels = colnames(xt))
	ss$.nn <- rep.int(reorder(factor(rownames(xt)), xt[[1]],
                                  FUN = mean,sort = sort), ncol(xt))
	se <- NULL
	if (!is.null(pv <- attr(xt, "postVar")))
	    se <- unlist(lapply(1:(dim(pv)[1]), function(i) sqrt(pv[i, i, ])))
	dotplot(.nn ~ values | ind, ss, se = se,
		prepanel = prepanel.ci, panel = panel.ci,
		xlab = NULL, main = mtit, ...)
    }
    setNames(lapply(names(x), f, ...), names(x))
}

##' @importFrom graphics plot
##' @S3method plot ranef.mer
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

##' @importFrom lattice qqmath
##' @S3method qqmath ranef.mer
qqmath.ranef.mer <- function(x, data, main = TRUE, ...)
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
    f <- function(nx) {
	xt <- x[[nx]]
        mtit <- if(main) nx # else NULL
	if (!is.null(pv <- attr(xt, "postVar")))
        {
	    d <- dim(pv)
	    se <- vapply(seq_len(d[1]), function(i) sqrt(pv[i, i, ]), numeric(d[3]))
	    nr <- nrow(xt)
	    nc <- ncol(xt)
	    ord <- unlist(lapply(xt, order)) + rep((0:(nc - 1)) * nr, each = nr)
	    rr <- 1:nr
	    ind <- gl(nc, nr, labels = names(xt))
	    xyplot(rep(qnorm((rr - 0.5)/nr), nc) ~ unlist(xt)[ord] | ind[ord],
		   se = se[ord], prepanel = prepanel.ci, panel = panel.ci,
		   scales = list(x = list(relation = "free")),
		   ylab = "Standard normal quantiles",
		   xlab = NULL, main = mtit, ...)
	} else {
	    qqmath(~values|ind, stack(xt),
		   scales = list(y = list(relation = "free")),
		   xlab = "Standard normal quantiles",
		   ylab = NULL, main = mtit, ...)
	}
    }
    sapply(names(x), f, simplify = FALSE)
}

##' @importFrom graphics plot
##' @S3method plot coef.mer
plot.coef.mer <- function(x, y, ...)
{
    ## remove non-varying columns from frames
    reduced <- lapply(x, function(el)
		      el[, !vapply(el, function(cc) all(cc == cc[1L]), NA)])
    plot.ranef.mer(reduced, ...)
}

##' @importFrom lattice dotplot
##' @S3method dotplot coef.mer
dotplot.coef.mer <- function(x, data, ...) {
    mc <- match.call()
    mc[[1]] <- as.name("dotplot.ranef.mer")
    eval(mc)
}

##' @importFrom stats weights
##' @S3method weights merMod
weights.merMod <- function(object, type = c("prior","working"), ...) {
    type <- match.arg(type)
    isPrior <- type == "prior"
    if(!isGLMM(object) && !isPrior)
        stop("working weights only available for GLMMs")
    res <- if(isPrior) object@resp$weights else object@pp$Xwts^2
    ## the working weights available through pp$Xwts should be
    ## equivalent to:
    ##     object@resp$weights*(object@resp$muEta()^2)/object@resp$variance()
    ## however, the unit tests in tests/glmmWeights.R suggest that this
    ## equivalence is approximate.  this may be fine, however, if the
    ## discrepancy is due to another instance of the general problem of
    ## reference class fields not being updated at the optimum, then this
    ## could cause real problems.  see for example:
    ## https://github.com/lme4/lme4/issues/166


    ## FIXME:  what to do about missing values (see stats:::weights.glm)?
    ## FIXME:  add unit tests
    return(res)
}

dim.merMod <- function(x) {
    getME(x, c("n", "p", "q", "p_i", "l_i", "q_i", "k", "m_i", "m"))
}

## Internal utility, only used in optwrap() :
##' @title Get the optimizer function and check it minimally
##' @param optimizer character string ( = function name) *or* function
getOptfun <- function(optimizer) {
    if (((is.character(optimizer) && optimizer == "optimx") ||
         deparse(substitute(optimizer)) == "optimx") &&
        !"package:optimx" %in% search())
        stop(shQuote("optimx")," package must be loaded in order to ",
             "use ",shQuote('optimizer="optimx"'))
    optfun <- if (is.character(optimizer)) {
	tryCatch(get(optimizer), error = function(e) NULL)
    } else optimizer
    if (is.null(optfun)) stop("couldn't find optimizer function ",optimizer)
    if (!is.function(optfun)) stop("non-function specified as optimizer")
    needArgs <- c("fn","par","lower","control")
    if (anyNA(match(needArgs, names(formals(optfun)))))
	stop("optimizer function must use (at least) formal parameters ",
	     paste(sQuote(needArgs), collapse = ", "))
    optfun
}

optwrap <- function(optimizer, fn, par, lower = -Inf, upper = Inf,
                    control = list(), adj = FALSE, calc.derivs = TRUE,
                    use.last.params = FALSE,
                    verbose = 0L)
{
    ## control must be specified if adj==TRUE;
    ##  otherwise this is a fairly simple wrapper
    optfun <- getOptfun(optimizer)
    optName <- if(is.character(optimizer)) optimizer
    else ## "good try":
        deparse(substitute(optimizer))[[1L]]

    lower <- rep(lower, length.out = length(par))
    upper <- rep(upper, length.out = length(par))

    if (adj)
        ## control parameter tweaks: only for second round in nlmer, glmer
        switch(optName,
               "bobyqa" = {
                   if(!is.numeric(control$rhobeg)) control$rhobeg <- 0.0002
                   if(!is.numeric(control$rhoend)) control$rhoend <- 2e-7
               },
               "Nelder_Mead" = {
                   if (is.null(control$xst))  {
                       thetaStep <- 0.1
                       nTheta <- length(environment(fn)$pp$theta)
                       betaSD <- sqrt(diag(environment(fn)$pp$unsc()))
                       control$xst <- 0.2* c(rep.int(thetaStep, nTheta),
                                             pmin(betaSD, 10))
                   }
                   if (is.null(control$xt)) control$xt <- control$xst*5e-4
               })
    switch(optName,
	   "bobyqa" = {
	       if(all(par == 0)) par[] <- 0.001  ## minor kludge
	       if(!is.numeric(control$iprint)) control$iprint <- min(verbose, 3L)
	   },
	   "Nelder_Mead" = control$verbose <- verbose,
	   ## otherwise:
	   if(verbose) warning(gettextf(
	       "'verbose' not yet passed to optimizer '%s'; consider fixing optwrap()",
					optName), domain = NA)
	   )
    arglist <- list(fn = fn, par = par, lower = lower, upper = upper, control = control)
    ## optimx: must pass method in control (?) because 'method' was previously
    ## used in lme4 to specify REML vs ML
    if (optName == "optimx") {
        if (is.null(method <- control$method))
            stop("must specify 'method' explicitly for optimx")
        arglist$control$method <- NULL
        arglist <- c(arglist, list(method = method))
    }
    ## FIXME: test!  effects of multiple warnings??
    ## may not need to catch warnings after all??
    curWarnings <- list()
    opt <- withCallingHandlers(do.call(optfun, arglist),
                               warning = function(w) {
                                   curWarnings <<- append(curWarnings,list(w$message))
                               })
    ## cat("***",unlist(tail(curWarnings,1)))
    ## FIXME: set code to warn on convergence !=0
    ## post-fit tweaking
    if (optName == "bobyqa") {
        opt$convergence <- opt$ierr
    }
    else if (optName == "optimx") {
	opt <- list(par = coef(opt)[1,],
		    fvalues = opt$value[1],
		    method = method,
		    conv = opt$convcode[1],
		    feval = opt$fevals + opt$gevals,
		    message = attr(opt,"details")[,"message"][[1]])
    }
    if (opt$conv != 0) {
        wmsg <- paste("convergence code",opt$conv,"from",optName)
        if (!is.null(opt$msg)) wmsg <- paste0(wmsg,": ",opt$msg)
        warning(wmsg)
        curWarnings <<- append(curWarnings,list(wmsg))
    }
    ## pp_before <- environment(fn)$pp
    ## save(pp_before,file="pp_before.RData")

    if (calc.derivs) {
        if (use.last.params) {
            ## +0 tricks R into doing a deep copy ...
            ## otherwise element of ref class changes!
            ## FIXME:: clunky!!
            orig_pars <- opt$par
            orig_theta <- environment(fn)$pp$theta+0
            orig_pars[seq_along(orig_theta)] <- orig_theta
        }
        if (verbose > 10) cat("computing derivatives\n")
        derivs <- deriv12(fn,opt$par,fx = opt$value)
        if (use.last.params) {
            ## run one more evaluation of the function at the optimized
            ##  value, to reset the internal/environment variables in devfun ...
            fn(orig_pars)
        }
    } else derivs <- NULL

    if (!use.last.params) {
        ## run one more evaluation of the function at the optimized
        ##  value, to reset the internal/environment variables in devfun ...
        fn(opt$par)
    }
    structure(opt, ## store all auxiliary information
              optimizer = optimizer,
              control   = control,
              warnings  = curWarnings,
              derivs    = derivs)
}

as.data.frame.VarCorr.merMod <- function(x,row.names = NULL,
                                         optional = FALSE,
                                         order = c("cov.last", "lower.tri"),
                                         ...)  {
    order <- match.arg(order)
    tmpf <- function(v,grp) {
        vcov <- c(diag(v), v[lt.v <- lower.tri(v, diag = FALSE)])
        sdcor <- c(attr(v,"stddev"),
                   attr(v,"correlation")[lt.v])
        nm <- rownames(v)
        n <- nrow(v)
        dd <- data.frame(grp = grp,
                         var1 = nm[c(seq(n), col(v)[lt.v])],
                         var2 = c(rep(NA,n), nm[row(v)[lt.v]]),
                         vcov,
                         sdcor,
                         stringsAsFactors = FALSE)
        if (order=="lower.tri") {
            ## reorder *back* to lower.tri order
            m <- matrix(NA,n,n)
            diag(m) <- seq(n)
            m[lower.tri(m)] <- (n+1):(n*(n+1)/2)
            dd <- dd[m[lower.tri(m, diag=TRUE)],]
        }
        dd
    }
    r <- do.call(rbind,
                 c(mapply(tmpf, x,names(x), SIMPLIFY = FALSE),
                   deparse.level = 0))
    if (attr(x,"useSc")) {
        ss <- attr(x,"sc")
        r <- rbind(r,data.frame(grp = "Residual",var1 = NA,var2 = NA,
                                vcov = ss^2,
                                sdcor = ss),
		   deparse.level = 0)
    }
    rownames(r) <- NULL
    r
}


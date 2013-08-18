##' @rdname modular
##' @name modular
##' @title Modular functions for mixed model fits
##' @aliases modular
##' @inheritParams lmer
##' @param \dots other potential arguments.
##' @details These functions make up the internal components of a(n) [gn]lmer fit.
##' \itemize{
##' \item \code{[g]lFormula} takes the arguments that would normally be passed to \code{[g]lmer},
##' checking for errors and processing the formula and data input to create
##' \item \code{mk(Gl|L)merDevfun} takes the output of the previous step (minus the \code{formula}
##' component) and creates a deviance function
##' \item \code{optimize(Gl|L)mer} takes a deviance function and optimizes over \code{theta} (or over \code{theta} and \code{beta}, if \code{stage} is set to 2 for \code{optimizeGlmer}
##' \item \code{updateGlmerDevfun} takes the first stage of a GLMM optimization (with \code{nAGQ=0},
##' optimizing over \code{theta} only) and produces a second-stage deviance function
##' \item \code{\link{mkMerMod}} takes the \emph{environment} of a deviance function, the results of
##' an optimization, a list of random-effect terms, a model frame, and a model all and produces a
##' \code{[g]lmerMod} object
##' }
##' @examples
##'
##' ### Fitting a linear mixed model in 4 modularized steps
##'
##' ## 1.  Parse the data and formula:
##' lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
##' names(lmod)
##' ## 2.  Create the deviance function to be optimized:
##' (devfun <- do.call(mkLmerDevfun, lmod))
##' ls(environment(devfun)) # the environment of devfun contains objects required for its evaluation
##' ## 3.  Optimize the deviance function:
##' opt <- optimizeLmer(devfun)
##' opt[1:3]
##' ## 4.  Package up the results:
##' mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
##'
##'
##' ### Same model in one line
##' lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##'
##'
##' ### Fitting a generalized linear mixed model in six modularized steps
##'
##' ## 1.  Parse the data and formula:
##' glmod <- glFormula(cbind(incidence, size - incidence) ~ period + (1 | herd),
##' data = cbpp, family = binomial)
##' names(glmod)
##' ## 2.  Create the deviance function for optimizing over theta:
##' (devfun <- do.call(mkGlmerDevfun, glmod))
##' ls(environment(devfun)) # the environment of devfun contains lots of info
##' ## 3.  Optimize over theta using a rough approximation (i.e. nAGQ = 0):
##' (opt <- optimizeGlmer(devfun))
##' ## 4.  Update the deviance function for optimizing over theta and beta:
##' (devfun <- updateGlmerDevfun(devfun, glmod$reTrms))
##' ## 5.  Optimize over theta and beta:
##' opt <- optimizeGlmer(devfun, stage=2)
##' opt[1:3]
##' ## 6.  Package up the results:
##' mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr)
##'
##'
##' ### Same model in one line
##' glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##' data = cbpp, family = binomial)
##'
NULL

### Small utilities to be used in lFormula() and glFormula()

doCheck <- function(x) {
    is.character(x) && !any(x == "ignore")
}

checkZrank <- function(Zt, n, ctrl, nonSmall = 1e6, allow.n=FALSE)
{
    stopifnot(is.list(ctrl), is.numeric(n), is.numeric(nonSmall))
    cstr <- "check.nobs.vs.rankZ"
    if (doCheck(cc <- ctrl[[cstr]])) { ## not NULL or "ignore"
	d <- dim(Zt)
	doTr <- d[1L] < d[2L] # Zt is "wide" => qr needs transpose(Zt)
	if(!(grepl("Small",cc) && prod(d) > nonSmall)) {
            rankZ <- rankMatrix(if(doTr) t(Zt) else Zt, method="qr",
                                sval = numeric(min(d)))
            if (allow.n) {
                unident <- n<rankZ
                cmp <- "<"
            } else {
                unident <- n<=rankZ
                cmp <- "<="
            }
            ## OR: cmp <- if (allow.n) "<" else "<="
            ##  if (do.call(cmp,list(n,rankZ)) ...
            if (unident) {
                wstr <- paste("number of observations",cmp,"rank(Z); variance-covariance matrix will be unidentifiable")
                switch(cc,
                       "warningSmall" =, "warning" = warning(wstr),
                       "stopSmall" =, "stop" = stop(wstr),
                       stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
            }
        }
    }
}

checkNlevels <- function(flist, n, ctrl, allow.n=FALSE)
{
    stopifnot(is.list(ctrl), is.numeric(n))
    nlevelVec <- unlist(lapply(flist, function(x) nlevels(droplevels(x)) ))
    ## Part 1 ----------------
    cstr <- "check.nlev.gtr.1"
    if (doCheck(cc <- ctrl[[cstr]]) && any(nlevelVec < 2)) {
	wstr <- "grouping factors must have > 1 sampled level"
	switch(cc,
	       "warning" = warning(wstr),
	       "stop" = stop(wstr),
	       stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
    }
    ## Part 2 ----------------
    cstr <- "check.nobs.vs.nlev"
    if (doCheck(cc <- ctrl[[cstr]])) {
        if (any(if(allow.n) nlevelVec > n else nlevelVec >= n))
            stop(gettextf(
                "number of levels of each grouping factor must be %s number of observations",
                if(allow.n) "<=" else "<"), domain=NA)
    }

    ## Part 3 ----------------
    cstr <- "check.nlev.gtreq.5"
    if (doCheck(cc <- ctrl[[cstr]]) && any(nlevelVec < 5)) {
	wstr <- "grouping factors with < 5 sampled levels may give unreliable estimates"
	switch(cc,
	       "warning" = warning(wstr),
	       "stop" = stop(wstr),
	       stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
    }
}



##' @rdname modular
##' @param control a list giving (for \code{[g]lFormula}) all options (see \code{\link{lmerControl}} for running the model;
##' (for \code{mkLmerDevfun,mkGlmerDevfun}) options for inner optimization step;
##' (for \code{optimizeLmer} and \code{optimize[Glmer}) control parameters for nonlinear optimizer (typically inherited from the \dots argument to \code{lmerControl})
##' @return \bold{lFormula, glFormula}: A list containing components,
##' \item{fr}{model frame}
##' \item{X}{fixed-effect design matrix}
##' \item{reTrms}{list containing information on random effects structure: result of \code{\link{mkReTrms}}}
##' \item{REML}{(lFormula only): logical flag: use restricted maximum likelihood? (Copy of argument.)}
##' @importFrom Matrix rankMatrix
##' @export
lFormula <- function(formula, data=NULL, REML = TRUE,
                     subset, weights, na.action, offset, contrasts = NULL,
                     control=lmerControl(), ...)
{
    control <- control$checkControl ## this is all we really need
    mf <- mc <- match.call()

    ignoreArgs <- c("start","verbose","devFunOnly","control")
    l... <- list(...)
    l... <- l...[!names(l...) %in% ignoreArgs]
    do.call("checkArgs",c(list("lmer"),l...))
    if (!is.null(list(...)[["family"]])) {
        ## lmer(...,family=...); warning issued within checkArgs
        mc[[1]] <- quote(lme4::glFormula)
        if (missing(control)) mc[["control"]] <- glmerControl()
        return(eval(mc, parent.frame()))
    }

    denv <- checkFormulaData(formula,data)
    mc$formula <- formula <- as.formula(formula,env=denv) ## substitute evaluated call
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substitute "|" by "+"
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    ## store full, original formula & offset
    attr(fr,"formula") <- formula
    attr(fr,"offset") <- mf$offset
    n <- nrow(fr)
    ## random effects and terms modules
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    checkNlevels(reTrms$ flist, n=n, control)
    checkZrank	(reTrms$ Zt,	n=n, control, nonSmall = 1e6)

    ## fixed-effects model matrix X - remove random effects from formula:
    fixedform <- formula
    fixedform[[3]] <- if(is.null(nb <- nobars(fixedform[[3]]))) 1 else nb
    mf$formula <- fixedform
    ## re-evaluate model frame to extract predvars component {FIXME? glFormula() does not..}
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"),"predvars.fixed") <- attr(attr(fixedfr,"terms"),"predvars")
    X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    p <- ncol(X)
    if ((rankX <- rankMatrix(X)) < p)
	stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))
    list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula)
}

## utility f'n for checking starting values
getStart <- function(start,lower,pred,returnVal=c("theta","all")) {
    returnVal <- match.arg(returnVal)
    ## default values
    theta <- pred$theta
    fixef <- pred$delb
    if (!is.null(start)) {
        if (is.numeric(start)) {
            theta <- start
        } else {
            if (!is.list(start)) stop("start must be a list or a numeric vector")
            if (!all(sapply(start,is.numeric))) stop("all elements of start must be numeric")
            if (length((badComp <- setdiff(names(start),c("theta","fixef"))))>0) {
                stop("incorrect components in start list: ",badComp)
            }
            if (!is.null(start$theta)) theta <- start$theta
            if (!is.null(start$fixef)) fixef <- start$fixef
        }
    }
    if (length(theta)!=length(pred$theta))
        stop("incorrect number of theta components (!=",length(pred$theta),")")
    if (length(fixef)!=length(pred$delb))
        stop("incorrect number of fixef components (!=",length(pred$delb),")")
    if (returnVal=="theta") theta else c(theta,fixef)
}

## update start
## should refactor this to
##  turn numeric start into start=list(theta=start) immediately ... ??
updateStart <- function(start,theta) {
    if (is.null(start)) return(NULL)
    if (is.numeric(start)) {
        start <- theta
    } else if (!is.null(start$theta)) start$theta <- theta
    start
}

##' @rdname modular
##' @param fr A model frame containing the variables needed to create an
##'   \code{\link{lmerResp}} or \code{\link{glmResp}} instance
##' @param X fixed-effects design matrix
##' @param reTrms information on random effects structure (see \code{\link{mkReTrms}})
##' @param REML (logical) fit restricted maximum likelihood model?
##' @param start starting values
##' @param verbose print output?
##' @return \bold{mkLmerDevfun, mkGlmerDevfun}: A function to calculate deviance
##' (or restricted deviance) as a function of the theta (random-effect) parameters
##' (for GlmerDevfun, of beta (fixed-effect) parameters as well).  These deviance
##' functions have an environment containing objects required for their evaluation.
##' CAUTION: The output object of \code{mk(Gl|L)merDevfun} is an \code{\link{environment}}
##' containing reference class objects (see \code{\link{ReferenceClasses}}, \code{\link{merPredD-class}},
##' \code{\link{lmResp-class}}), which behave in ways that may surprise many users. For example, if the
##' output of \code{mk(Gl|L)merDevfun} is naively copied, then modifications to the original will
##' also appear in the copy (and vice versa). To avoid this behavior one must make a deep copy
##' (see \code{\link{ReferenceClasses}} for details).
##' \cr
##' \cr
##' @export
mkLmerDevfun <- function(fr, X, reTrms, REML = TRUE, start = NULL, verbose=0, control=lmerControl(), ...) {

    ## FIXME: make sure verbose gets handled properly
    #if (missing(fr)) {
    ## reconstitute frame
    #}
    ## pull necessary arguments for making the model frame out of ...
    p <- ncol(X) # maybe also do rank check on X here??
    rho <- new.env(parent=parent.env(environment()))
    rho$pp <- do.call(merPredD$new, c(reTrms[c("Zt","theta","Lambdat","Lind")],
                                      n=nrow(X), list(X=X)))
    REMLpass <- if(REML) p else 0L
    if(missing(fr)) rho$resp <- mkRespMod(REML = REMLpass, ...)
    else rho$resp <- mkRespMod(fr, REML = REMLpass)
    ## note:  REML does double duty as rank of X and a flag for using REML
    ## maybe this should be mentioned in the help file for mkRespMod??
    ## currently that help file says REML is logical
    devfun <- mkdevfun(rho, 0L, verbose, control)
    theta <- getStart(start,reTrms$lower,rho$pp)
    devfun(rho$pp$theta) # one evaluation to ensure all values are set
    rho$lower <- reTrms$lower # SCW:  in order to be more consistent with mkLmerDevfun
    return(devfun) # this should pass the rho environment implicitly
}


##' @rdname modular
##' @inheritParams lmer
##' @inheritParams lmerControl
##' @param devfun a deviance function, as generated by \code{\link{mkLmerDevfun}}
##' @return \bold{optimizeLmer}: Results of an optimization.
##' \cr
##' \cr
##' @export
optimizeLmer <- function(devfun,
                         optimizer="Nelder_Mead",
                         restart_edge=FALSE,
                         start = NULL,
                         verbose = 0L,
                         control = list()) {
    verbose <- as.integer(verbose)
    rho <- environment(devfun)
    opt <- optwrap(optimizer,
                   devfun,
                   getStart(start,rho$lower,rho$pp),
                   lower=rho$lower,
                   control=control,
                   adj=FALSE, verbose=verbose)

    if (restart_edge) {
        ## FIXME: should we be looking at rho$pp$theta or opt$par
        ##  at this point???  in koller example (for getData(13)) we have
        ##   rho$pp$theta=0, opt$par=0.08
        if (length(bvals <- which(rho$pp$theta==rho$lower))>0) {
            ## *don't* use numDeriv -- cruder but fewer dependencies, no worries
            ##  about keeping to the interior of the allowed space
            theta0 <- new("numeric",rho$pp$theta) ## 'deep' copy ...
            d0 <- devfun(theta0)
            btol <- 1e-5  ## FIXME: make user-settable?
            bgrad <- sapply(bvals,
                            function(i) {
                                bndval <- rho$lower[i]
                                theta <- theta0
                                theta[i] <- bndval+btol
                                (devfun(theta)-d0)/btol
                            })
            ## what do I need to do to reset rho$pp$theta to original value???
            devfun(theta0) ## reset rho$pp$theta after tests
            ## FIXME: allow user to specify ALWAYS restart if on boundary?
            if (any(bgrad<0)) {
                if (verbose) message("some theta parameters on the boundary, restarting")
                opt <- optwrap(optimizer,
                               devfun,
                               opt$par,
                               lower=rho$lower, control=control,
                               adj=FALSE, verbose=verbose)
            }
        }
    }
    return(opt)
}

## TODO: remove any arguments that aren't actually used by glFormula (same for lFormula)
## TODO(?): lFormula() and glFormula()  are very similar: merge or use common baseFun()
##' @rdname modular
##' @inheritParams glmer
##' @export
glFormula <- function(formula, data=NULL, family = gaussian,
                      subset, weights, na.action, offset,
                      contrasts = NULL, mustart, etastart,
                      control=glmerControl(), ...) {
    ## FIXME: does start= do anything? test & fix

    control <- control$checkControl ## this is all we really need
    mf <- mc <- match.call()
    ## extract family, call lmer for gaussian
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if( is.function(family)) family <- family()
    if (isTRUE(all.equal(family, gaussian()))) {
        mc[[1]] <- quote(lme4::lFormula)
        mc["family"] <- NULL            # to avoid an infinite loop
        return(eval(mc, parent.frame()))
    }
    if (family$family %in% c("quasibinomial", "quasipoisson", "quasi"))
        stop('"quasi" families cannot be used in glmer')

    ignoreArgs <- c("start","verbose","devFunOnly","optimizer", "control", "nAGQ")
    l... <- list(...)
    l... <- l...[!names(l...) %in% ignoreArgs]
    do.call("checkArgs",c(list("glmer"),l...))

    denv <- checkFormulaData(formula,data)
    mc$formula <- formula <- as.formula(formula,env=denv)    ## substitute evaluated version

    m <- match(c("data", "subset", "weights", "na.action", "offset",
                 "mustart", "etastart"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substitute "|" by "+"
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    ## store full, original formula & offset
    attr(fr,"formula") <- formula
    attr(fr,"offset") <- mf$offset
    n <- nrow(fr)
    ## random effects and terms modules
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    ## TODO: allow.n = !useSc {see FIXME below}
    checkNlevels(reTrms$ flist, n=n, control, allow.n=TRUE)
    checkZrank	(reTrms$ Zt,	n=n, control, nonSmall = 1e6, allow.n=TRUE)

    ## FIXME: adjust test for families with estimated scale parameter:
    ##   useSc is not defined yet/not defined properly?
    ##  if (useSc && maxlevels == n)
    ##          stop("number of levels of each grouping factor must be",
    ##                "greater than number of obs")

    ## fixed-effects model matrix X - remove random parts from formula:
    form <- formula
    form[[3]] <- if(is.null(nb <- nobars(form[[3]]))) 1 else nb
    X <- model.matrix(form, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    p <- ncol(X)
    if ((rankX <- rankMatrix(X)) < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))

    list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula)
}

##' @rdname modular
##' @export
mkGlmerDevfun <- function(fr, X, reTrms, family, nAGQ = 1L, verbose = 0L,
                          control=glmerControl(), ...){
    stopifnot(length(nAGQ <- as.integer(nAGQ)) == 1L,
            nAGQ >= 0L,
            nAGQ <= 25L)
    verbose <- as.integer(verbose)
    rho             <- as.environment(list(verbose=verbose, tolPwrss=control$tolPwrss,
                                           compDev=control$compDev))
    parent.env(rho) <- parent.frame()
    rho$pp          <- do.call(merPredD$new,
                               c(reTrms[c("Zt","theta","Lambdat","Lind")],
                                 n=nrow(X), list(X=X)))
    if (missing(fr)) rho$resp <- mkRespMod(family=family, ...)
    else rho$resp             <- mkRespMod(fr, family=family)
    if (length(unique(rho$resp$y)) < 2L)
        stop("Response is constant - cannot fit the model")
    rho$verbose     <- as.integer(verbose)
    ## initialize (from mustart)
    .Call(glmerLaplace, rho$pp$ptr(), rho$resp$ptr(), 0L, control$tolPwrss, verbose)
    rho$lp0         <- rho$pp$linPred(1) # each pwrss opt begins at this eta
    rho$pwrssUpdate <- glmerPwrssUpdate
    rho$lower       <- reTrms$lower     # not needed in rho?
    devfun <- mkdevfun(rho, 0L, verbose, control)
                                        #if (devFunOnly && !nAGQ) return(devfun)
    return(devfun) # this should pass the rho environment implicitly
}



##' @rdname modular
##' @param nAGQ number of Gauss-Hermite quadrature points
##' @param stage optimization stage (1: nAGQ=0, optimize over theta only; 2: nAGQ possibly >0, optimize over theta and beta)
##' @export
optimizeGlmer <- function(devfun,
                          optimizer="bobyqa",
                          restart_edge=FALSE,
                          verbose = 0L,
                          control = list(),
                          nAGQ = 1L,
                          stage = 1,
                          start = NULL) {
    ## FIXME: do we need nAGQ here?? or can we clean up?
    verbose <- as.integer(verbose)
    rho <- environment(devfun)
    if (stage==1) {
        start <- getStart(start, lower=rho$lower, pred=rho$pp, "theta")
        adj <- FALSE
    } else { ## stage == 2
        start <- getStart(start, lower=rho$lower, pred=rho$pp, returnVal="all")
        adj <- TRUE
        if (missing(optimizer)) optimizer <- "Nelder_Mead"  ## BMB: too clever?
    }
    opt <- optwrap(optimizer, devfun, start, rho$lower,
                   control=control, adj=adj, verbose=verbose)
    if (stage==1) {
        rho$control <- attr(opt,"control")
        rho$nAGQ <- nAGQ
    } else {  ## stage == 2
        rho$resp$setOffset(rho$baseOffset)
    }
    ## FIXME: implement this ...
    if (restart_edge) stop("restart_edge not implemented for optimizeGlmer yet")
    return(opt)
}

## only do this function if nAGQ > 0L
##' @rdname modular
##' @export
updateGlmerDevfun <- function(devfun, reTrms, nAGQ = 1L){
    rho <- environment(devfun)
    rho$nAGQ       <- nAGQ
    rho$lower      <- c(rho$lower, rep.int(-Inf, length(rho$pp$beta0)))
    rho$lp0        <- rho$pp$linPred(1)
    rho$dpars      <- seq_along(rho$pp$theta)
    rho$baseOffset <- rho$resp$offset + 0 # forcing a copy (!)
    rho$GQmat      <- GHrule(nAGQ)
    rho$fac        <- reTrms$flist[[1]]
    if (nAGQ > 1L) {
        if (length(reTrms$flist) != 1L || length(reTrms$cnms[[1]]) != 1L)
            stop("nAGQ > 1 is only available for models with a single, scalar random-effects term")
    }
    devfun <- mkdevfun(rho, nAGQ)  # does this attach rho to devfun??
    return(devfun)
}

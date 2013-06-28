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

doCheck <- function(x) {
    !is.null(x) && x!="ignore"
}

##' @rdname modular
##' @param control a list giving (for \code{[g]lFormula}) all options (see \code{\link{lmerControl}} for running the model;
##' (for \code{mkLmerDevfun,mkGlmerDevfun}) options for inner optimization step;
##' (for \code{optimizeLmer} and \code{optimize[Glmer}) control parameters for nonlinear optimizer (typically inherited from the \dots argument to \code{lmerControl})
##' @return \bold{lFormula, glFormula}: A list containing components,
##' \item{fr}{model frame}
##' \item{X}{fixed-effect design matrix}
##' \item{reTrms}{list containing information on random effects structure: result of \code{\link{mkMerMod}}}
##' \item{REML}{(lFormula only): logical flag: use restricted maximum likelihood? (Copy of argument.)}
##' @export
lFormula <- function(formula, data=NULL, REML = TRUE, 
                     subset, weights, na.action, offset, contrasts = NULL,
                     control=lmerControl(), ...)
{

    mf <- mc <- match.call()

    control <- control$checkControl ## this is all we really need
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
    fr.form <- subbars(formula) # substituted "|" by "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    ## random effects and terms modules
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    nlevelVec <- unlist(lapply(reTrms$flist, function(x) nlevels(droplevels(x)) ))
    cstr <- "check.numlev.gtr.1"
    if (doCheck(cc <- control[[cstr]]) && any(nlevelVec<2)) {
        wstr <- "grouping factors must have > 1 sampled level"
        switch(cc,warning=warning(wstr),stop=stop(wstr),
               stop(paste0("unknown check level for '",cstr,"'")))
    }
    if (any(nlevelVec >= nrow(fr)))
        stop("number of levels of each grouping factor must be ",
             "< number of observations")
    cstr <- "check.numlev.gtreq.5"
    if (doCheck(cc <- control[[cstr]]) && any(nlevelVec<5)) {
        wstr <- "grouping factors with < 5 sampled levels may give unreliable estimates"
        switch(cc,warning=warning(wstr),stop=stop(wstr),
               stop(paste0("unknown check level for '",cstr,"'")))
    }
    cstr <- "check.numobs.vs.rankZ"
    if (doCheck(cc <- control[[cstr]]) &&                   ## not NULL or "ignore"
        !(grepl("Small",cc) && prod(dim(reTrms$Zt))>1e6)    ## not "*Small" and large Z mat
        && (nrow(fr) <= (rankZ <- Matrix::rankMatrix(reTrms$Zt))))
    {
        wstr <- "number of observations <= rank(Z) ; variance-covariance matrix is likely to be unidentifiable"
        switch(cc,warningSmall=,warning=warning(wstr),stopSmall=,stop=stop(wstr),
               stop(paste0("unknown check level for '",cstr,"'")))
    }
    ## fixed-effects model matrix X - remove random effects from formula:
    fixedform <- formula
    fixedform[[3]] <- if(is.null(nb <- nobars(fixedform[[3]]))) 1 else nb
    mf$formula <- fixedform
    ## re-evaluate model frame to extract predvars component
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"),"predvars.fixed") <- attr(attr(fixedfr,"terms"),"predvars")
    X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    p <- ncol(X)
    if ((rankX <- Matrix::rankMatrix(X)) < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))
    return(list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula))
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
    ## FIXME: should we apply start elsewhere? what about reTrms$theta?
    if (!is.null(start)) {
        rho$pp$setTheta(unlist(start))
    }
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
                         verbose = 0L,
                         control = list()) {
    verbose <- as.integer(verbose)
    rho <- environment(devfun)
    opt <- optwrap(optimizer,
                   devfun, rho$pp$theta, lower=rho$lower, control=control,
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
            ## FIXME: opt$fval is wrong
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
                               devfun, opt$par,
                               lower=rho$lower, control=control$optControl,
                               adj=FALSE, verbose=verbose)
            }
        }
    }
    return(opt)
}

## TODO:  remove any arguments that aren't actually used by glFormula (same for lFormula)
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
                                        # extract family, call lmer for gaussian
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
    fr.form <- subbars(formula) # substitute "|" for "+" -
    environment(fr.form) <- environment(formula)
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
                                        # random-effects module
    reTrms <- mkReTrms(findbars(formula[[3]]), fr)
    nlevelVec <- unlist(lapply(reTrms$flist, function(x) nlevels(droplevels(x)) ))
    if (any(nlevelVec==1)) stop("grouping factors must have > 1 sampled level")
    if ((maxlevels <- max(nlevelVec)) > nrow(fr))
        stop("number of levels of each grouping factor must be",
             ">= number of obs")
    ## FIXME: duplicated code between lFormula, glFormula
    if (any(nlevelVec<5))  warning("grouping factors with < 5 sampled levels may give unreliable estimates")
    cstr <- "check.numobs.vs.rankZ"
    ## FIXME: adjust test for families with estimated scale parameter:
    ##   useSc is not defined yet/not defined properly?
    if (doCheck(cc <- control[[cstr]]) &&                   ## not NULL or "ignore"
        !(grepl(cc,"Small") && prod(dim(reTrms$Zt))>1e6)    ## not "*Small" and large Z mat
        && (nrow(fr) < (rankZ <- Matrix::rankMatrix(reTrms$Zt))))    ## test
    {
        wstr <- "number of observations<rank(Z); variance-covariance matrix will be unidentifiable"
        switch(cc,warningSmall=,warning=warning(wstr),stopSmall=,stop=stop(wstr),
               stop(paste0("unknown check level for '",cstr,"'")))
    }
    ##  if (useSc && maxlevels == nrow(fr))
    ##          stop("number of levels of each grouping factor must be",
    ##                "greater than number of obs")
    
    ## fixed-effects model matrix X - remove random parts from formula:
    form <- formula
    form[[3]] <- if(is.null(nb <- nobars(form[[3]]))) 1 else nb
    X <- model.matrix(form, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    p <- ncol(X)
    
    if ((rankX <- Matrix::rankMatrix(X)) < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))
    
    out <- list(fr = fr, X = X, reTrms = reTrms, family = family)
    attr(out, "formula") <- formula
    return(out)
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
                          stage = 1) {
    ## FIXME: do we need nAGQ here?? or can we clean up?
    verbose <- as.integer(verbose)
    rho <- environment(devfun)
    if (stage==1) {
        start <- rho$pp$theta
        adj <- FALSE
    } else { ## stage == 2
        start <- c(rho$pp$theta, rho$pp$delb)
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
    ## FIXME: should we worry about restart_edge?
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

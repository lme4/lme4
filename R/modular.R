#### --> ../man/modular.Rd
####     ==================

### Small utilities to be used in lFormula() and glFormula()

doCheck <- function(x) {
    is.character(x) && !any(x == "ignore")
}

RHSForm <- function(formula) {
    formula[[length(formula)]]
}

`RHSForm<-` <- function(formula,value) {
    formula[[length(formula)]] <- value
    formula
}

##' Original formula, minus response ( = '~ <RHS>') :
noLHSform <- function(formula) {
    if (length(formula)==2) formula else formula[-2]
}

##' @param cstr name of control being set
##' @param val value of control being set
checkCtrlLevels <- function(cstr, val, smallOK=FALSE) {
    bvals <- c("stop","warning","ignore")
    if (smallOK) bvals <- outer(bvals,c("","Small"),paste0)
    if (!is.null(val) && !val %in% bvals)
        stop("invalid control level ",sQuote(val)," in ",cstr,": valid options are {",
             paste(sapply(bvals,sQuote),collapse=","),"}")
    invisible(NULL)
}

## general identifiability checker, used both in checkZdim and checkZrank
wmsg <- function(n,cmp.val,allow.n,msg1="",msg2="",msg3="") {
    if (allow.n) {
        unident <- n < cmp.val
        cmp <- "<"
        rstr <- ""
    } else {
        unident <- n <= cmp.val
        cmp <- "<="
        rstr <- " and the residual variance (or scale parameter)"
    }
    ## %s without spaces intentional (don't want an extra space if the
    ## message component is empty)
    wstr <- sprintf("%s (=%d) %s %s (=%d)%s; the random-effects parameters%s are probably unidentifiable",msg1,n,cmp,msg2,cmp.val,msg3,rstr)
    list(unident=unident,wstr=wstr)
}

checkZdims <- function(Ztlist, n, ctrl, allow.n=FALSE) {
    ## Ztlist: list of Zt matrices - one for each r.e. term
    ## n: no. observations
    ## allow.n: allow as many random-effects as there are observations
    ## for each term?
    ##
    ## For each r.e. term, test if Z has more columns than rows to detect
    ## unidentifiability:
    stopifnot(is.list(Ztlist), is.numeric(n))
    cstr <- "check.nobs.vs.nRE"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    term.names <- names(Ztlist)
    rows <- vapply(Ztlist, nrow, numeric(1L))
    cols <- vapply(Ztlist, ncol, numeric(1L))
    stopifnot(all(cols == n))
    if (doCheck(cc)) {
        for(i in seq_along(Ztlist)) {
            ww <- wmsg(cols[i],rows[i],allow.n,"number of observations",
                       "number of random effects",
                       sprintf(" for term (%s)",term.names[i]))
            if(ww$unident) {
            switch(cc,
                   "warning" = warning(ww$wstr,call.=FALSE),
                   "stop" = stop(ww$wstr,call.=FALSE),
                   stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
            }
        }
    }
}


checkZrank <- function(Zt, n, ctrl, nonSmall = 1e6, allow.n=FALSE)
{
    stopifnot(is.list(ctrl), is.numeric(n), is.numeric(nonSmall))
    cstr <- "check.nobs.vs.rankZ"
    if (doCheck(cc <- ctrl[[cstr]])) { ## not NULL or "ignore"
        checkCtrlLevels(cstr, cc, smallOK=TRUE)
	d <- dim(Zt)
	doTr <- d[1L] < d[2L] # Zt is "wide" => qr needs transpose(Zt)
	if(!(grepl("Small",cc) && prod(d) > nonSmall)) {
            rankZ <- rankMatrix(if(doTr) t(Zt) else Zt, method="qr",
                                sval = numeric(min(d)))
            ww <- wmsg(n,rankZ,allow.n,"number of observations","rank(Z)")
            if (ww$unident) {
                switch(cc,
                       "warningSmall" =, "warning" = warning(ww$wstr,call.=FALSE),
                       "stopSmall" =, "stop" = stop(ww$wstr,call.=FALSE),
                       stop(gettextf("unknown check level for '%s'", cstr),
                            domain=NA))
            }
        }
    }
}

## check scale of non-dummy columns of X, both
## against each other and against 1 (implicit scale of theta parameters)?
## (shouldn't matter for lmer models?)
## TODO: check for badly centred models?
## TODO: check scale of Z columns?
## What should the rules be?  try to find problematic columns
##   and rescale them?  scale+center?  Just scale or scale+center
##   all numeric columns?
## 
checkScaleX <- function(X,  kind="warning", tol=1e3) {
    cstr <- "check.scaleX"
    kinds <- eval(formals(lmerControl)[["check.scaleX"]])
    if (!kind %in% kinds) stop(sprintf("unknown check-scale option: %s",kind))
    if (is.null(kind) || kind == "ignore") return(X)
    ## else :
    cont.cols <- apply(X,2,function(z) !all(z %in% c(0,1)))
    col.sd <- apply(X[,cont.cols, drop=FALSE],2,sd)
    sdcomp <- outer(col.sd,col.sd,"/")
    logcomp <- abs(log(sdcomp[lower.tri(sdcomp)]))
    logsd <- abs(log(col.sd))
    wmsg <- "Some predictor variables are on very different scales: consider rescaling"
    if (any(c(logcomp,logsd) > log(tol))) {
        if (kind %in% c("warning","stop")) {
            switch(kind, "warning" = warning(wmsg, call.=FALSE),
                   "stop" = stop(wmsg, call.=FALSE))
        } else {
            ## mimic scale() because we don't want to make a copy in
            ##  order to retrieve the center/scale
            X[,cont.cols] <- sweep(X[,cont.cols,drop=FALSE],2,col.sd,"/")
            attr(X,"scaled:scale") <- setNames(col.sd,colnames(X)[cont.cols])
            wmsg <- "Some predictor variables on very different scales: auto-rescaled (results NOT adjusted)"
            if (kind=="warn+rescale") warning(wmsg,call.=FALSE)
        }
    }
    X
}

checkNlevels <- function(flist, n, ctrl, allow.n=FALSE)
{
    stopifnot(is.list(ctrl), is.numeric(n))
    nlevelVec <- unlist(lapply(flist, function(x) nlevels(droplevels(x)) ))
    ## Part 1 ----------------
    cstr <- "check.nlev.gtr.1"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    if (doCheck(cc) && any(nlevelVec < 2)) {
	wstr <- "grouping factors must have > 1 sampled level"
	switch(cc,
	       "warning" = warning(wstr,call.=FALSE),
	       "stop" = stop(wstr,call.=FALSE),
	       stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
    }
    ## Part 2 ----------------
    cstr <- "check.nobs.vs.nlev"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    if (doCheck(cc)) {
        if (any(if(allow.n) nlevelVec > n else nlevelVec >= n))
            stop(gettextf(
                "number of levels of each grouping factor must be %s number of observations",
                if(allow.n) "<=" else "<"), domain=NA)
    }

    ## Part 3 ----------------
    cstr <- "check.nlev.gtreq.5"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    if (doCheck(cc) && any(nlevelVec < 5)) {
	wstr <- "grouping factors with < 5 sampled levels may give unreliable estimates"
	switch(cc,
	       "warning" = warning(wstr,call.=FALSE),
	       "stop" = stop(wstr,call.=FALSE),
	       stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
    }
}

##' Coefficients (columns) are dropped from a design matrix to
##' ensure that it has full rank.
##'
##'  Redundant columns of the design matrix are identified with the
##'  LINPACK implementation of the \code{\link{qr}} decomposition and
##'  removed. The returned design matrix will have \code{qr(X)$rank}
##'  columns.
##'
##' (Note: I have lifted this function from the ordinal (soon Rufus)
##' package and modified it slightly./rhbc)
##'
##' @title Ensure full rank design matrix
##' @param X a design matrix, e.g., the result of
##' \code{\link{model.matrix}} possibly of less than full column rank,
##' i.e., with redundant parameters.
##' @param silent should a message not be issued if X is column rank
##' deficient?
##' @return The design matrix \code{X} without redundant columns.
##' @seealso \code{\link{qr}} and \code{\link{lm}}
##' @author Rune Haubo Bojesen Christensen (drop.coef()); Martin Maechler
chkRank.drop.cols <- function(X, kind, tol = 1e-7, method = "qr.R") {
    ## Test and match arguments:
    stopifnot(is.matrix(X))
    kinds <- eval(formals(lmerControl)[["check.rankX"]])
    ## c("message+drop.cols", "ignore",
    ##   "silent.drop.cols", "warn+drop.cols", "stop.deficient"),

    if(kind == "ignore") return(X)
    ## else :
    p <- ncol(X)
    if (kind == "stop.deficient") {
        if ((rX <- rankMatrix(X, tol=tol, method=method)) < p)
            stop(gettextf("the fixed-effects model matrix is column rank deficient (rank(X) = %d < %d = p); the fixed effects will be jointly unidentifiable",
                          rX, p),
                 call. = FALSE)
    } else {
        ## kind is one of "message+drop.cols", "silent.drop.cols", "warn+drop.cols"
        ## --> consider to drop extraneous columns: "drop.cols":

        ## Perform the qr-decomposition of X using LINPACK method,
        ## as we need the "good" pivots (and the same as lm()):
        ## FIXME: strongly prefer rankMatrix(X, method= "qr.R")
        qr.X <- qr(X, tol = tol, LAPACK = FALSE)
        rnkX <- qr.X$rank
        if (rnkX == p)
            return(X) ## return X if X has full column rank
        if (kind != "silent.drop.cols") { ## message about no. dropped columns:
	    msg <- sprintf(ngettext(p - rnkX,
			"fixed-effect model matrix is rank deficient so dropping %d column / coefficient",
			"fixed-effect model matrix is rank deficient so dropping %d columns / coefficients"),
			   p - rnkX)
            (if(kind == "warn+drop.cols") warning else message)(msg, domain = NA)
        }
        ## Save properties of X
        contr <- attr(X, "contrasts")
        asgn <- attr(X, "assign")

        ## Return the columns correponding to the first qr.x$rank pivot
        ## elements of X:
        keep <- qr.X$pivot[seq_len(rnkX)]
        X <- X[, keep, drop = FALSE]
	if (rankMatrix(X, tol=tol, method=method) < ncol(X))
            stop(gettextf("Dropping columns failed to produce full column rank design matrix"),
                 call. = FALSE)

        ## Re-assign relevant attributes:
        if(!is.null(contr)) attr(X, "contrasts") <- contr
        if(!is.null(asgn))  attr(X, "assign")    <- asgn[keep]
    }
    X
}

## NA predict and restore rownames of original data if necessary
napredictx <- function(x,...) {
    res <- napredict(x)
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

    cstr <- "check.formula.LHS"
    checkCtrlLevels(cstr,control[[cstr]])
    denv <- checkFormulaData(formula,data,checkLHS=(control$check.formula.LHS=="stop"))
    #mc$formula <- formula <- as.formula(formula,env=denv) ## substitute evaluated call
    formula <- as.formula(formula,env=denv)
    ## as.formula ONLY sets environment if not already explicitly set ...
    ## ?? environment(formula) <- denv
    # get rid of || terms so update() works as expected
    RHSForm(formula) <- expandDoubleVerts(RHSForm(formula))
    mc$formula <- formula

    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substitute "|" by "+"
    environment(fr.form) <- environment(formula)
    ## (DRY! copied from glFormula)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in ...),
    ## so they have to be put there ...
    for (i in c("weights", "offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i,get(i,parent.frame()),environment(fr.form))
    }
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    ## store full, original formula & offset
    attr(fr,"formula") <- formula
    attr(fr,"offset") <- mf$offset
    n <- nrow(fr)
    ## random effects and terms modules
    reTrms <- mkReTrms(findbars(RHSForm(formula)), fr)
    checkNlevels(reTrms$flist, n=n, control)
    checkZdims(reTrms$Ztlist, n=n, control, allow.n=FALSE)
    if (any(is.na(reTrms$Zt))) {
        stop("NA in Z (random-effects model matrix): ",
             "please use ",
             shQuote("na.action='na.omit'"),
             " or ",
             shQuote("na.action='na.exclude'"))
    }
    checkZrank(reTrms$Zt, n=n, control, nonSmall = 1e6)

    ## fixed-effects model matrix X - remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- if(is.null(nb <- nobars(RHSForm(fixedform)))) 1 else nb
    mf$formula <- fixedform
    ## re-evaluate model frame to extract predvars component
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"),"predvars.fixed") <-
        attr(attr(fixedfr,"terms"),"predvars")
    X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    ## backward compatibility (keep no longer than ~2015):
    if(is.null(rankX.chk <- control[["check.rankX"]]))
        rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    if(is.null(scaleX.chk <- control[["check.scaleX"]]))
        scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    X <- checkScaleX(X, kind=scaleX.chk)

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
            if (length((badComp <- setdiff(names(start),c("theta","fixef")))) > 0) {
                stop("incorrect components in start list: ",badComp)
            }
            if (!is.null(start$theta)) theta <- start$theta
            noFixef <- is.null(start$fixef)
            noBeta <- is.null(start$beta)
            if (!noFixef) {
                fixef <- start$fixef
                if (!noBeta) {
                    message("Starting values for fixed effects coefficients",
                            "specified through both 'fixef' and 'beta',",
                            "only 'fixef' used")
                }
            } else if(!noBeta) {
                fixef <- start$beta
            }
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
    ## note: REML does double duty as rank of X and a flag for using
    ## REML maybe this should be mentioned in the help file for
    ## mkRespMod??  currently that help file says REML is logical.  a
    ## consequence of this double duty is that it is impossible to fit
    ## a model with no fixed effects using REML.
    devfun <- mkdevfun(rho, 0L, verbose, control)
    theta <- getStart(start,reTrms$lower,rho$pp)
    if (length(rho$resp$y) > 0)  ## only if non-trivial y
        devfun(rho$pp$theta) # one evaluation to ensure all values are set
    rho$lower <- reTrms$lower # SCW:  in order to be more consistent with mkLmerDevfun
    return(devfun) # this should pass the rho environment implicitly
}


##' @param devfun a deviance function, as generated by \code{\link{mkLmerDevfun}}
##' @return \bold{optimizeLmer}: Results of an optimization.
optimizeLmer <- function(devfun,
                         optimizer=    formals(lmerControl)$optimizer,
                         restart_edge= formals(lmerControl)$restart_edge,
                         boundary.tol = formals(lmerControl)$boundary.tol,
                         start = NULL,
                         verbose = 0L,
                         control = list(),
                         ...) {
    verbose <- as.integer(verbose)
    rho <- environment(devfun)
    opt <- optwrap(optimizer,
                   devfun,
                   getStart(start,rho$lower,rho$pp),
                   lower=rho$lower,
                   control=control,
                   adj=FALSE, verbose=verbose,
                   ...)


    if (restart_edge) {
        ## FIXME: should we be looking at rho$pp$theta or opt$par
        ##  at this point???  in koller example (for getData(13)) we have
        ##   rho$pp$theta=0, opt$par=0.08
        if (length(bvals <- which(rho$pp$theta==rho$lower)) > 0) {
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
            if (any(bgrad < 0)) {
                if (verbose) message("some theta parameters on the boundary, restarting")
                opt <- optwrap(optimizer,
                               devfun,
                               opt$par,
                               lower=rho$lower, control=control,
                               adj=FALSE, verbose=verbose,
                               ...)
            }
        }
    }
    if (boundary.tol > 0) {
        opt <- check.boundary(rho,opt,devfun,boundary.tol)
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

    cstr <- "check.formula.LHS"
    checkCtrlLevels(cstr,control[[cstr]])

    denv <- checkFormulaData(formula, data,
			     checkLHS = (control$check.formula.LHS=="stop"))
    mc$formula <- formula <- as.formula(formula,env=denv)    ## substitute evaluated version

    m <- match(c("data", "subset", "weights", "na.action", "offset",
                 "mustart", "etastart"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fr.form <- subbars(formula) # substitute "|" by "+"
    environment(fr.form) <- environment(formula)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in ...),
    ## so they have to be put there ...
    for (i in c("weights", "offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i,get(i,parent.frame()),environment(fr.form))
    }
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    ## store full, original formula & offset
    attr(fr,"formula") <- formula
    attr(fr,"offset") <- mf$offset
    n <- nrow(fr)
    ## random effects and terms modules
    reTrms <- mkReTrms(findbars(RHSForm(formula)), fr)
    ## TODO: allow.n = !useSc {see FIXME below}
    checkNlevels(reTrms$ flist, n=n, control, allow.n=TRUE)
    checkZdims(reTrms$Ztlist, n=n, control, allow.n=TRUE)
    checkZrank(reTrms$ Zt, n=n, control, nonSmall = 1e6, allow.n=TRUE)

    ## FIXME: adjust test for families with estimated scale parameter:
    ##   useSc is not defined yet/not defined properly?
    ##  if (useSc && maxlevels == n)
    ##          stop("number of levels of each grouping factor must be",
    ##                "greater than number of obs")

    ## fixed-effects model matrix X - remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- if(is.null(nb <- nobars(RHSForm(fixedform))))
        1 else nb
    mf$formula <- fixedform
    ## re-evaluate model frame to extract predvars component
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"),"predvars.fixed") <-
        attr(attr(fixedfr,"terms"),"predvars")
    X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    ## backward compatibility (keep no longer than ~2015):
    if(is.null(rankX.chk <- control[["check.rankX"]]))
        rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    if(is.null(scaleX.chk <- control[["check.scaleX"]]))
        scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    X <- checkScaleX(X, kind=scaleX.chk)

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
    ## allow trivial y
    if (length(y <- rho$resp$y) > 0) {
        if (length(unique(y)) < 2L)
            stop("Response is constant - cannot fit the model")
        rho$verbose     <- as.integer(verbose)

        ## initialize (from mustart)
        .Call(glmerLaplace, rho$pp$ptr(), rho$resp$ptr(), 0L, control$tolPwrss, verbose)
        rho$lp0         <- rho$pp$linPred(1) # each pwrss opt begins at this eta
        rho$pwrssUpdate <- glmerPwrssUpdate
    }
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
                          boundary.tol = formals(glmerControl)$boundary.tol,
                          verbose = 0L,
                          control = list(),
                          nAGQ = 1L,
                          stage = 1,
                          start = NULL,
                          ...) {
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
                   control=control, adj=adj, verbose=verbose,
                   ...)
    if (stage==1) {
        rho$control <- attr(opt,"control")
        rho$nAGQ <- nAGQ
    } else {  ## stage == 2
        rho$resp$setOffset(rho$baseOffset)
    }
    ## FIXME: implement this ...
    if (restart_edge) stop("restart_edge not implemented for optimizeGlmer yet")
    if (boundary.tol > 0) {
        opt <- check.boundary(rho,opt,devfun,boundary.tol)
    }
    return(opt)
}

check.boundary <- function(rho,opt,devfun,boundary.tol) {
    bdiff <- rho$pp$theta-rho$lower[seq_along(rho$pp$theta)]
    if (any(edgevals <- 0 < bdiff & bdiff < boundary.tol)) {
        ## try sucessive "close-to-edge parameters" to see
        ## if we can improve by setting them equal to the boundary
        pp <- opt$par
        for (i in which(edgevals)) {
            tmppar <- pp
            tmppar[i] <- rho$lower[i]
            if (devfun(tmppar) < opt$fval) pp[i] <- tmppar[i]
        }
        opt$par <- pp
        opt$fval <- devfun(opt$par)
        ## re-run to reset (whether successful or not)
        ## FIXME: derivatives, Hessian etc. (and any other
        ## opt messages) *not* recomputed
    }
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

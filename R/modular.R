#### --> ../man/modular.Rd
####     ==================

### Small utilities to be used in lFormula() and glFormula()

doCheck <- function(x) {
    is.character(x) && !any(x == "ignore")
}

RHSForm <- function(form,as.form=FALSE) {
    rhsf <- form[[length(form)]]
    if (as.form) reformulate(deparse(rhsf)) else rhsf
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
    bvals <- c("message","warning","stop","ignore")
    if (smallOK) bvals <- outer(bvals, c("","Small"), paste0)
    if (!is.null(val) && !val %in% bvals)
        stop("invalid control level ",sQuote(val)," in ",cstr,": valid options are {",
             paste(sapply(bvals,sQuote),collapse=","),"}")
    invisible(NULL)
}

## general identifiability checker, used both in checkZdim and checkZrank
wmsg <- function(n, cmp.val, allow.n, msg1="", msg2="", msg3="") {
    if (allow.n) { ## allow  n == cmp.val
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
    wstr <- sprintf("%s (=%d) %s %s (=%d)%s; the random-effects parameters%s are probably unidentifiable",
                    msg1, n, cmp,msg2, cmp.val,msg3,                    rstr)
    list(unident=unident, wstr=wstr)
}

##' For each r.e. term, test if Z has more columns than rows to detect
##' unidentifiability:
##' @title
##' @param Ztlist list of Zt matrices - one for each r.e. term
##' @param n no. observations
##' @param ctrl
##' @param allow.n allow as many random-effects as there are observations
##' for each term?
##' @return possibly empty character string with warning messages
checkZdims <- function(Ztlist, n, ctrl, allow.n=FALSE) {
    stopifnot(is.list(Ztlist), is.numeric(n))
    cstr <- "check.nobs.vs.nRE"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    term.names <- names(Ztlist)
    rows <- vapply(Ztlist, nrow, 1L)
    cols <- vapply(Ztlist, ncol, 1L)
    stopifnot(all(cols == n))
    if (doCheck(cc)) {
        unique(unlist(lapply(seq_along(Ztlist), function(i) {
            ww <- wmsg(cols[i], rows[i], allow.n, "number of observations",
                       "number of random effects",
                       sprintf(" for term (%s)", term.names[i]))
            if(ww$unident) {
                switch(cc,
                       "warning" = warning(ww$wstr, call.=FALSE),
                       "stop"    = stop   (ww$wstr, call.=FALSE),
                       stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
                ww$wstr
            } else character()
        }))) ## -> possibly empty vector of error messages
    } else character()
}

##' @importFrom Matrix rankMatrix
checkZrank <- function(Zt, n, ctrl, nonSmall = 1e6, allow.n=FALSE)
{
    stopifnot(is.list(ctrl), is.numeric(n), is.numeric(nonSmall))
    cstr <- "check.nobs.vs.rankZ"
    if (doCheck(cc <- ctrl[[cstr]])) { ## not NULL or "ignore"
        checkCtrlLevels(cstr, cc, smallOK=TRUE)
        d <- dim(Zt)
        doTr <- d[1L] < d[2L] # Zt is "wide" => qr needs transpose(Zt)
        if(!(grepl("Small",cc) && prod(d) > nonSmall)) {
            rankZ <- rankMatrix(if(doTr) t(Zt) else Zt, method="qr")
            ww <- wmsg(n, rankZ, allow.n, "number of observations", "rank(Z)")
            if(is.na(rankZ)) {
                cc <- "stop"
                ww <-
                    list(unident = TRUE,
                         wstr = sub("^.*;",
                                    "rank(Z) is NA: invalid random effect factors?",
                                    ww$wstr))
            }
            if (ww$unident) {
                switch(cc,
                       "warningSmall" =, "warning" = warning(ww$wstr,call.=FALSE),
                       "stopSmall" =, "stop" = stop(ww$wstr,call.=FALSE),
                       stop(gettextf("unknown check level for '%s'", cstr),
                            domain=NA))
                ww$wstr
            } else character()
        } else character()
    } else character()
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
    ## cstr <- "check.scaleX"
    kinds <- eval(formals(lmerControl)[["check.scaleX"]])
    if (!kind %in% kinds) stop(gettextf("unknown check-scale option: %s",kind))
    if (is.null(kind) || kind == "ignore") return(X)
    ## else :
    cont.cols <- apply(X,2,function(z) !all(z %in% c(0,1)))
    col.sd <- apply(X[,cont.cols, drop=FALSE], 2L, sd)
    sdcomp <- outer(col.sd,col.sd,"/")
    logcomp <- abs(log(sdcomp[lower.tri(sdcomp)]))
    logsd <- abs(log(col.sd))
    if (any(c(logcomp,logsd) > log(tol))) {
        wmsg <- "Some predictor variables are on very different scales:"
        if (kind %in% c("warning","stop")) {
            wmsg <- paste(wmsg, "consider rescaling")
            switch(kind,
                   "warning" = warning(wmsg, call.=FALSE),
                   "stop" = stop(wmsg, call.=FALSE))
        } else {
            wmsg <- paste(wmsg, "auto-rescaled (results NOT adjusted)")
            ## mimic scale() because we don't want to make a copy in
            ##  order to retrieve the center/scale
            X[,cont.cols] <- sweep(X[,cont.cols,drop=FALSE],2,col.sd,"/")
            attr(X,"scaled:scale") <- setNames(col.sd,colnames(X)[cont.cols])
            if (kind == "warn+rescale") warning(wmsg, call.=FALSE)
        }
    } else
        wmsg <- character()
    structure(X, msgScaleX = wmsg)
}


##' @title Check that grouping factors have at least 2 and </<= nobs(.) levels
##' @param flist mkReTrms(.)$flist
##' @param n
##' @param ctrl
##' @param allow.n
##' @return
checkNlevels <- function(flist, n, ctrl, allow.n=FALSE)
{
    stopifnot(is.list(ctrl), is.numeric(n))
    nlevelVec <- vapply(flist, function(x) nlevels(factor(x, exclude=NA)), 1)
    ## Part 1 ----------------
    cstr <- "check.nlev.gtr.1"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    if (doCheck(cc) && any(nlevelVec < 2)) {
        wstr <- "grouping factors must have > 1 sampled level"
        switch(cc,
               "warning" = warning(wstr,call.=FALSE),
               "stop"    =    stop(wstr,call.=FALSE),
               ## FIXME: should never get here since we have checkCtrLevels test above?
               stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
    } else wstr <- character()
    ## Part 2 ----------------
    cstr <- "check.nobs.vs.nlev"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    if (doCheck(cc) && any(if(allow.n) nlevelVec > n else nlevelVec >= n)) {
        ## figure out which factors are the problem?
        w <- if (allow.n) which(nlevelVec>n) else which(nlevelVec>=n)
        bad_facs <- names(nlevelVec)[w]
        wst2 <- gettextf(
            "number of levels of each grouping factor must be %s number of observations",
            if(allow.n) "<=" else "<")
        wst2 <- paste0(wst2," (problems: ",paste(bad_facs,collapse=", "),")")
        switch(cc,
               "warning" = warning(wst2, call.=FALSE),
               "stop"    =    stop(wst2, call.=FALSE)
               ## shouldn't reach here
               )
    } else wst2 <- character()

    ## Part 3 ----------------
    cstr <- "check.nlev.gtreq.5"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    if (doCheck(cc) && any(nlevelVec < 5)) {
        wst3 <- "grouping factors with < 5 sampled levels may give unreliable estimates"
        switch(cc,
               "warning" = warning(wst3, call.=FALSE),
               "stop"    = stop   (wst3, call.=FALSE),
               stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
    } else wst3 <- character()
    ## return:
    c(wstr, wst2, wst3) ## possibly == character(0)
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
##' @importFrom Matrix rankMatrix
##' @author Rune Haubo Bojesen Christensen (drop.coef()); Martin Maechler
chkRank.drop.cols <- function(X, kind, tol = 1e-7, method = "qr") {
    ## Test and match arguments:
    stopifnot(is.matrix(X)) # i.e., *not* sparse
    kinds <- eval(formals(lmerControl)[["check.rankX"]])
    if (!kind %in% kinds) stop(gettextf("undefined option for 'kind': %s", kind))
    ## c("message+drop.cols", "ignore",
    ##   "silent.drop.cols", "warn+drop.cols", "stop.deficient"),

    if(kind == "ignore") return(X)
    ## else :
    p <- ncol(X)
    if (kind == "stop.deficient") {
        if ((rX <- rankMatrix(X, tol=tol, method=method)) < p)
            stop(gettextf(sub("\n +", "\n",
                  "the fixed-effects model matrix is column rank deficient (rank(X) = %d < %d = p);
                   the fixed effects will be jointly unidentifiable"),
                          rX, p), call. = FALSE)
    } else {
        ## kind is one of "message+drop.cols", "silent.drop.cols", "warn+drop.cols"
        ## --> consider to drop extraneous columns: "drop.cols":

        ## Perform the qr-decomposition of X using LINPACK method,
        ## as we need the "good" pivots (and the same as lm()):
        ## this rankMatrix(X, method="qrLINPACK"):  FIXME?  rankMatrix(X, method= "qr.R")
        qr.X <- qr(X, tol = tol, LAPACK = FALSE)
        rnkX <- qr.X$rank
        if (rnkX == p)
            return(X) ## return X if X has full column rank
        ## else:

        ## message about no. dropped columns:
        msg <- sprintf(ngettext(p - rnkX,
                "fixed-effect model matrix is rank deficient so dropping %d column / coefficient",
                "fixed-effect model matrix is rank deficient so dropping %d columns / coefficients"),
                       p - rnkX)
        if (kind != "silent.drop.cols")
            (if(kind == "warn+drop.cols") warning else message)(msg, domain = NA)
        ## Save properties of X
        contr <- attr(X, "contrasts")
        asgn <- attr(X, "assign")

        ## Return the columns correponding to the first qr.x$rank pivot
        ## elements of X:
        keep <- qr.X$pivot[seq_len(rnkX)]
        dropped.names <- colnames(X[,-keep,drop=FALSE])
        X <- X[, keep, drop = FALSE]
        if (rankMatrix(X, tol=tol, method=method) < ncol(X))
            stop(gettextf("Dropping columns failed to produce full column rank design matrix"),
                 call. = FALSE)

        ## Re-assign relevant attributes:
        if(!is.null(contr)) attr(X, "contrasts") <- contr
        if(!is.null(asgn))  attr(X, "assign")    <- asgn[keep]
        attr(X, "msgRankdrop") <- msg
        attr(X, "col.dropped") <- setNames(qr.X$pivot[(rnkX+1L):p],
                                           dropped.names)
    }
    X
}

# check that response is not constant
checkResponse <- function(y, ctrl) {
  stopifnot(is.list(ctrl))
  cstr <- "check.response.not.const"
  checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
  if (doCheck(cc) && length(unique(y)) < 2L) {
    wstr <- "Response is constant"
    switch(cc,
           "warning" = warning(wstr, call.=FALSE),
           "stop"   =    stop(wstr, call.=FALSE),
           stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
  } else character()
}

##' Extract all warning msgs from a merMod object
.merMod.msgs <- function(x) {
    ## currently only those found with 'X' :
    aX <- attributes(x@pp$X)
    wmsgs <- grep("^msg", names(aX))
    if(any(has.msg <- nchar(Xwmsgs <- unlist(aX[wmsgs])) > 0))
        Xwmsgs[has.msg]
    else
        character()
}


## NA predict and restore rownames of original data if necessary
## napredictx <- function(x,...) {
##     res <- napredict(x)
## }

##' @rdname modular
##' @param control a list giving (for \code{[g]lFormula}) all options (see \code{\link{lmerControl}} for running the model;
##' (for \code{mkLmerDevfun,mkGlmerDevfun}) options for inner optimization step;
##' (for \code{optimizeLmer} and \code{optimize[Glmer}) control parameters for nonlinear optimizer (typically inherited from the \dots argument to \code{lmerControl})
##' @return \bold{lFormula, glFormula}: A list containing components,
##' \item{fr}{model frame}
##' \item{X}{fixed-effect design matrix}
##' \item{reTrms}{list containing information on random effects structure: result of \code{\link{mkReTrms}}}
##' \item{REML}{(lFormula only): logical flag: use restricted maximum likelihood? (Copy of argument.)}
##' @export
lFormula <- function(formula, data=NULL, REML = TRUE,
                     subset, weights, na.action, offset, contrasts = NULL,
                     control=lmerControl(), ...)
{
    control <- control$checkControl ## this is all we really need
    mf <- mc <- match.call()

    dontChk <- c("start", "verbose", "devFunOnly")
    dots <- list(...)
    do.call(checkArgs, c(list("lmer"), dots[!names(dots) %in% dontChk]))
    if (!is.null(dots[["family"]])) {
        ## lmer(...,family=...); warning issued within checkArgs
        mc[[1]] <- quote(lme4::glFormula)
        if (missing(control)) mc[["control"]] <- glmerControl()
        return(eval(mc, parent.frame()))
    }

    cstr <- "check.formula.LHS"
    checkCtrlLevels(cstr,control[[cstr]])
    denv <- checkFormulaData(formula, data,
                             checkLHS = control$check.formula.LHS == "stop")
    #mc$formula <- formula <- as.formula(formula,env=denv) ## substitute evaluated call
    formula <- as.formula(formula, env=denv)
    ## as.formula ONLY sets environment if not already explicitly set.
    ## ?? environment(formula) <- denv
    # get rid of || terms so update() works as expected
    RHSForm(formula) <- reformulas::expandDoubleVerts(RHSForm(formula))
    mc$formula <- formula

    ## (DRY! copied from glFormula)
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    fr.form <- reformulas::subbars(formula) # substitute "|" by "+"
    environment(fr.form) <- environment(formula)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in '...'),
    ## so they have to be put there:
    for (i in c("weights", "offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i,get(i,parent.frame()),environment(fr.form))
    }
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    if (nrow(fr) == 0L) stop("0 (non-NA) cases")
    ## convert character vectors to factor (defensive)
    fr <- factorize(fr.form, fr, char.only=TRUE)
    ## store full, original formula & offset
    attr(fr,"formula") <- formula
    attr(fr,"offset") <- mf$offset
    n <- nrow(fr)
    ## random effects and terms modules
    reTrms <- reformulas::mkReTrms(reformulas::findbars(RHSForm(formula)), fr)
    wmsgNlev <- checkNlevels(reTrms$flist, n=n, control)
    wmsgZdims <- checkZdims(reTrms$Ztlist, n=n, control, allow.n=FALSE)
    if (anyNA(reTrms$Zt)) {
        stop("NA in Z (random-effects model matrix): ",
             "please use ",
             shQuote("na.action='na.omit'"),
             " or ",
             shQuote("na.action='na.exclude'"))
    }
    wmsgZrank <- checkZrank(reTrms$Zt, n=n, control, nonSmall = 1e6)

    ## fixed-effects model matrix X - remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- reformulas::nobars(RHSForm(fixedform))
    mf$formula <- fixedform
    ## re-evaluate model frame to extract predvars component
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"), "predvars.fixed") <-
        attr(attr(fixedfr,"terms"), "predvars")
    ## so we don't have to fart around retrieving which vars we need
    ##  in model.frame(.,fixed.only=TRUE)
    attr(attr(fr,"terms"), "varnames.fixed") <- names(fixedfr)

    ## ran-effects model frame (for predvars)
    ## important to COPY formula (and its environment)?
    ranform <- formula
    RHSForm(ranform) <- reformulas::subbars(RHSForm(reOnly(formula)))
    mf$formula <- ranform
    ranfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"), "predvars.random") <-
        attr(terms(ranfr), "predvars")

    ## FIXME: shouldn't we have this already in the full-frame predvars?
    X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    ## backward compatibility (keep no longer than ~2015):
    if(is.null(rankX.chk <- control[["check.rankX"]]))
        rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    if(is.null(scaleX.chk <- control[["check.scaleX"]]))
        scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    X <- checkScaleX(X, kind=scaleX.chk)

    list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula,
         wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))
}

## utility f'n for checking starting values
getStart <- function(start, pred, returnVal = c("theta","all")) {
    returnVal <- match.arg(returnVal)
    doFixef <- returnVal == "all"
    ## default values
    theta <- pred$theta
    fixef <- pred$delb
    if (is.numeric(start)) {
        theta <- start
    } else if (is.list(start)) {
        if (!all(vapply(start, is.numeric, NA)))
            stop("all elements of start must be numeric")
        if (length((badComp <- setdiff(names(start), c("theta","fixef")))) > 0)
            stop("incorrect components in start list: ", badComp)
        if (!is.null(start$theta)) theta <- start$theta
        if (doFixef) {
            noBeta <- is.null(start$beta)
            if(!is.null(sFE <- start$fixef) || !noBeta) {
                fixef <-
                    if(!is.null(sFE)) {
                        if(!noBeta)
                            message("Starting values for fixed effects coefficients",
                                    "specified through both 'fixef' and 'beta',",
                                    "only 'fixef' used")
                        ## FIXME? accumulating heuristic evidence for drop1()-like case
                        if((p <- length(fixef)) < length(sFE) && p == ncol(pred$X) &&
                           is.character(ns <- names(sFE)) &&
                           all((cnX <- dimnames(pred$X)[[2L]]) %in% ns))
                            ## take "matching" fixef[] only
                            sFE[cnX]
                        else
                            sFE
                    }
                    else if(!noBeta)
                        start$beta
            }
            if (length(fixef)!=length(pred$delb))
                stop("incorrect number of fixef components (!=",length(pred$delb),")")
        }
    }
    else if (!is.null(start))
        stop("'start' must be NULL, a numeric vector or named list of such vectors")
    if (!is.null(start) && length(theta) != length(pred$theta))
        stop("incorrect number of theta components (!=",length(pred$theta),")")

    if(doFixef) c(theta, fixef) else theta
}

## update start
## should refactor this to
##  turn numeric start into start=list(theta=start) immediately ??
updateStart <- function(start, theta) {
    if (is.numeric(start)) {
        theta
    } else {
        if (!is.null(start$theta))
            start$theta <- theta
        start
    }
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
mkLmerDevfun <- function(fr, X, reTrms, REML = TRUE, start = NULL,
                         verbose = 0, control = lmerControl(), ...)
{
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
    rho$resp <-
        if(missing(fr))
             mkRespMod(    REML = REMLpass, ...)
        else mkRespMod(fr, REML = REMLpass)
    ## FIXME / note: REML does double duty as rank of X and a flag for using
    ## REML maybe this should be mentioned in the help file for
    ## mkRespMod??  currently that help file says REML is logical.  a
    ## consequence of this double duty is that it is impossible to fit
    ## a model with no fixed effects using REML (MM: ==> FIXME)
    ## devfun <- mkdevfun(rho, 0L, verbose=verbose, control=control)

    ## prevent R CMD check false pos. warnings (in this function only):
    pp <- resp <- NULL
    rho$lmer_Deviance <- lmer_Deviance
    devfun <- function(theta)
        .Call(lmer_Deviance, pp$ptr(), resp$ptr(), as.double(theta))
    environment(devfun) <- rho

    # if all random effects are of the form 1|f and starting values not
    # otherwise provided (and response variable is present, i.e. not doing
    # a simulation) then compute starting values
    if (is.null(start) &&
        all(reTrms$cnms == "(Intercept)") &&
        length(reTrms$flist) == length(reTrms$lower) &&
        !is.null(y <- model.response(fr))) {
        v <- sapply(reTrms$flist, function(f) var(ave(y, f)))
        v.e <- var(y) - sum(v)
        if (!is.na(v.e) && v.e > 0) {
            v.rel <- v / v.e
            if (all(v.rel >= reTrms$lower^2)) rho$pp$setTheta(sqrt(v.rel))
        }
    }

    ## theta <- getStart(start, rho$pp)
    ## ^^^^^ unused / obfuscation? should the above be rho$pp$setTheta(.) ?
    ## MM: commenting it did not break any of our checks
    if (length(rho$resp$y) > 0)  ## only if non-trivial y
        devfun(rho$pp$theta) # one evaluation to ensure all values are set
    rho$lower <- reTrms$lower # to be more consistent with mkGlmerDevfun
    devfun # this should pass the rho environment implicitly
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
                   getStart(start, rho$pp),
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
            theta0 <- new("numeric",rho$pp$theta) ## 'deep' copy:
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
            if (any(is.na(bgrad))) {
                warning("some gradient components are NA near boundaries, skipping boundary check")
                return(opt)
            } else {
                if (any(bgrad < 0)) {
                    if (verbose) message("some theta parameters on the boundary, restarting")
                    opt <- optwrap(optimizer,
                                   devfun,
                                   opt$par,
                                   lower=rho$lower, control=control,
                                   adj=FALSE, verbose=verbose,
                                   ...)
                }
            } ## bgrad not NA
        }
    } ## if restart.edge
    if (boundary.tol > 0)
        check.boundary(rho, opt, devfun, boundary.tol)
    else
        opt
}

## TODO: remove any arguments that aren't actually used by glFormula (same for lFormula)
## TODO(?): lFormula() and glFormula()  are very similar: merge or use common baseFun()
##' @rdname modular
##' @inheritParams glmer
##' @export
glFormula <- function(formula, data=NULL, family = gaussian,
                      subset, weights, na.action, offset,
                      contrasts = NULL, start, mustart, etastart,
                      control = glmerControl(), ...) {
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

    dontChk <- c("verbose", "devFunOnly", "optimizer", "nAGQ")
    dots <- list(...)
    do.call(checkArgs, c(list("glmer"), dots[!names(dots) %in% dontChk]))

    cstr <- "check.formula.LHS"
    checkCtrlLevels(cstr, control[[cstr]])

    denv <- checkFormulaData(formula, data,
                             checkLHS = control$check.formula.LHS == "stop")
    mc$formula <- formula <- as.formula(formula, env = denv)    ## substitute evaluated version

    ## DRY ...
    m <- match(c("data", "subset", "weights", "na.action", "offset",
                 "mustart", "etastart"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    fr.form <- reformulas::subbars(formula) # substitute "|" by "+"
    environment(fr.form) <- environment(formula)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in '...'),
    ## so they have to be put there:
    for (i in c("weights", "offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i, get(i, parent.frame()), environment(fr.form))
    }
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    ## convert character vectors to factor (defensive)
    fr <- factorize(fr.form, fr, char.only = TRUE)
    ## store full, original formula & offset
    attr(fr,"formula") <- formula
    attr(fr,"offset") <- mf$offset
    ## attach starting coefficients to model frame so we can
    ##  pass them through to mkRespMod -> family()$initialize ...
    if (!missing(start) && is.list(start)) {
        fixef <- start$fixef %||% start$beta
        attr(fr,"start") <- fixef
    }
    n <- nrow(fr)
    ## random effects and terms modules
    reTrms <- reformulas::mkReTrms(reformulas::findbars(RHSForm(formula)), fr)
    ## TODO: allow.n = !useSc {see FIXME below}
    wmsgNlev <- checkNlevels(reTrms$ flist, n = n, control, allow.n = TRUE)
    wmsgZdims <- checkZdims(reTrms$Ztlist, n = n, control, allow.n = TRUE)
    wmsgZrank <- checkZrank(reTrms$ Zt, n = n, control, nonSmall = 1e6, allow.n = TRUE)

    ## FIXME: adjust test for families with estimated scale parameter:
    ##   useSc is not defined yet/not defined properly?
    ##  if (useSc && maxlevels == n)
    ##          stop("number of levels of each grouping factor must be",
    ##                "greater than number of obs")

    ## fixed-effects model matrix X - remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- reformulas::nobars(RHSForm(fixedform))
    mf$formula <- fixedform
    ## re-evaluate model frame to extract predvars component
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"),"predvars.fixed") <-
        attr(attr(fixedfr,"terms"),"predvars")

    ## ran-effects model frame (for predvars)
    ## important to COPY formula (and its environment)?
    ranform <- formula
    RHSForm(ranform) <- reformulas::subbars(RHSForm(reOnly(formula)))
    mf$formula <- ranform
    ranfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"), "predvars.random") <-
        attr(terms(ranfr), "predvars")

    X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet
    ## backward compatibility (keep no longer than ~2015):
    if(is.null(rankX.chk <- control[["check.rankX"]]))
        rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    if(is.null(scaleX.chk <- control[["check.scaleX"]]))
        scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    X <- checkScaleX(X, kind=scaleX.chk)

    list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
         wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))
}

##' @rdname modular
##' @export
mkGlmerDevfun <- function(fr, X, reTrms, family, nAGQ = 1L, verbose = 0L,
                          maxit = 100L, control = glmerControl(), ...) {
    stopifnot(length(nAGQ <- as.integer(nAGQ)) == 1L,
              0L <= nAGQ, nAGQ <= 25L)
    verbose <- as.integer(verbose)
    maxit   <- as.integer(maxit)
    rho <- list2env(list(verbose=verbose, maxit=maxit,
                         tolPwrss= control$tolPwrss,
                         compDev = control$compDev),
                    parent = parent.frame())
    rho$pp <- do.call(merPredD$new,
                      c(reTrms[c("Zt","theta","Lambdat","Lind")],
                        n=nrow(X), list(X=X)))
    rho$resp <- if (missing(fr))
        mkRespMod(family=family, ...)
    else
        mkRespMod(fr, family=family)
    nAGQinit <- if(control$nAGQ0initStep) 0L else 1L
    ## allow trivial y
    if (length(y <- rho$resp$y) > 0) {
        checkResponse(y, control$checkControl)
        rho$verbose <- as.integer(verbose)

        ## initialize (from mustart)
        .Call(glmerLaplace, rho$pp$ptr(), rho$resp$ptr(), nAGQinit,
              control$tolPwrss, maxit, verbose)
        rho$lp0         <- rho$pp$linPred(1) # each pwrss opt begins at this eta
        rho$pwrssUpdate <- glmerPwrssUpdate
    }
    rho$lower <- reTrms$lower     # not needed in rho?
    mkdevfun(rho, nAGQinit, maxit=maxit, verbose=verbose, control=control)
    ## this should pass the rho environment implicitly
}



##' @rdname modular
##' @param nAGQ number of Gauss-Hermite quadrature points
##' @param stage optimization stage (1: nAGQ=0, optimize over theta only; 2: nAGQ possibly >0, optimize over theta and beta)
##' @export
optimizeGlmer <- function(devfun,
                          optimizer = if(stage == 1) "bobyqa" else "Nelder_Mead",
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
    if (stage == 1) {
        start <- getStart(start, rho$pp, "theta")
        adj <- FALSE
    } else { ## stage == 2
        start <- getStart(start, rho$pp, "all")
        adj <- TRUE
    }
    opt <- optwrap(optimizer, devfun, start, rho$lower,
                   control=control, adj=adj, verbose=verbose,
                   ...)
    if (stage == 1) {
        rho$control <- attr(opt,"control")
        rho$nAGQ <- nAGQ
    } else {  ## stage == 2
        rho$resp$setOffset(rho$baseOffset)
    }
    if (restart_edge) ## FIXME: implement this ...
        stop("restart_edge not implemented for optimizeGlmer yet")

    if (boundary.tol > 0) {
        opt <- check.boundary(rho, opt, devfun, boundary.tol)
        if(stage != 1) rho$resp$setOffset(rho$baseOffset)
    }
    opt
}

check.boundary <- function(rho,opt,devfun,boundary.tol) {
    bdiff <- rho$pp$theta - rho$lower[seq_along(rho$pp$theta)]
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
    if (nAGQ > 1L) {
        if (length(reTrms$flist) != 1L || length(reTrms$cnms[[1]]) != 1L)
            stop("nAGQ > 1 is only available for models with a single, scalar random-effects term")
    }
    rho <- environment(devfun)
    rho$nAGQ       <- nAGQ
    rho$lower      <- c(rho$lower, rep.int(-Inf, length(rho$pp$beta0)))
    rho$lp0        <- rho$pp$linPred(1)
    rho$dpars      <- seq_along(rho$pp$theta)
    rho$baseOffset <- forceCopy(rho$resp$offset) # forcing a copy (!)
    rho$GQmat      <- GHrule(nAGQ)
    rho$fac        <- reTrms$flist[[1]]
    mkdevfun(rho, nAGQ)  # does this attach rho to devfun??
}

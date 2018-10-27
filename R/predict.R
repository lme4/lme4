##' test for no-random-effect specification: TRUE if NA or ~0, other
##' possibilities are NULL or a non-trivial formula
isRE <- function(re.form) {
    isForm <- inherits(re.form, "formula")
    (is.null(re.form) || isForm || !is.na(re.form)) &&
    (!isForm || length(re.form) != 2L || !identical(re.form[[2L]], 0))
}

if(FALSE) {
isRE(NA)   ## FALSE
isRE(~0)   ##  "
isRE(~y+x) ## TRUE
isRE(NULL) ##  "
isRE(~0+x) ##  "
}

##' Random Effects formula only
reOnly <- function(f, response=FALSE) {
    reformulate(paste0("(", vapply(findbars(f), safeDeparse, ""), ")"),
                response = if(response && length(f)==3L) f[[2]])
}

reFormHack <- function(re.form,ReForm,REForm,REform) {
    warnDeprec <- function(name)
        warning(gettextf("'%s' is deprecated; use '%s' instead", name, "re.form"),
                call.=FALSE, domain=NA)
    if (!missing(ReForm)) {
        warnDeprec("ReForm")
        return(ReForm)
    }
    if (!missing(REForm)) {
        warnDeprec("REForm")
        return(REForm)
    }
    if (!missing(REform)) {
        warnDeprec("REform")
        return(REform)
    }
    re.form
}

## '...' may contain fixed.only=TRUE, random.only=TRUE, ..
get.orig.levs <- function(object, FUN=levels, newdata=NULL, ...) {
    Terms <- terms(object,...)
    mf <- model.frame(object, ...)
    isFac <- vapply(mf, is.factor, FUN.VALUE=TRUE)
    ## ignore response variable
    isFac[attr(Terms,"response")] <- FALSE
    mf <- mf[isFac]
    orig_levs <- if(any(isFac)) lapply(mf, FUN) # else NULL
    ## if necessary (allow.new.levels ...) add in new levels
    if (!is.null(newdata)) {
        for (n in names(mf)) {
            orig_levs[[n]] <- c(orig_levs[[n]],
                                setdiff(unique(as.character(newdata[[n]])),orig_levs[[n]]))
        }
    }
    ## more clues about factor-ness of terms
    if(!is.null(orig_levs)) attr(orig_levs,"isFac") <- isFac
    orig_levs
}


##' Force new parameters into a merMod object
##' (should this be an S3 method?):   params(object) <- newvalue
##' @param object
##' @param params a list of the form specified by \code{start} in
##' \code{\link{lmer}} or \code{\link{glmer}} (i.e. a list containing
##' theta and/or beta; maybe eventually further parameters ...
##' What updating do we have to do in order to make the resulting object
##' consistent/safe?
##' TODO: What kind of checking of input do we have to do?  (Abstract from
##' lmer/glmer code ...)
##' TODO: make sure this gets updated when the parameter structure changes
##' from (theta, beta) to alpha=(theta, beta, phi)
##' @param inplace logical specifying if object should be modified in place; not yet
##' @param subset logical; needs to be true, if only parts of params are to be reset
##' @examples
##' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##' fm1M <- setParams(fm1,list(theta=rep(1,3)))
##' getME(fm1M,"theta")
##' getME(fm1,"theta") ## check that original didn't get messed up
##' ## check that @resp and @pp are the only reference class slots ...
##' sapply(slotNames(fm1),function(x) class(slot(fm1,x)))
setParams <- function(object, params, inplace=FALSE, subset=FALSE) {
    pNames <- c("beta","theta")
    if (object@devcomp$dims["useSc"]) pNames <- c(pNames, "sigma")
    if (!is.list(params) || length(setdiff(names(params),pNames)) > 0)
        stop("params should be specifed as a list with elements from ",
             "{",paste(shQuote(pNames),collapse=", "),"}")
    if (!subset && length(setdiff(pNames,names(params))) > 0) {
        warning("some parameters not specified in setParams()")
    }
    nbeta <- length(object@pp$beta(1))
    ntheta <- length(object@pp$theta)
    if (!is.null(beta <- params$beta) && length(beta)!=nbeta)
        stop("length mismatch in beta (",length(beta),
             "!=",nbeta,")")
    if (!is.null(theta <- params$theta) && length(theta)!=ntheta)
        stop("length mismatch in theta (",length(theta),
             "!=",ntheta,")")
    matchNames <- function(x,tn,vecname="theta") {
        if (!is.null(pn <- names(x))) {
            if (!setequal(pn,tn)) {
                ## pn.not.tn <- setdiff(pn,tn)
                ## tn.not.pn <- setdiff(tn,pn)
                ## TO DO: more detail?
                stop("mismatch between ",shQuote(vecname)," parameter vector names and internal names (",
                     paste(tn,collapse=","),")")
            }
            x <- x[tn]  ## reorder
        } else {
            message(vecname," parameter vector not named: assuming same order as internal vector")
        }
        x
    }
    theta <- matchNames(theta,tnames(object),"theta")
    beta <- matchNames(beta,colnames(getME(object,"X")),"beta")
    sigma <- params$sigma
    if(inplace) {
        stop("modification in place (copy=FALSE) not yet implemented")
    } else { ## copy :
        newObj <- object
        ## make a copy of the reference class slots to
        ##  decouple them from the original object
        newObj@pp <- newObj@pp$copy()
        newObj@resp <- newObj@resp$copy()
        if (!is.null(beta)) {
            newObj@pp$setBeta0(beta)
            newObj@beta <- beta
        }
        if (!is.null(theta)) {
            ## where does theta live and how do I set it?
            ## (1) in .@theta
            ## (2) in .@pp$theta
            newObj@theta <- theta
            newObj@pp$setTheta(theta)
        }
        if (!is.null(sigma)) {
            snm <- if (object@devcomp$dims[["REML"]]) "sigmaREML" else "sigmaML"
            newObj@devcomp[["cmp"]][snm] <- sigma
        }
        return(newObj)
    }
}


##' Make new random effect terms from specified object and new data,
##' possibly omitting some random effect terms
##' @param object fitted model object
##' @param newdata (optional) data frame containing new data
##' @param re.form formula specifying random effect terms to include (NULL=all, ~0)
##' @param na.action
##'
##' @note Hidden; _only_ used (twice) in this file
mkNewReTrms <- function(object, newdata, re.form=NULL, na.action=na.pass,
                        allow.new.levels=FALSE,
                        sparse = max(lengths(orig.random.levs)) > 100)
{
    ## construct (fixed) model frame in order to find out whether there are
    ## missing data/what to do about them
    ## need rfd to inherit appropriate na.action; need grouping
    ## variables as well as any covariates that are included
    ## in RE terms
    ## FIXME: mfnew is new data frame, rfd is processed new data
    ##        why do we need both/what is each doing/how do they differ?
    ##        rfd is *only* used in mkReTrms
    ##        mfnew is *only* used for its na.action attribute (!) [fixed only]
    ##        using model.frame would mess up matrix-valued predictors (GH #201)
    fixed.na.action <- NULL
    if (is.null(newdata)) {
        rfd <- mfnew <- model.frame(object)
        fixed.na.action <- attr(mfnew,"na.action")
    } else {
        if (!identical(na.action,na.pass)) {
            ## only need to re-evaluate for NAs if na.action != na.pass
            mfnew <- model.frame(delete.response(terms(object,fixed.only=TRUE)),
                                 newdata, na.action=na.action)
            fixed.na.action <- attr(mfnew,"na.action")
        }
        ## make sure we pass na.action with new data
        ## it would be nice to do something more principled like
        ## rfd <- model.frame(~.,newdata,na.action=na.action)
        ## but this adds complexities (stored terms, formula, etc.)
        ## that mess things up later on ...
        ## rfd <- na.action(get_all_vars(delete.response(terms(object,fixed.only=FALSE)), newdata))
        newdata.NA <- newdata
        if (!is.null(fixed.na.action)) {
            newdata.NA <- newdata.NA[-fixed.na.action,]
        }
        tt <- delete.response(terms(object,random.only=TRUE))
        orig.random.levs <- get.orig.levs(object, random.only=TRUE, newdata=newdata.NA)
        orig.random.cntr <- get.orig.levs(object, random.only=TRUE,
                                          FUN=contrasts, sparse=sparse)
        ## need to let NAs in RE components go through -- they're handled downstream
        if (inherits(re.form,"formula")) {
            ## We use the RE terms *from the original model fit* to construct
            ##  the model frame. This is good for preserving predvars information,
            ##  getting interactions constructed correctly, etc etc etc,
            ##  but can fail if a partial RE specification is used and some of the variables
            ##  in the original RE form are missing from 'newdata' ...
            ##  Fill them in as necessary. Filling in NA is OK - these vars won't actually
            ##  be used later ...
            pv <- attr(tt,"predvars")
            for (i in 2:(length(pv))) {
                missvars <- setdiff(all.vars(pv[[i]]),all.vars(re.form))
                for (mv in missvars) {
                    newdata.NA[[mv]] <- NA
                }
            }
        }

        ## see comments about why suppressWarnings() is needed below ...
        rfd <- suppressWarnings(
            model.frame(tt, newdata.NA, na.action=na.pass, xlev=orig.random.levs))
        ## restore contrasts (why???)
        ## find *factor* variables involved in terms (left-hand side of RE formula): reset their contrasts
        termvars <- unique(unlist(lapply(findbars(formula(object,random.only=TRUE)),
                                  function(x) all.vars(x[[2]]))))
        for (fn in Reduce(intersect, list(
                              names(orig.random.cntr), termvars, names(rfd)))) {
            ## a non-factor grouping variable *may* sneak in here via simulate(...)
            if (!is.factor(rfd[[fn]])) rfd[[fn]] <- factor(rfd[[fn]])
            contrasts(rfd[[fn]]) <- orig.random.cntr[[fn]]
        }

        if (!is.null(fixed.na.action))
            attr(rfd,"na.action") <- fixed.na.action
        ##
        ## ## need terms to preserve info about spline/orthog polynomial bases
        ## attr(rfd,"terms") <- terms(object)
        ## ## ... but variables list messes things up; can we fix it?
        ## vlist <- lapply(all.vars(terms(object)), as.name)
        ## attr(attr(rfd,"terms"),"variables") <-  as.call(c(quote(list), vlist))
        ##
        ## take out variables that appear *only* in fixed effects
        ## all.v <- all.vars(delete.response(terms(object,fixed.only=FALSE)))
        ## ran.v <- vapply(findbars(formula(object)),all.vars,"")
        ## fix.v <- all.vars(delete.response(terms(object,fixed.only=TRUE)))
        ## rfd <- model.frame(delete.response(terms(object,fixed.only=FALSE)),
        ## newdata,na.action=na.action)
    }
    if (inherits(re.form, "formula")) {
        ## DROP values with NAs in fixed effects
        if (length(fixed.na.action) > 0) {
            newdata <- newdata[-fixed.na.action,]
        }
        ## note: mkReTrms automatically *drops* unused levels
        ReTrms <- mkReTrms(findbars(re.form[[2]]), rfd)
        ## update Lambdat (ugh, better way to do this?)
        ReTrms <- within(ReTrms,Lambdat@x <- unname(getME(object,"theta")[Lind]))
        if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA)))
            stop("NAs are not allowed in prediction data",
                 " for grouping variables unless allow.new.levels is TRUE")
        ns.re <- names(re <- ranef(object))
        nRnms <- names(Rcnms <- ReTrms$cnms)
        if (!all(nRnms %in% ns.re))
            stop("grouping factors specified in re.form that were not present in original model")
        new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
        ## fill in/delete levels as appropriate
        re_x <- Map(function(r,n) levelfun(r,n,
                                           allow.new.levels=allow.new.levels),
                    re[names(new_levels)],
                    new_levels)
        ## pick out random effects values that correspond to
        ##  random effects incorporated in re.form ...
        ## NB: Need integer indexing, as nRnms can be duplicated: (age|Subj) + (sex|Subj) :
        re_new <- lapply(seq_along(nRnms), function(i) {
            rname <- nRnms[i]
            if (!all(Rcnms[[i]] %in% names(re[[rname]])))
                stop("random effects specified in re.form that were not present in original model")
            re_x[[rname]][,Rcnms[[i]]]
        })
        re_new <- unlist(lapply(re_new, t))  ## must TRANSPOSE RE matrices before unlisting
        ## FIXME? use vapply(re_new, t, FUN_VALUE=????)
    }
    Zt <- ReTrms$Zt
    attr(Zt, "na.action") <- attr(re_new, "na.action") <- fixed.na.action
    list(Zt=Zt, b=re_new, Lambdat = ReTrms$Lambdat)
}

##' @param x a random effect (i.e., data frame with rows equal to levels, columns equal to terms
##' @param n vector of new levels
levelfun <- function(x,nl.n,allow.new.levels=FALSE) {
    ## 1. find and deal with new levels
    if (!all(nl.n %in% rownames(x))) {
        if (!allow.new.levels) stop("new levels detected in newdata")
        ## create an all-zero data frame corresponding to the new set of levels ...
        nl.n.comb <- union(nl.n, rownames(x))
        newx <- as.data.frame(matrix(0, nrow=length(nl.n.comb), ncol=ncol(x),
                                     dimnames=list(nl.n.comb, names(x))))
        ## then paste in the matching RE values from the original fit/set of levels
        newx[rownames(x),] <- x
        x <- newx
    }
    ## 2. find and deal with missing old levels
    ## ... these should have been dropped when making the Z matrices
    ##     etc. in mkReTrms, so we'd better drop them here to match ...
    if (!all(r.inn <- rownames(x) %in% nl.n)) {
        x <- x[r.inn,,drop=FALSE]
    }
    return(x)
}

##'
##' \code{\link{predict}} method for \code{\linkS4class{merMod}} objects
##'
##' @title Predictions from a model at new data values
##' @param object a fitted model object
##' @param newdata data frame for which to evaluate predictions
##' @param newparams new parameters to use in evaluating predictions
##' @param re.form formula for random effects to condition on.  If \code{NULL},
##' include all random effects; if \code{NA} or \code{~0},
##' include no random effects
##' @param terms a \code{\link{terms}} object - not used at present
##' @param type character string - either \code{"link"}, the default,
##'    or \code{"response"} indicating the type of prediction object returned
##' @param allow.new.levels (logical) if FALSE (default), then any new levels
##'    (or NA values) detected in \code{newdata} will trigger an error; if TRUE, then
##'    the prediction will use the unconditional (population-level)
##'    values for data with previously unobserved levels (or \code{NA}s)
##' @param na.action function determining what should be done with missing values for fixed effects in \code{newdata}. The default is to predict \code{NA}: see \code{\link{na.pass}}.
##' @param ... optional additional parameters.  None are used at present.
##' @return a numeric vector of predicted values
##' @note There is no option for computing standard errors of predictions because it is difficult to define an efficient method that incorporates uncertainty in the variance parameters; we recommend \code{\link{bootMer}} for this task.
##' @examples
##' (gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 |herd), cbpp, binomial))
##' str(p0 <- predict(gm1))            # fitted values
##' str(p1 <- predict(gm1,re.form=NA))  # fitted values, unconditional (level-0)
##' newdata <- with(cbpp, expand.grid(period=unique(period), herd=unique(herd)))
##' str(p2 <- predict(gm1,newdata))    # new data, all RE
##' str(p3 <- predict(gm1,newdata,re.form=NA)) # new data, level-0
##' str(p4 <- predict(gm1,newdata,re.form=~(1|herd))) # explicitly specify RE
##' @method predict merMod
##' @export
predict.merMod <- function(object, newdata=NULL, newparams=NULL,
                           re.form=NULL,
                           ReForm,
                           REForm,
                           REform,
                           random.only=FALSE,
                           terms=NULL, type=c("link","response"),
                           allow.new.levels=FALSE, na.action=na.pass, ...) {
    ## FIXME: appropriate names for result vector?
    ## FIXME: make sure behaviour is entirely well-defined for NA in grouping factors

    ## Dealing with NAs:
    ## we might need to distinguish among
    ##  (i) NAs in original data and in new data
    ##  (ii) na.action possibilities (exclude, fail, omit, pass)
    ##  (iii) na.action setting in original fit and in predict()
    ##  (iii) NAs in (fixed effect) predictors vs RE grouping variables
    ##  (iv) setting of allow.new.level

    ## NAs in original data (in the fixed or random effects)
    ## may lead to a model frame within the
    ## fitted object that is missing rows; if na.exclude was used,
    ## these will need to be reconstituted in the prediction.
    ##
    ## For the most part, 'na.action's used at the predict stage
    ## (i.e. for newdata) will work on NAs *in the fixed effects*
    ## without further intervention; 'na.pass' will automatically
    ## produce NA values in the output, so 'na.exclude' is not really
    ## necessary (but might get specified anyway)
    ##
    ## In the random effects, NAs in newdata will give a population-level
    ## prediction if allow.new.levels is TRUE; if it's FALSE they give
    ## an error (although it could be argued that in that case they
    ## should follow 'na.action' instead ...)

    re.form <- reFormHack(re.form,ReForm,REForm,REform)

    if (length(list(...)) > 0) warning("unused arguments ignored")

    type <- match.arg(type)
    if (!is.null(terms))
        stop("terms functionality for predict not yet implemented")
    if (!is.null(newparams))
        object <- setParams(object,newparams)

    if (is.null(newdata) && is.null(re.form) &&
        is.null(newparams) && !random.only) {
        ## raw predict() call, just return fitted values
        ##   (inverse-link if appropriate)
        if (isLMM(object) || isNLMM(object)) {
            ## make sure we do *NOT* have NAs in fitted object
            pred <- na.omit(fitted(object))
        } else { ## inverse-link
            pred <-  switch(type,response=object@resp$mu, ## == fitted(object),
                            link=object@resp$eta)
            if (is.null(nm <- rownames(model.frame(object))))
                nm <- seq_along(pred)
            names(pred) <- nm
        }
        fit.na.action <- NULL
        ## flow jumps to end for na.predict
    } else { ## newdata and/or re.form and/or newparams and/or random.only specified
        fit.na.action <- attr(object@frame,"na.action")  ## original NA action

        nobs <- if (is.null(newdata)) nrow(object@frame) else nrow(newdata)
        pred <- rep(0,nobs)

        if (!random.only) {
            X <- getME(object, "X")
            X.col.dropped <- attr(X, "col.dropped")
            ## modified from predict.glm ...
            if (is.null(newdata)) {
                ## Use original model 'X' matrix and offset
                ## orig. offset: will be zero if there are no matches ...
                offset <- model.offset(model.frame(object))
                if (is.null(offset)) offset <- 0
            } else {  ## new data specified
                ## evaluate new fixed effect
                RHS <- formula(substitute(~R,
                         list(R=RHSForm(formula(object,fixed.only=TRUE)))))
                ## https://github.com/lme4/lme4/issues/414
                ## contrasts are not relevant in random effects;
                ##  model.frame.default warns about dropping contrasts
                ##  if (1) xlev is specified and (2) any factors in
                ##  original data frame had contrasts set

                ## alternative solution: drop contrasts manually
                ## (could assign to a new variable newdata2 for safety,
                ## but I don't think newdata
                ##  is used downstream in this function?)

                ## isFacND <- which(vapply(newdata, is.factor, FUN.VALUE = TRUE))
                ## for (j in isFacND) {
                ##    attr(newdata[[j]], "contrasts") <- NULL
                ## }

                orig.fixed.levs <- get.orig.levs(object, fixed.only=TRUE)
                mfnew <- suppressWarnings(
                    model.frame(delete.response(terms(object,fixed.only=TRUE)),
                                newdata,
                                na.action = na.action, xlev = orig.fixed.levs))

                X <- model.matrix(RHS, data=mfnew,
                                  contrasts.arg=attr(X,"contrasts"))
                ## hack to remove unused interaction levels?
                ## X <- X[,colnames(X0)]

                offset <- 0 # rep(0, nrow(X))
                tt <- terms(object)
                if (!is.null(off.num <- attr(tt, "offset"))) {
                    for (i in off.num)
                        offset <- offset + eval(attr(tt,"variables")[[i + 1]], newdata)
                }
                ## FIXME?: simplify(no need for 'mfnew'): can this be different from 'na.action'?
                fit.na.action <- attr(mfnew,"na.action")
                ## only need to drop if new data specified ...
                if(is.numeric(X.col.dropped) && length(X.col.dropped) > 0)
                    X <- X[, -X.col.dropped, drop=FALSE]
            }
            pred <- drop(X %*% fixef(object))

            ## FIXME:: need to unname()  ?
            ## FIXME: is this redundant??
            ## if (!is.null(frOffset <- attr(object@frame,"offset")))
            ##     offset <- offset + eval(frOffset, newdata)
            pred <- pred+offset

        } ## end !(random.only)

        if (isRE(re.form)) {
            if (is.null(re.form))
                re.form <- reOnly(formula(object)) # RE formula only
            rfd <- if (is.null(newdata)) object@frame else newdata
            newRE <- mkNewReTrms(object, rfd, re.form, na.action=na.action,
                                 allow.new.levels=allow.new.levels)
            REvals <- base::drop(as(newRE$b %*% newRE$Zt, "matrix"))
            pred <- pred + REvals
            if (random.only) {
                fit.na.action <- attr(newRE$Zt,"na.action")
            }
        }
        if (isGLMM(object) && type=="response") {
            pred <- object@resp$family$linkinv(pred)
        }
    }  ## newdata/newparams/re.form

    ## fill in NAs as appropriate:
    ##   if NAs were detected in original model fit, OR in updated model frame construction
    ## but DON'T double-NA if raw prediction in the first place
    if (is.null(newdata)) {
        fit.na.action <- attr(model.frame(object),"na.action")
        if (!missing(na.action)) {
            ## hack to override action where explicitly specified
            if (!is.null(fit.na.action))
                class(fit.na.action) <- class(attr(na.action(NA),"na.action"))
        }
    }
    napredict(fit.na.action, pred)
} # end {predict.merMod}


##' Simulate responses from the model represented by a fitted model object
##'
simulate.formula <- function(object, nsim = 1, seed = NULL, family,
                             weights=NULL, offset=NULL, ...) {
    ## N.B. *must* name all arguments so that 'object' is missing in .simulateFun()
    .simulateFun(formula=object, nsim=nsim, seed=seed,
                 family=family, weights=weights, offset=offset, ...)
}

simulate.merMod <- function(object, nsim = 1, seed = NULL, use.u = FALSE,
                            re.form=NA, ReForm, REForm, REform,
                            newdata=NULL, newparams=NULL,
                            family=NULL,
                            allow.new.levels=FALSE, na.action=na.pass, ...) {
    mc <- match.call()
    mc[[1]] <- quote(lme4::.simulateFun)
    eval(mc, parent.frame(1L))
}

.simulateFun <- function(object, nsim = 1, seed = NULL, use.u = FALSE,
                         re.form=NA, ReForm, REForm, REform,
                         newdata=NULL, newparams=NULL,
                         formula=NULL,family=NULL,
                         weights=NULL,
                         offset=NULL,
                         allow.new.levels=FALSE,
                         na.action=na.pass,
                         cond.sim=TRUE,
                         ...) {

    nullWts <- FALSE

    if (is.null(weights)) {
        if (is.null(newdata))
            weights <- weights(object)
        else {

            nullWts <- TRUE # this flags that 'weights' wasn't supplied by the user
            weights <- rep(1,nrow(newdata))
        }
    }
    if (missing(object)) {
        if (is.null(formula) || is.null(newdata) || is.null(newparams)) {
            stop("if ",sQuote("object")," is missing, must specify all of ",
                 sQuote("formula"),", ",sQuote("newdata"),", and ",
                 sQuote("newparams"))
        }

        ## construct fake-fitted object from data, params
        ## copied from glm(): DRY; this all stems from the
        ## original sin of handling family=gaussian as a special
        ## case
        if (is.character(family))
            family <- get(family, mode = "function", envir = parent.frame())
        if (is.function(family))
            family <- family()
        if (is.null(family) ||
            (family$family=="gaussian" && family$link=="identity")) {
            lmod <- lFormula(formula,newdata,
                             weights=weights,
                             offset=offset,
                             control=lmerControl(check.formula.LHS="ignore"))
            devfun <- do.call(mkLmerDevfun, lmod)
            object <- mkMerMod(environment(devfun),
                               ## (real parameters will be filled in later)
                               opt = list(par=NA,fval=NA,conv=NA),
                               lmod$reTrms, fr = lmod$fr)
        } else {
            glmod <- glFormula(formula,newdata,family=family,
                               weights=weights,
                               offset=offset,
                               control=glmerControl(check.formula.LHS="ignore"))
            devfun <- do.call(mkGlmerDevfun, glmod)
            object <- mkMerMod(environment(devfun),
                               ## (real parameters will be filled in later)
                               opt = list(par=NA,fval=NA,conv=NA),
                               glmod$reTrms, fr = glmod$fr)
        }
        ## would like to do this:
        ## so predict() -> fitted() -> set default names will work
        ## instead we have a special case in fitted()
        ## object@resp$mu <- rep(NA_real_,nrow(model.frame(object)))
    }
    stopifnot((nsim <- as.integer(nsim[1])) > 0,
              is(object, "merMod"))
    if (!is.null(newparams)) {
        object <- setParams(object,newparams)
    }

    ## need to save this before we reset re.form
    re.form.miss <- missing(re.form)
    re.form <- reFormHack(re.form,ReForm,REForm,REform)

    if (!missing(use.u)) {
        if (!re.form.miss) {
            stop("should specify only one of ",sQuote("use.u"),
                 " and ",sQuote("re.form"))
        }
        re.form <- if (use.u) NULL else ~0
    }
    if (is.null(re.form)) { # formula w/o response
        re.form <- noLHSform(formula(object))
    }
    if(!is.null(seed)) set.seed(seed)
    if(!exists(".Random.seed", envir = .GlobalEnv))
        runif(1) # initialize the RNG if necessary
    RNGstate <- .Random.seed

    sigma <- sigma(object)
    ## OBSOLETE: no longer use X?
    ## n <- nrow(X <- getME(object, "X"))
    ## link <- if (isGLMM(object)) "response"

    ## predictions, conditioned as specified, on link scale
    ## previously: do **NOT** use na.action as specified here (inherit
    ##     from object instead, for consistency)
    ## now: use na.omit, because we have to match up
    ##    with whatever is done in mkNewReTrms
    etapred <- predict(object, newdata=newdata, re.form=re.form,
                       type="link", na.action=na.omit,
                       allow.new.levels=allow.new.levels)
    n <- length(etapred)

    ## now add random components:
    ##  only the ones we did *not* condition on

    ## compre.form <- noLHSform(formula(object))
    ## construct RE formula ONLY: leave out fixed terms,
    ##   which might have loose terms like offsets in them ...

    ##' combine unary or binary operator + arguments (sugar for 'substitute')
    makeOp <- function(x,y,op=NULL) {
        if (is.null(op)) {  ## unary
            substitute(OP(X),list(X=x,OP=y))
        } else substitute(OP(X,Y), list(X=x,OP=op,Y=y))
    }

    compReForm <- reOnly(formula(object))
    if (isRE(re.form)) {
        rr <- reOnly(re.form)[[2]] ## expand RE and strip ~
        ftemplate <- substitute(.~.-XX, list(XX=rr))
        compReForm <- update.formula(compReForm,ftemplate)[-2]
        ## update, then delete LHS
    }

    ## (1) random effect(s)
    sim.reff <- if (!is.null(findbars(compReForm))) {
        newRE <- mkNewReTrms(object, newdata, compReForm,
                             na.action=na.action,
                             allow.new.levels=allow.new.levels)
        ## this *can* justifiably happen, if we are using mkNewReTrms
        ## in the context of predicting/simulating with a non-trivial
        ## re.form ...
        ## <obsolete> paranoia ...
        ## <obsolete> stopifnot(!is.null(newdata) ||
        ##       isTRUE(all.equal(newRE$Lambdat,getME(object,"Lambdat"))))
        U <- t(newRE$Lambdat %*% newRE$Zt) # == Z Lambda
        u <- rnorm(ncol(U)*nsim)
        ## UNSCALED random-effects contribution:
        as(U %*% matrix(u, ncol = nsim), "matrix")
    } else 0

    val <- if (isLMM(object)) {
          ## result will be matrix  n x nsim :
          etapred + sigma * (sim.reff +
                               ## residual contribution:
                               if (cond.sim)
                                   matrix(rnorm(n * nsim), ncol = nsim)
                               else 0)
    } else if (isGLMM(object)) {
        ## GLMM
        ## n.b. DON'T scale random-effects (???)
        etasim <- etapred+sim.reff
        family <- normalizeFamilyName(object@resp$family)
        musim <- family$linkinv(etasim) #-> family$family == "negative.binomial" if(NB)
        ## ntot <- length(musim) ## FIXME: or could be dims["n"]?
        ##

        if (family$family=="binomial" && is.matrix(r <- model.response(object@frame))) {

            # unless the user passed in new weights, take them from the response matrix
            # e.g. cbind(incidence, size-incidence) ~ ...
            if(nullWts) weights <- rowSums(r)
        }

        if (is.null(sfun <- simfunList[[family$family]]) &&
            is.null(family$simulate))
            stop("simulation not implemented for family",
                 family$family)
        ## don't rely on automatic recycling
        if (cond.sim) {
             val <- sfun(object, nsim=1, ftd = rep_len(musim, n*nsim),
                         wts = weights)
        } else {
             val  <- rep_len(musim, n*nsim)
        }
        ## split results into nsims: need special case for binomial matrix/factor responses
        if (family$family=="binomial" && is.matrix(r <- model.response(object@frame))) {
            lapply(split(val[[1]], gl(nsim, n, 2 * nsim * n)), matrix,
                          ncol = 2, dimnames = list(NULL, colnames(r)))
        } else if (family$family=="binomial" && is.factor(val[[1]])) {
            split(val[[1]], gl(nsim,n))
        } else split(val, gl(nsim,n))
    } else
        stop("simulate method for NLMMs not yet implemented")

    ## from src/library/stats/R/lm.R
    if(!is.list(val)) {
        dim(val) <- c(n, nsim)
        val <- as.data.frame(val)
    } else class(val) <- "data.frame"
    names(val) <- paste("sim", seq_len(nsim), sep="_")
    ## have not yet filled in NAs, so need to use names of fitted
    ## object NOT including values with NAs
    f <- fitted(object)
    nm <- names(f)[!is.na(f)]
    ## unnamed input, *or* simulation from new data ...
    if (length(nm) == 0) {
        nm <- as.character(seq(n))
    } else if (!is.null(newdata)) {
        nm <- rownames(newdata)
    }
    row.names(val) <- nm

    fit.na.action <- attr(model.frame(object), "na.action")

    if (!missing(na.action) &&  !is.null(fit.na.action)) {
        ## retrieve name of na.action type ("omit", "exclude", "pass")
        class.na.action <- class(attr(na.action(NA), "na.action"))
        if (class.na.action != class(fit.na.action)) {
            ## hack to override action where explicitly specified
            class(fit.na.action) <- class.na.action
        }
    }

    nafun <- function(x) { x[] <- apply(x,
                                        2L,
                                        napredict,
                                        omit = fit.na.action); x }
    val <- if (is.matrix(val[[1]])) {
        ## have to handle binomial response matrices differently --
        ## fill in NAs as appropriate in *both* columns
        structure(lapply(val, nafun),
                  ## have to put this back into a (weird) data frame again,
                  ## carefully (should do the napredict stuff
                  ## earlier, so we don't have to redo this transformation!)
                  class = "data.frame")
    } else {
        as.data.frame(lapply(val, napredict, omit=fit.na.action))
    }

    ## reconstruct names: first get rid of NAs, then refill them
    ## as appropriate based on fit.na.action (which may be different
    ## from the original model's na.action spec)
    nm2 <-
        if (is.null(newdata))
            names(napredict(na.omit(f), omit=fit.na.action))
        else
            rownames(napredict(newdata, omit=fit.na.action))
    if (length(nm2) > 0)
        row.names(val) <- nm2

    structure(val,
              ## as.data.frame(lapply(...)) blows away na.action attribute,
              ##  so we have to re-assign here
              na.action = fit.na.action,
              seed = RNGstate)
}## .simulateFun()

########################
## modified from stats/family.R
## TODO: the $simulate methods included with R families by default
## are not sufficiently flexible to be re-used by lme4.
## these are modified by:
## (1) adding a 'ftd' argument for the fitted values
##     that defaults to fitted(object), to allow more flexibility
##     e.g. in conditioning on or marginalizing over random effects
##     (fitted(object) can be produced from predict.merMod() with
##     alternative parameters rather than being extracted directly
##     from the fitted objects -- this allows simulation with new
##     parameters or new predictor variables
## (2) modifying wts from object$prior.weights to weights(object)
## (3) adding wts as an argument
##
##
## these can be incorporated by overwriting the simulate()
## components, or calling them
##
gaussian_simfun <- function(object, nsim, ftd=fitted(object),
                            wts=weights(object)) {

    if (any(wts != 1)) warning("ignoring prior weights")
    rnorm(nsim*length(ftd), ftd, sd=sigma(object))
}

binomial_simfun <- function(object, nsim, ftd=fitted(object),
                            wts=weights(object)) {
    n <- length(ftd)
    ntot <- n*nsim
    if (any(wts %% 1 != 0))
        stop("cannot simulate from non-integer prior.weights")
    ## Try to figure out if the original data were
    ## proportions, a factor or a two-column matrix
    if (!is.null(m <- model.frame(object))) {
        y <- model.response(m)
        if(is.factor(y)) {
            ## ignore weights
            yy <- factor(1+rbinom(ntot, size = 1, prob = ftd),
                         labels = levels(y))
            split(yy, rep(seq_len(nsim), each = n))
        } else if(is.matrix(y) && ncol(y) == 2) {
            yy <- vector("list", nsim)
            for (i in seq_len(nsim)) {
                Y <- rbinom(n, size = wts, prob = ftd)
                YY <- cbind(Y, wts - Y)
                colnames(YY) <- colnames(y)
                yy[[i]] <- YY
            }
            yy
        } else
            rbinom(ntot, size = wts, prob = ftd)/wts
    } else rbinom(ntot, size = wts, prob = ftd)/wts
}

poisson_simfun <- function(object, nsim, ftd=fitted(object),
                           wts=weights(object)) {
        ## A Poisson GLM has dispersion fixed at 1, so prior weights
        ## do not have a simple unambiguous interpretation:
        ## they might be frequency weights or indicate averages.
        wts <- weights(object)
        if (any(wts != 1)) warning("ignoring prior weights")
        rpois(nsim*length(ftd), ftd)
    }


##' FIXME: need a gamma.shape.merMod method in order for this to work.
##'        (see initial shot at gamma.shape.merMod below)
Gamma_simfun <- function(object, nsim, ftd=fitted(object),
                         wts=weights(object)) {
    if (any(wts != 1)) message("using weights to scale shape parameter")
    ## used to use gamma.shape(), but sigma() is more general
    ## (wouldn't work *outside* of the merMod context though)
    shape <- sigma(object)*wts
    rgamma(nsim*length(ftd), shape = shape, rate = shape/ftd)
}

gamma.shape.merMod <- function(object, ...) {
    if(family(object)$family != "Gamma")
        stop("Can not fit gamma shape parameter because Gamma family not used")

    y <- getME(object, "y")
    mu <- getME(object, "mu")
    w <- weights(object)
                                        # Sec 8.3.2 (MN)
    L <- w*(log(y/mu)-((y-mu)/mu))
    dev <- -2*sum(L)
                                        # Eqs. between 8.2 & 8.3 (MN)
    Dbar <- dev/length(y)
    structure(list(alpha = (6+2*Dbar)/(Dbar*(6+Dbar)),
                   SE = NA), # FIXME: obtain standard error
              class = "gamma.shape")
}


## FIXME: include without inducing SuppDists dependency?
## inverse.gaussian_simfun <- function(object, nsim, ftd=fitted(object)) {
##     if(is.null(tryCatch(loadNamespace("SuppDists"),
##                         error = function(e) NULL)))
##         stop("need CRAN package 'SuppDists' for the 'inverse.gaussian' family")
##     wts <- weights(object)
##     if (any(wts != 1)) message("using weights as inverse variances")
##     SuppDists::rinvGauss(nsim * length(ftd), nu = ftd,
##                          lambda = wts/summary(object)$dispersion)
## }

## in the original MASS version, .Theta is assigned into the environment
## (triggers a NOTE in R CMD check)
## modified from @aosmith16 GH contribution

negative.binomial_simfun <- function (object, nsim, ftd = fitted(object),
                                          wts=weights(object))
{

    if (any(wts != 1))
        warning("ignoring prior weights")
    theta <- getNBdisp(object)
    rnbinom(nsim * length(ftd), mu = ftd, size = theta)
}


simfunList <- list(gaussian = gaussian_simfun,
                   binomial = binomial_simfun,
                   poisson  = poisson_simfun,
                   Gamma    = Gamma_simfun,
                   negative.binomial = negative.binomial_simfun)

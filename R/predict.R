##' test for no-random-effect specification: TRUE if NA or ~0, other
##' possibilities are NULL or a non-trivial formula
noReForm <- function(re.form) {
    (!is.null(re.form) && !is(re.form,"formula") && is.na(re.form)) ||
        (is(re.form,"formula") && length(re.form)==2 && identical(re.form[[2]],0))
}

reOnly <- function(f) {
    reformulate(paste0("(",sapply(findbars(f),deparse),")"))
}

reFormHack <- function(re.form,ReForm,REForm,REform) {
    if (!missing(ReForm)) {
        message(shQuote("re.form")," is now preferred to ",shQuote("ReForm"))
        return(ReForm)
    }
    if (!missing(REForm)) {
        message(shQuote("re.form")," is now preferred to ",shQuote("REForm"))
        return(REForm)
    }
    if (!missing(REform)) {
        message(shQuote("re.form")," is now preferred to ",shQuote("REform"))
        return(REform)
    }
    re.form
}

## noReForm(NA)   ## TRUE
## noReForm(~0)   ## TRUE
## noReForm(~y+x)
## noReForm(NULL)
## noReForm(~0+x)

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
    if (useSc <- object@devcomp$dims["useSc"]) pNames <- c(pNames,"sigma")
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
##' @param ReForm backward compatibility argument
##' @param na.action
##'
mkNewReTrms <- function(object, newdata, re.form=NULL, ReForm,
                        na.action=na.pass,
                        allow.new.levels=FALSE) {
    ## construct (fixed) model frame in order to find out whether there are
    ## missing data/what to do about them
    re.form <- reFormHack(re.form,ReForm)  ## back-compatibility
    mfnew <- if (is.null(newdata)) {
        model.frame(object)  ## FIXME: check ...
    } else {
        Terms <- terms(object,fixed.only=TRUE)
        model.frame(delete.response(Terms),newdata, na.action=na.action)
    }
    if (is(re.form,"formula")) {
        ## DROP values with NAs in fixed effects
        if (length(fit.na.action <- attr(mfnew,"na.action"))>0) {
            newdata <- newdata[-fit.na.action,]
        }
        re <- ranef(object)
        rfd <- if(is.null(newdata)) object@frame else newdata
        ReTrms <- mkReTrms(findbars(re.form[[2]]), rfd)
        if (!allow.new.levels &&
            any(sapply(ReTrms$flist,function(x) any(is.na(x)))))
            stop("NAs are not allowed in prediction data",
                 " for grouping variables unless allow.new.levels is TRUE")
        unames <- unique(sort(names(ReTrms$cnms)))  ## FIXME: same as names(ReTrms$flist) ?
        ## convert numeric grouping variables to factors as necessary
        ## must use all.vars() for examples
        ## for (i in all.vars(RHSForm(re.form))) {
        ## TO DO: should restrict attention to grouping factors only
        getgrpvars <- function(x) all.vars(x[[3]])
        all.grp.vars <- unique(unlist(lapply(findbars(re.form),getgrpvars)))
        for (i in all.grp.vars) {
            if (!is.matrix(rfd[[i]]))
                rfd[[i]] <- factor(rfd[[i]])
        }
        Rfacs <- lapply(setNames(unames,unames),
                        function(x) factor(eval(parse(text=x),envir=rfd)))
        new_levels <- lapply(Rfacs,function(x) levels(droplevels(factor(x))))
        ## FIXME: should this be unique(as.character(x)) instead?
        ##   (i.e., what is the proper way to protect against non-factors?)
        levelfun <- function(x,n) {
            ## find and deal with new levels
            nl.n <- new_levels[[n]]
            if (!all(nl.n %in% rownames(x))) {
                if (!allow.new.levels) stop("new levels detected in newdata")
                ## create an all-zero data frame corresponding to the new set of levels ...
                newx <- as.data.frame(matrix(0, nrow=length(nl.n), ncol=ncol(x),
                                             dimnames=list(nl.n, names(x))))
                ## then paste in the matching RE values from the original fit/set of levels
                newx[rownames(x),] <- x
                x <- newx
            }
            ## find and deal with missing old levels
	    if (!all(r.inn <- rownames(x) %in% nl.n)) x[r.inn,,drop=FALSE] else x
        }
        ## fill in/delete levels as appropriate
        re_x <- mapply(levelfun, re, names(re), SIMPLIFY=FALSE)
        re_new <- list()
        if (any(!names(ReTrms$cnms) %in% names(re)))
            stop("grouping factors specified in re.form that were not present in original model")
        ## pick out random effects values that correspond to
        ##  random effects incorporated in re.form ...
        for (i in seq_along(ReTrms$cnms)) {
            rname <- names(ReTrms$cnms)[i]
            if (any(!ReTrms$cnms[[rname]] %in% names(re[[rname]])))
                stop("random effects specified in re.form",
                     " that were not present in original model")
            re_new[[i]] <- re_x[[rname]][,ReTrms$cnms[[rname]]]
        }
        re_newvec <- unlist(lapply(re_new,t))  ## must TRANSPOSE RE matrices before unlisting
    }
    Zt <- ReTrms$Zt
    attr(Zt,"na.action") <-
        attr(re_newvec,"na.action") <- attr(mfnew,"na.action")
    list(Zt=Zt,b=re_newvec)
}

##'
##' \code{\link{predict}} method for \code{\linkS4class{merMod}} objects
##'
##' @title Predictions from a model at new data values
##' @param object a fitted model object
##' @param newdata data frame for which to evaluate predictions
##' @param newparams new parameters to use in evaluating predictions
##' @param newX new design matrix to use in evaluating predictions
##' (alternative to \code{newdata})
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
predict.merMod <- function(object, newdata=NULL, newparams=NULL, newX=NULL,
                           re.form=NULL,
                           ReForm,
                           REForm,
                           REform,
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
    
    if (length(list(...)>0)) warning("unused arguments ignored")

    
    fit.na.action <- NULL
    type <- match.arg(type)
    if (!is.null(terms))
        stop("terms functionality for predict not yet implemented")
    if (!is.null(newparams)) {
        object <- setParams(object,newparams)
    }
    if ((rawPred <- is.null(newdata) &&
                    is.null(re.form) && is.null(newparams))) { 
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
        ## flow jumps to end for na.predict
    } else { ## newdata and/or re.form and/or newparams specified
        X_orig <- getME(object, "X")
        ## modified from predict.glm ...
        if (is.null(newdata)) {
            ## get original model matrix and offset
            X <- X_orig
            fit.na.action <- attr(object@frame,"na.action")  ## original NA action
            ## orig. offset: will be zero if there are no matches ...
            offset <- rowSums(object@frame[grepl("offset\\(.*\\)",
                                                 names(object@frame))])

        } else {  ## new data specified
            ## evaluate new fixed effect
            RHS <- formula(substitute(~R,
                              list(R=RHSForm(formula(object,fixed.only=TRUE)))))
            Terms <- terms(object,fixed.only=TRUE)
            isFac <- vapply(mf <- model.frame(object,
                                              fixed.only=TRUE),
                            is,"factor",FUN.VALUE=TRUE)
            ## ignore response variable
            isFac[attr(Terms,"response")] <- FALSE
            orig_levs <- if (length(isFac)==0) NULL else lapply(mf[isFac],levels)
            X <- model.matrix(RHS, mfnew <- model.frame(delete.response(Terms),
                                                        newdata,
                                                        na.action=na.action,
                                                        xlev=orig_levs),
                              contrasts.arg=attr(X_orig,"contrasts"))
            offset <- rep(0, nrow(X))
            tt <- terms(object)
            if (!is.null(off.num <- attr(tt, "offset"))) {
                for (i in off.num)
                    offset <- offset + eval(attr(tt,"variables")[[i + 1]], newdata)
            }
            fit.na.action <- attr(mfnew,"na.action")
        }
        pred <- drop(X %*% fixef(object))
        ## FIXME:: need to unname()  ?
        ## FIXME: is this redundant??
        ## if (!is.null(frOffset <- attr(object@frame,"offset")))
        ##     offset <- offset + eval(frOffset, newdata)
        pred <- pred+offset
        if (!noReForm(re.form)) {
            if (is.null(re.form)) {
		re.form <- reOnly(formula(object)) # RE formula only
                ReTrms <- mkReTrms(findbars(re.form), newdata)
            }
            newRE <- mkNewReTrms(object,newdata,re.form,na.action=na.action,
                                 allow.new.levels=allow.new.levels)
            pred <- pred + drop(as.matrix(newRE$b %*% newRE$Zt))
        }
        if (isGLMM(object) && type=="response") {
            pred <- object@resp$family$linkinv(pred)
        }
    }  ## newdata/newparams/re.form
    ## fill in NAs as appropriate: if NAs were detected in
    ## original model fit, OR in updated model frame construction
    ## but DON'T double-NA if raw prediction in the first place
    old.fit.na.action <- attr(model.frame(object),"na.action")
    if (!is.null(fit.na.action) ||  ## NA used in new construction
        (!is.null(fit.na.action <- old.fit.na.action))) {
        if (!missing(na.action)) {
            ## hack to override action where explicitly specified
            class(fit.na.action) <- class(attr(na.action(NA),"na.action"))
        }
        pred <- napredict(fit.na.action,pred)
    }
    pred
}

##' @importFrom stats simulate
NULL
##' Simulate responses from the model represented by a fitted model object
##'
##' @title Simulate responses from a \code{\linkS4class{merMod}} object
##' @param object a fitted model object
##' @param nsim positive integer scalar - the number of responses to simulate
##' @param seed an optional seed to be used in \code{set.seed} immediately
##'     before the simulation so as to generate a reproducible sample.
##' @param ... optional additional arguments, none are used at present
##' @examples
##' ## test whether fitted models are consistent with the
##' ##  observed number of zeros in CBPP data set:
##' gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'              data = cbpp, family = binomial)
##' gg <- simulate(gm1,1000)
##' zeros <- sapply(gg,function(x) sum(x[,"incidence"]==0))
##' plot(table(zeros))
##' abline(v=sum(cbpp$incidence==0),col=2)
##' @method simulate merMod
##' @export
##'
simulate.formula <- function(object, nsim = 1, seed = NULL, family, weights=NULL, offset=NULL, ...) {
    ## N.B. *must* name all arguments so that 'object' is missing in
    ## .simulateFun
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
                         allow.new.levels=FALSE, na.action=na.pass, ...) {

    if (missing(object)) {
        if (is.null(formula) || is.null(newdata) || is.null(newparams)) {
            stop("if ",sQuote("object")," is missing, must specify all of ",
                 sQuote("formula"),", ",sQuote("newdata"),", and ",
                 sQuote("newparams"))
        }

        if (is.null(weights)) weights <- rep(1,nrow(newdata))
        ## construct fake-fitted object from data, params
        ## copied from glm(): DRY; this all stems from the
        ## original sin of handling family=gaussian as a special
        ## case
        if (is.character(family))
            family <- get(family, mode = "function", envir = parent.frame())
        if (is.function(family))
            family <- family()
        if (is.null(family) || family$family=="gaussian") {
            lmod <- lFormula(formula,newdata,
                             weights=weights,
                             offset=offset,
                             control=lmerControl(check.formula.LHS="ignore"))
            devfun <- do.call(mkLmerDevfun, lmod)
            object <- mkMerMod(environment(devfun),
                               ## (real parameters will be filled in later)
                               opt=list(par=NA,fval=NA,conv=NA),
                               lmod$reTrms,fr=lmod$fr)
        } else {
            glmod <- glFormula(formula,newdata,family=family,
                               weights=weights,
                               offset=offset,
                               control=glmerControl(check.formula.LHS="ignore"))
            devfun <- do.call(mkGlmerDevfun, glmod)
            object <- mkMerMod(environment(devfun),
                               ## (real parameters will be filled in later)
                               opt=list(par=NA,fval=NA,conv=NA),
                               glmod$reTrms,fr=glmod$fr)
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
    n <- nrow(X <- getME(object, "X"))
    link <- if (isGLMM(object)) "response"

    ## predictions, conditioned as specified, on link scale
    ## previously: do **NOT** use na.action as specified here (inherit
    ##     from object instead, for consistency)
    ## now: use na.omit, because we have to match up
    ##    with whatever is done in mkNewReTrms
    etapred <- predict(object, newdata=newdata, re.form=re.form,
                       type="link", na.action=na.omit)

    ## now add random components:
    ##  only the ones we did *not* condition on

    ## compre.form <- noLHSform(formula(object))
    ## construct RE formula ONLY: leave out fixed terms,
    ##   which might have loose terms like offsets in them ...
    fb <- findbars(formula(object))
    pfun <- function(x) paste("(",paste(deparse(x),collapse=" "),")")
    compReForm <- reformulate(sapply(fb,pfun))


    if (!noReForm(re.form)) {
        rr <- re.form[[length(re.form)]] ## RHS of formula
        ftemplate <- substitute(.~.-XX, list(XX=rr))
        compReForm <- update.formula(compReForm,ftemplate)
    }

    ## (1) random effect(s)
    if (!is.null(findbars(compReForm))) {
        newRE <- mkNewReTrms(object,newdata,compReForm,
                             na.action=na.action,
                             allow.new.levels=allow.new.levels)
        U <- t(getME(object, "Lambdat") %*% newRE$Zt)
        u <- rnorm(ncol(U)*nsim)
        sim.reff <- ## UNSCALED random-effects contribution:
            as(U %*% matrix(u, ncol = nsim), "matrix")
    } else sim.reff <- 0
    if (isLMM(object)) {
        ## result will be matrix  n x nsim :
        val <- etapred + sigma * (sim.reff +
                                  ## residual contribution:
                                  matrix(rnorm(n * nsim), ncol = nsim))
    } else if (isGLMM(object)) {
        ## GLMM
        ## n.b. DON'T scale random-effects (???)
        etasim <- etapred+sim.reff
        family <- object@resp$family
        musim <- family$linkinv(etasim)
        ntot <- length(musim) ## FIXME: or could be dims["n"]?
        ## FIXME: is it possible to leverage family$simulate ... ???
        ##
        if (is.null(sfun <- simfunList[[family$family]]) &&
            is.null(family$simulate))
            stop("simulation not implemented for family",
                 family$family)
        ## don't rely on automatic recycling
        musim <- rep(musim,length.out=n*nsim)
        val <- sfun(object,
                    nsim=1,
                    ftd=musim)
        ## split results into nsims: need special case for binomial matrix/factor responses
        if (family$family=="binomial" && is.matrix(r <- model.response(object@frame))) {
            val <- lapply(split(val[[1]], gl(nsim, n, 2 * nsim * n)), matrix,
                          ncol = 2, dimnames = list(NULL, colnames(r)))
        } else if (family$family=="binomial" && is.factor(val[[1]])) {
            val <- split(val[[1]],gl(nsim,n))
        } else val <- split(val,gl(nsim,n))
    } else {
        stop("simulate method for NLMMs not yet implemented")
    }
    ## from src/library/stats/R/lm.R
    if(!is.list(val)) {
        dim(val) <- c(n, nsim)
        val <- as.data.frame(val)
    }  else class(val) <- "data.frame"
    names(val) <- paste("sim", seq_len(nsim), sep="_")
    ## have not yet filled in NAs, so need to use names of fitted
    ## object NOT including values with NAs
    f <- fitted(object)
    nm <- names(f)[!is.na(f)]  
    if (length(nm)==0) nm <- seq(n)
    row.names(val) <- nm

    fit.na.action <- attr(model.frame(object),"na.action")

    if (!missing(na.action) &&  !is.null(fit.na.action)) {
        ## retrieve name of na.action type ("omit", "exclude", "pass")
        class.na.action <- class(attr(na.action(NA),"na.action"))
        if (class.na.action != class(fit.na.action)) {
            ## hack to override action where explicitly specified
            class(fit.na.action) <- class.na.action
        }
    }

    if (is.matrix(val[[1]])) {
        ## have to handle binomial response matrices differently --
        ## fill in NAs as appropriate in *both* columns
        val <- lapply(val,function(x) { apply(x,2,napredict,
                                       omit=fit.na.action) })
        ## have to put this back into a (weird) data frame again,
        ## carefully (should do the napredict stuff
        ## earlier, so we don't have to redo this transformation!)
        class(val) <- "data.frame"
    } else {
        val <- as.data.frame(lapply(val,napredict, omit=fit.na.action))
    }

    ## reconstruct names: first get rid of NAs, then refill them
    ## as appropriate based on fit.na.action (which may be different
    ## from the original model's na.action spec)
    if (length(nm2 <- names(napredict(na.omit(f),
                                        omit=fit.na.action)))>0)
        row.names(val) <- nm2

    ## as.data.frame(lapply(...)) blows away na.action attribute,
    ##  so we have to re-assign here
    attr(val,"na.action") <- fit.na.action

    attr(val, "seed") <- RNGstate
    val
}

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
##
##
## these can be incorporated by overwriting the simulate()
## components, or calling them
##
gaussian_simfun <- function(object, nsim, ftd=fitted(object)) {
    wts <- weights(object)
    if (any(wts != 1)) warning("ignoring prior weights")
    rnorm(nsim*length(ftd), ftd, sd=sigma(object))
}

binomial_simfun <- function(object, nsim, ftd=fitted(object)) {
    n <- length(ftd)
    ntot <- n*nsim
    wts <- weights(object)
    if (any(wts %% 1 != 0))
        stop("cannot simulate from non-integer prior.weights")
    ## Try to fathom out if the original data were
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

poisson_simfun <- function(object, nsim, ftd=fitted(object)) {
        ## A Poisson GLM has dispersion fixed at 1, so prior weights
        ## do not have a simple unambiguous interpretation:
        ## they might be frequency weights or indicate averages.
        wts <- weights(object)
        if (any(wts != 1)) warning("ignoring prior weights")
        rpois(nsim*length(ftd), ftd)
    }


##' FIXME: need a gamma.shape.merMod method in order for this to work.
##'        (see initial shot at gamma.shape.merMod below)
Gamma_simfun <- function(object, nsim, ftd=fitted(object)) {
    wts <- weights(object)
    if (any(wts != 1)) message("using weights as shape parameters")
    ## ftd <- fitted(object)
    shape <- MASS::gamma.shape(object)$alpha * wts
    rgamma(nsim*length(ftd), shape = shape, rate = shape/ftd)
}

gamma.shape.merMod <- function(object, ...) {
    if(!(family(object)$family == "Gamma"))
        stop("Can not fit gamma shape parameter ",
             "because Gamma family not used")

    y <- getME(object, "y")
    mu <- getME(object, "mu")
    w <- weights(object)
                                        # Sec 8.3.2 (MN)
    L <- w*(log(y/mu)-((y-mu)/mu))
    dev <- -2*sum(L)
                                        # Eqs. between 8.2 & 8.3 (MN)
    Dbar <- dev/length(y)
    alpha <- (6+2*Dbar)/(Dbar*(6+Dbar))
                                        # FIXME: obtain standard error
    res <- list(alpha = alpha, SE = NA)
    class(res) <- "gamma.shape"
    res
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
negative.binomial_simfun <- function (object, nsim, ftd=fitted(object))
{
    stop("not implemented yet")
    ## val <- rnbinom(nsim * length(ftd), mu=ftd, size=.Theta)
}

simfunList <- list(gaussian=gaussian_simfun,
                binomial=binomial_simfun,
                poisson=poisson_simfun,
                Gamma=Gamma_simfun,
                negative.binomial=negative.binomial_simfun)

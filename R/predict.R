##' test for no-random-effect specification: TRUE if NA or ~0, other
##' possibilities are NULL or a non-trivial formula
noReForm <- function(ReForm) {
    (!is.null(ReForm) && !is(ReForm,"formula") && is.na(ReForm)) ||
        (is(ReForm,"formula") && length(ReForm)==2 && identical(ReForm[[2]],0))
}
## noReForm(NA)   ## TRUE
## noReForm(~0)   ## TRUE
## noReForm(~y+x)
## noReForm(NULL)
## noReForm(~0+x)

##' Force new parameters into a merMod object
##' (should this be an S3 method?)
##' @param params a list of the form specified by \code{start} in
##' \code{\link{lmer}} or \code{\link{glmer}} (i.e. a list containing
##' theta and/or beta; maybe eventually further parameters ...
##' What updating do we have to do in order to make the resulting object
##' consistent/safe?
##' TODO: What kind of checking of input do we have to do?  (Abstract from
##' lmer/glmer code ...)
##' TODO: make sure this gets updated when the parameter structure changes
##' from (theta, beta) to alpha=(theta, beta, phi)
setParams <- function(object,params,copy=TRUE) {
    if (copy) {
        newObj <- object
        ## make a copy of the reference class slots to
        ##  decouple them from the original object
        newObj@pp <- newObj@pp$copy()
        newObj@resp <- newObj@resp$copy()
        if (!is.null(beta <- params$beta)) newObj@pp$setBeta0(beta)
        if (!is.null(theta <- params$theta)) {
            ## where does theta live and how do I set it?
            ## (1) in .@theta
            ## (2) in .@pp$theta
            newObj@theta <- theta
            newObj@pp$setTheta(theta)
        }
        return(newObj)
    } else stop("modification in place (copy=FALSE) not yet implemented")
}

##' @examples
##' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
##' fm1M <- setParams(fm1,list(theta=rep(1,3)))
##' getME(fm1M,"theta")
##' getME(fm1,"theta") ## check that original didn't get messed up
##' ## check that @resp and @pp are the only reference class slots ...
##' sapply(slotNames(fm1),function(x) class(slot(fm1,x)))

##' Make new random effect terms from specified object and new data,
##' possibly omitting some random effect terms
##' @param object fitted model object
##' @param newdata (optional) data frame containing new data
##' @param ReForm formula specifying random effect terms to include (NULL=all, ~0)
##' @param na.action
##'
mkNewReTrms <- function(object,newdata,ReForm=NULL, na.action=na.pass,
                        allow.new.levels=FALSE) {
    ## construct (fixed) model frame in order to find out whether there are
    ## missing data/what to do about them
    mfnew <- if (is.null(newdata)) {
        model.frame(object)  ## FIXME: check ...
    } else {
        Terms <- terms(object,fixed.only=TRUE)
        model.frame(delete.response(Terms),newdata, na.action=na.action)
    }
    if (is(ReForm,"formula")) {
        ## DROP values with NAs in fixed effects
        if (length(fit.na.action <- attr(mfnew,"na.action"))>0) {
            newdata <- newdata[-fit.na.action,]
        }
        re <- ranef(object)
        rfd <- if(is.null(newdata)) object@frame else newdata 
        ReTrms <- mkReTrms(findbars(ReForm[[2]]), rfd)
        if (!allow.new.levels &&
            any(sapply(ReTrms$flist,function(x) any(is.na(x)))))
            stop("NAs are not allowed in prediction data",
                 " for grouping variables unless allow.new.levels is TRUE")
        unames <- unique(sort(names(ReTrms$cnms)))  ## FIXME: same as names(ReTrms$flist) ?
        ## convert numeric grouping variables to factors as necessary
        for (i in all.vars(ReForm[[2]])) {
            rfd[[i]] <- factor(rfd[[i]])
        }
        Rfacs <- setNames(lapply(unames,
                                 function(x) factor(eval(parse(text=x),envir=rfd))),
                          unames)
        new_levels <- lapply(Rfacs,function(x) levels(droplevels(factor(x))))
        ## FIXME: should this be unique(as.character(x)) instead?
        ##   (i.e., what is the proper way to protect against non-factors?)
        levelfun <- function(x,n) {
            ## find and deal with new levels
            if (any(!new_levels[[n]] %in% rownames(x))) {
                if (!allow.new.levels) stop("new levels detected in newdata")
                ## create an all-zero data frame corresponding to the new set of levels ...
                newx <- as.data.frame(matrix(0,nrow=length(new_levels[[n]]),
                                             ncol=ncol(x),
                                             dimnames=list(new_levels[[n]],
                                                 names(x))))
                ## then paste in the matching RE values from the original fit/set of levels
                newx[rownames(x),] <- x
                x <- newx
            }
            ## find and deal with missing old levels
            if (any(!rownames(x) %in% new_levels[[n]])) {
                x <- x[rownames(x) %in% new_levels[[n]],,drop=FALSE]
            }
            x
        }
        ## fill in/delete levels as appropriate
        re_x <- mapply(levelfun,re,names(re),SIMPLIFY=FALSE)
        re_new <- list()
        if (any(!names(ReTrms$cnms) %in% names(re)))
            stop("grouping factors specified in ReForm that were not present in original model")
        ## pick out random effects values that correspond to
        ##  random effects incorporated in ReForm ...
        for (i in seq_along(ReTrms$cnms)) {
            rname <- names(ReTrms$cnms)[i]
            if (any(!ReTrms$cnms[[rname]] %in% names(re[[rname]])))
                stop("random effects specified in ReForm",
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
##' @param ReForm formula for random effects to condition on.  If \code{NULL},
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
##' str(p1 <- predict(gm1,ReForm=NA))  # fitted values, unconditional (level-0)
##' newdata <- with(cbpp, expand.grid(period=unique(period), herd=unique(herd)))
##' str(p2 <- predict(gm1,newdata))    # new data, all RE
##' str(p3 <- predict(gm1,newdata,ReForm=NA)) # new data, level-0
##' str(p4 <- predict(gm1,newdata,ReForm=~(1|herd))) # explicitly specify RE
##' @method predict merMod
##' @export
predict.merMod <- function(object, newdata=NULL, newparams=NULL, newX=NULL,
                           ReForm=NULL,
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
    
    if (length(list(...)>0)) warning("unused arguments ignored")
    type <- match.arg(type)
    if (!is.null(terms))
        stop("terms functionality for predict not yet implemented")
    if (!is.null(newparams)) {
        object <- setParams(object,newparams)
    }
    if (is.null(newdata) && is.null(ReForm) && is.null(newparams)) {
        ## raw predict() call, just return fitted values
        ##   (inverse-link if appropriate)
        if (isLMM(object) || isNLMM(object)) {
            pred <- fitted(object)
        } else { ## inverse-link
            pred <-  switch(type,response=object@resp$mu, ## == fitted(object),
                            link=object@resp$eta)
        }
        if (!is.null(fit.na.action <- attr(model.frame(object),"na.action"))) {
            pred <- napredict(fit.na.action,pred)
        }
        return(pred)
    } else { ## newdata and/or ReForm and/or newparams specified
        X_orig <- getME(object, "X")
        if (is.null(newdata)) {
            X <- X_orig
            fit.na.action <- attr(object@frame,"na.action")  ## original NA action
        } else {
            ## evaluate new fixed effect
            RHS <- formula(object,fixed.only=TRUE)[-2]
            Terms <- terms(object,fixed.only=TRUE)
            X <- model.matrix(RHS, mfnew <- model.frame(delete.response(Terms),
                                                        newdata,
                                                        na.action=na.action),
                              contrasts.arg=attr(X_orig,"contrasts"))
            fit.na.action <- attr(mfnew,"na.action")
        }
        pred <- drop(X %*% fixef(object))
        ## modified from predict.glm ...
        offset <- rep(0, nrow(X))
        tt <- terms(object)
        ## FIXME:: need to unname()  ?
        if (!is.null(off.num <- attr(tt, "offset"))) {
            for (i in off.num)
                offset <- offset + eval(attr(tt,"variables")[[i + 1]], newdata)
        }
        ## FIXME: is this redundant??
        if (!is.null(frOffset <- attr(object@frame,"offset")))
            offset <- offset + eval(frOffset, newdata)
        pred <- pred+offset
        if (!noReForm(ReForm)) {
            if (is.null(ReForm)) {
                ## original formula, minus response
                ReForm <- formula(object)[-2]  
            }
            newRE <- mkNewReTrms(object,newdata,ReForm,na.action=na.action,
                                 allow.new.levels=allow.new.levels)
            pred <- pred + drop(as.matrix(newRE$b %*% newRE$Zt))
        }            
        if (isGLMM(object) && type=="response") {
            pred <- object@resp$family$linkinv(pred)
        }
        ## fill in NAs as appropriate
        if (!is.null(fit.na.action)) {
            pred <- napredict(fit.na.action,pred)
        }
        return(pred)
    } ## newdata/newparams/ReForm
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
##' @param use.u (logical) if \code{TRUE}, generate a simulation conditional on the current
##' random-effects estimates; if \code{FALSE} generate new Normally distributed random-effects values
##' @param ReForm (formula): specifies which random-effects estimates to condition on; if specified, overrides \code{use.u}
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
simulate.merMod <- function(object, nsim = 1, seed = NULL, use.u = FALSE,
                            ReForm=NA, newdata=NULL, newparams=NULL,
                            allow.new.levels=FALSE, ...) {
    stopifnot((nsim <- as.integer(nsim[1])) > 0,
              is(object, "merMod"))
    if (!is.null(newparams)) {
        object <- setParams(object,newparams)
    }
    if (!missing(use.u)) {
        if (!missing(ReForm)) {
            stop("should specify only one of ",sQuote("use.u"),
                 " and ",sQuote("ReForm"))
        }
        ReForm <- if (use.u) NULL else ~0
    }
    if (is.null(ReForm)) {
        ## original formula, minus response
        ReForm <- formula(object)[-2]  
    }
    if(!is.null(seed)) set.seed(seed)
    if(!exists(".Random.seed", envir = .GlobalEnv))
        runif(1) # initialize the RNG if necessary
    RNGstate <- .Random.seed

    sigma <- sigma(object)
    n <- nrow(X <- getME(object, "X"))
    if (is.null(nm <- names(fitted(object)))) nm <- seq(n)
    link <- if (isGLMM(object)) "response"

    ## predictions, conditioned as specified, on link scale
    etapred <- predict(object, newdata=newdata, ReForm=ReForm,
                       type="link")

    ## now add random components:
    ##  only the ones we did *not* condition on

    compReForm <- formula(object)[-2]

    if (!noReForm(ReForm)) {
        rr <- ReForm[[length(ReForm)]] ## RHS of formula
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
        val <- switch(family$family,
                      poisson=rpois(ntot,lambda=musim),
                      binomial={
                          w <- weights(object)
                          Y <- rbinom(ntot,prob=musim,size=w)
                          resp <- model.response(object@frame)
                          if (!is.matrix(resp)) {  ## bernoulli, or weights specified
                              if (is.factor(resp)) {
                                  if (any(weights(object)!=1)) stop("non-uniform weights with factor response??")
                                  f <- factor(levels(resp)[Y+1],levels=levels(resp))
                                  split(f, rep(seq_len(nsim), each = n))
                              } else {
                                  Y/w
                              }
                          } else {
                              ## FIXME: should "N-size" (column 2) be named?
                              ## copying structures from stats/R/family.R
                              nresp <- nrow(resp)
                              YY <- cbind(Y, w - Y)
                              yy <- lapply(split(YY,gl(nsim,nresp,2*nsim*nresp)),
                                           matrix, ncol=2,
                                           dimnames=list(NULL,colnames(resp)))
                              names(yy) <- paste("sim",seq_along(yy),sep="_")
                              yy
                          }
                      },
                      stop("simulation not implemented for family",
                           family$family))
    } else {
        stop("simulate method for NLMMs not yet implemented")
    }
    ## from src/library/stats/R/lm.R
    if(!is.list(val)) {
        dim(val) <- c(n, nsim)
        val <- as.data.frame(val)
    }
    else
        class(val) <- "data.frame"
    names(val) <- paste("sim", seq_len(nsim), sep="_")
    row.names(val) <- nm
    attr(val, "seed") <- RNGstate
    val
}

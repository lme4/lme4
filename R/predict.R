##' \code{\link{predict}} method for \code{\linkS4class{merMod}} objects
##'
##' @title Predictions from a model at new data values
##' @param object a fitted model object
##' @param newdata data frame for which to evaluate predictions
##' @param REform formula for random effects to include.  If NULL,
##'    include all random effects; if NA, include no random effects
##' @param terms a \code{\link{terms}} object - not used at present
##' @param type character string - either \code{"link"}, the default,
##'    or \code{"response"} indicating the type of prediction object returned
##' @param allow.new.levels (logical) if FALSE (default), then any new levels
##'    (or NA values) detected in \code{newdata} will trigger an error; if TRUE, then
##'    the prediction will use the unconditional (population-level)
##'    values for data with previously unobserved levels (or NAs)
##' @param na.action function determining what should be done with missing values for fixed effects in \code{newdata}. The default is to predict \code{NA}: see \code{\link{na.pass}}.
##' @param ... optional additional parameters.  None are used at present.
##' @return a numeric vector of predicted values
##' @note There is no option for computing standard errors of predictions because it is difficult to define an efficient method that incorporates uncertainty in the variance parameters; we recommend \code{\link{bootMer}} for this task.
##' @examples
##' (gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 |herd), cbpp, binomial))
##' str(p0 <- predict(gm1))            # fitted values
##' str(p1 <- predict(gm1,REform=NA))  # fitted values, unconditional (level-0)
##' newdata <- with(cbpp, expand.grid(period=unique(period), herd=unique(herd)))
##' str(p2 <- predict(gm1,newdata))    # new data, all RE
##' str(p3 <- predict(gm1,newdata,REform=NA)) # new data, level-0
##' str(p4 <- predict(gm1,newdata,REform=~(1|herd))) # explicitly specify RE
##' @method predict merMod
##' @export
predict.merMod <- function(object, newdata=NULL, newparams=NULL, REform=NULL,
                           terms=NULL, type=c("link","response"),
                           allow.new.levels=FALSE, na.action=na.pass, ...) {
    ## FIXME: appropriate names for result vector?
    ## FIXME: make sure behaviour is entirely well-defined for NA in grouping factors

    if (length(list(...)>0)) warning("unused arguments ignored")
    if (isLMM(object) && !missing(type)) warning("type argument ignored for linear mixed models")
    fit.na.action <- attr(object@frame,"na.action")
    type <- match.arg(type)
    if (!is.null(terms)) stop("terms functionality for predict not yet implemented")
    ## FIXME/WARNING: how do we/can we do this in an eval-safe way???
    form_orig <- formula(object)
    if (is.null(newdata) && is.null(REform)) {
        ## raw predict() call, just return fitted values (inverse-link if appropriate)
        if (isLMM(object) || isNLMM(object)) {
            pred <- fitted(object)
        } else { ## inverse-link
            pred <-  switch(type,response=object@resp$mu, ## fitted(object),
                            link=object@resp$eta)  ## fixme: getME() ?
        }
        if (!is.null(fit.na.action)) {
            pred <- napredict(fit.na.action,pred)
        }
        return(pred)
   
      } else { ## newdata and/or REform specified
        X_orig <- getME(object, "X")
        if (is.null(newdata)) {
            X <- X_orig
        } else {
            ## evaluate new fixed effect
            RHS <- formula(object,fixed.only=TRUE)[-2]
            Terms <- terms(object,fixed.only=TRUE)
            X <- model.matrix(RHS, mfnew <- model.frame(delete.response(Terms),
                                                        newdata, na.action=na.action),
                              contrasts.arg=attr(X_orig,"contrasts"))
        }
        pred <- drop(X %*% fixef(object))
        ## modified from predict.glm ...
        offset <- rep(0, nrow(X))
        tt <- terms(object)
        ## FIXME:: need to unname()  ?
        if (!is.null(off.num <- attr(tt, "offset"))) {
            for (i in off.num) offset <- offset + eval(attr(tt,"variables")[[i + 1]], newdata)
        }
        ## FIXME: is this redundant??
        if (!is.null(frOffset <- attr(object@frame,"offset")))
            offset <- offset + eval(frOffset, newdata)
        pred <- pred+offset
        if (is.null(REform)) {
            REform <- form_orig[-2]
        }
        ## FIXME: ??? can't apply is.na() to a 'language' object?
        ##  what's the appropriate test?
        if (is.language(REform)) {
            na.action.name <- deparse(match.call()$na.action) ## ugh
            if (!is.null(newdata) && na.action.name %in% c("na.exclude","na.omit")) {
                ## strip NAs from data for random-effects matrix construction
                if (length(nadrop <- attr(mfnew,"na.action"))>0) {
                    newdata <- newdata[-nadrop,]
                }
            }
            re <- ranef(object)
            ## ok? -- newdata used even though it was just tested for null
            if(is.null(newdata)) rfd <- object@frame else rfd <- newdata # get data for REform
            ReTrms <- mkReTrms(findbars(REform[[2]]), rfd)
            if (!allow.new.levels && any(sapply(ReTrms$flist,function(x) any(is.na(x)))))
                stop("NAs are not allowed in prediction data for grouping variables unless allow.new.levels is TRUE")
            unames <- unique(sort(names(ReTrms$cnms)))  ## FIXME: same as names(ReTrms$flist) ?
            ## convert numeric grouping variables to factors as necessary
            for (i in all.vars(REform[[2]])) {
                newdata[[i]] <- factor(newdata[[i]])
            }
            Rfacs <- setNames(lapply(unames,function(x) factor(eval(parse(text=x),envir=newdata))),
                              unames)
            new_levels <- lapply(Rfacs,function(x) levels(droplevels(factor(x))))
            ## FIXME: should this be unique(as.character(x)) instead?
            ##   (i.e., what is the proper way to protect against non-factors?)
            levelfun <- function(x,n) {
                ## find and deal with new levels
                if (any(!new_levels[[n]] %in% rownames(x))) {
                    if (!allow.new.levels) stop("new levels detected in newdata")
                    ## create an all-zero data frame corresponding to the new set of levels ...
                    newx <- as.data.frame(matrix(0,nrow=length(new_levels[[n]]),ncol=ncol(x),
                                                 dimnames=list(new_levels[[n]],names(x))))
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
                stop("grouping factors specified in REform that were not present in original model")
            ## pick out random effects values that correspond to
            ##  random effects incorporated in REform ...
            for (i in seq_along(ReTrms$cnms)) {
                rname <- names(ReTrms$cnms)[i]
                if (any(!ReTrms$cnms[[rname]] %in% names(re[[rname]])))
                    stop("random effects specified in REform that were not present in original model")
                re_new[[i]] <- re_x[[rname]][,ReTrms$cnms[[rname]]]
            }
            re_newvec <- unlist(lapply(re_new,t))  ## must TRANSPOSE RE matrices before unlisting
            if(!is.null(newdata)) pred <- pred + drop(as.matrix(re_newvec %*% ReTrms$Zt))
        } ## predictions with REform!=NA
        if (isGLMM(object) && type=="response") {
            pred <- object@resp$family$linkinv(pred)
        }
        ## fill in NAs as appropriate
        if (is.null(newdata) && !is.null(fit.na.action)) {
            pred <- napredict(fit.na.action,pred)
        } else {
            pred <- napredict(na.action,pred)
        }
        return(pred)
    }
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
                            REform=NA, ...) {
    stopifnot((nsim <- as.integer(nsim[1])) > 0,
	      is(object, "merMod"))
    if (!missing(use.u) && !missing(REform)) stop("should specify only one of ",sQuote("use.u"),
                                                  " and ",sQuote("REform"))
    if (!missing(use.u)) REform <- if (use.u) NULL else NA
    if(!is.null(seed)) set.seed(seed)
    if(!exists(".Random.seed", envir = .GlobalEnv))
	runif(1) # initialize the RNG if necessary
    RNGstate <- .Random.seed

    sigma <- sigma(object)
    n <- nrow(X <- getME(object, "X"))
    if (is.null(nm <- names(fitted(object)))) nm <- seq(n)


    pp <- predict(object, newdata=newdata, REform
    # fixed-effect contribution
    etasim.fix <- as.vector(X %*% getME(object, "beta"))
    if (length(offset <- getME(object,"offset"))>0) {
      etasim.fix <- etasim.fix+offset
    }
    U <- getME(object, "Z") %*% getME(object, "Lambda")
    u <- if (use.u) {
        rep(getME(object, "u"), nsim)/sigma  ## ??? u is 'spherized' but not scaled ???
    } else {
        rnorm(ncol(U)*nsim)
    }
    etasim.reff <- ## UNSCALED random-effects contribution:
        as(U %*% matrix(u, ncol = nsim), "matrix")
    if (is(object@resp,"lmerResp")) {
      ## result will be matrix  n x nsim :
      val <- etasim.fix + sigma * (etasim.reff +
        ## residual contribution:
        matrix(rnorm(n * nsim), ncol = nsim))
    } else if (is(object@resp,"glmResp")) {
      ## GLMM
      ## n.b. DON'T scale random-effects (???)
      	      etasim <- etasim.fix+etasim.reff
              ## FIXME:: try to avoid @call ...
	      family <- object@call$family
	      if(is.symbol(family)) family <- as.character(family)
	      if(is.character(family))
		  family <- get(family, mode = "function", envir = parent.frame(2))
	      if(is.function(family)) family <- family()
              if(is.language(family)) family <- eval(family)
	      if(is.null(family$family)) stop("'family' not recognized")
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

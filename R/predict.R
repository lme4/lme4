##' \code{\link{predict}} method for \code{\linkS4class{merMod}} objects
##'
##' @title Predictions from a model at new data values
##' @param object a fitted model object
##' @param newdata data frame for which to evaluate predictions
##' @param REform formula for random effects to include.  If NULL,
##'    include all random effects; if NA, include no random effects
##' @param terms a \code{\link{terms}} object - not used at present
##' @param type character string - either \code{"link"}, the default,
##'    or \code{"response"} indicating the type of prediction object returned.
##' @param allow.new.levels (logical) if FALSE, then any new levels
##'    detected in \code{newdata} will trigger an error; if TRUE, then
##'    the prediction will use the unconditional (population-level)
##'    values for data with previously unobserved levels
##' @param ... optional additional parameters.  None are used at present.
##' @return a numeric vector of predicted values
##' @note explain why there is no option for computing standard errors of predictions!
##' @note offsets not yet handled
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
predict.merMod <- function(object, newdata=NULL, REform=NULL,
                           terms=NULL, type=c("link","response"),
                           allow.new.levels=FALSE, ...) {
    ## FIXME: appropriate names for result vector?
    if (any(getME(object,"offset")!=0)) stop("offsets not handled yet")  ## FIXME
    type <- match.arg(type)
    if (!is.null(terms)) stop("terms functionality for predict not yet implemented")
    X_orig <- getME(object, "X")
    ## FIXME/WARNING: how do we do this in an eval-safe way???
    form_orig <- eval(object@call$formula,parent.frame())
    if (is.null(newdata) && is.null(REform)) {
        if (is(object@resp,"lmerResp")) return(fitted(object))
        return(switch(type,response=object@resp$mu, ## fitted(object),
                      link=object@resp$eta))  ## fixme: getME() ?
    } else {
        if (is.null(newdata)) {
            X <- X_orig
        } else {
            form <- form_orig
            form[[3]] <- if(is.null(nb <- nobars(form[[3]]))) 1 else nb
            RHS <- form[-2]
            X <- model.matrix(RHS, newdata, contrasts.arg=attr(X_orig,"contrasts"))
        }
        pred <- drop(X %*% fixef(object))
        if (is.null(REform)) {
            REform <- form_orig[-2]
        }
        ## FIXME: ??? can't apply is.na() to a 'language' object?
        ##  what's the appropriate test?
        if (is.language(REform)) {
            re <- ranef(object)
            ## 
            ReTrms <- mkReTrms(findbars(REform[[2]]),newdata)
            new_levels <- lapply(newdata[unique(sort(names(ReTrms$cnms)))],levels)
            re_x <- mapply(function(x,n) {
                if (any(!new_levels[[n]] %in% rownames(x))) {
                    if (!allow.new.levels) stop("new levels detected in newdata")
                    ## create an all-zero data frame corresponding to the new set of levels ...
                    newx <- as.data.frame(matrix(0,nrow=length(new_levels[[n]]),ncol=ncol(x),
                                                 dimnames=list(new_levels[[n]],names(x))))
                    ## then paste in the matching RE values from the original fit/set of levels
                    newx[rownames(x),] <- x
                    x <- newx
                }
                x
            },
                           re,names(re),SIMPLIFY=FALSE)
            ## separate random effects from orig model into individual columns
            re_List <- do.call(c,lapply(re_x,as.list))
            re_names <- names(re_List)
            z_names <- mapply(paste,names(ReTrms$cnms),ReTrms$cnms,MoreArgs=list(sep="."))
            ## pick out random effects values that correspond to
            ##  random effects incorporated in REform ...
            ## FIXME: more tests for possible things going wrong here?
            m <- match(z_names,re_names)
            if (any(is.na(m)))
                stop("random effects specified in REform that were not present in original model")
            re_new <- unlist(re_List[m])
            pred <- pred + drop(as.matrix(re_new %*% ReTrms$Zt))
        } ## REform provided
    } ## predictions with new data or new REform
    ## FIXME: would like to have an isGLMM() accessor for merMod objects?
    if (is(object@resp,"glmResp") && type=="response") {
        pred <- object@resp$family$linkinv(pred)
    }
    return(pred)
}
  
    

##' @export
lmer1 <- function(formula, data=NULL, REML = TRUE, sparseX = FALSE,
                 control = list(), start = NULL,
                 verbose = 0L, subset, weights, na.action, offset,
                 contrasts = NULL, devFunOnly=FALSE,
                 optimizer="Nelder_Mead", ...)
{
    verbose <- as.integer(verbose)
    restart <- TRUE ## FIXME; set default elsewhere?
    if (!is.null(control$restart)) {
        restart <- control$restart
        control$restart <- NULL
    }

    mf <- mc <- match.call()
    checkArgs("lmer",sparseX,...)
    if (!is.null(list(...)[["family"]])) {
        ## lmer(...,family=...); warning issued within checkArgs
        mc[[1]] <- as.name("glmer")
        return(eval(mc, parent.frame()) )
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
                                        # fixed-effects model matrix X
    fixedform <- formula
    fixedform[[3]] <- if(is.null(nb <- nobars(fixedform[[3]]))) 1 else nb
    mf$formula <- fixedform
                                        # re-evaluate model frame to extract predvars component
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"),"predvars.fixed") <- attr(attr(fixedfr,"terms"),"predvars")
    X <- model.matrix(fixedform, fixedfr, contrasts)
    p <- ncol(X)
    if ((rankX <- rankMatrix(X)) < p)
        stop(gettextf("rank of X = %d < ncol(X) = %d", rankX, p))
    if (!length(bb <- findbars(formula[[3]])))
        stop("No random effects terms specified in formula")
    createRSC(bb, X, fr)
}

#' Create an RSC object from the formula and the model.frame
#' 
#' @param bb list of random-effects terms (ie. result of findbars)
#' @param X the fixed-effects model matrix
#' @param fr the model frame for the mixed-effects model
#' @return an RSC object
#' @examples
#' fm1 <- lmer1(diameter ~ (1|plate) + (1|sample), Penicillin)
#' str(fm1)
#' with(Penicillin, RSCupdate(fm1, diameter))
#' @export
createRSC <- function(bb, X, fr) {
    stopifnot(is.list(bb), all(sapply(bb, is.language)),
              is.matrix(X), (p <- ncol(X)) > 0L,
              is.data.frame(fr))
    n <- nrow(X)
    names(bb) <- unlist(lapply(bb, function(x) deparse(x[[3]])))
    rhs <- lapply(bb, function(x) factor(eval(x[[3]], fr))) # factors
    nl  <- sapply(rhs, function(x) length(levels(x)))       # number of levels
    if (any(diff(nl) > 0)) {            # order terms by decreasing nl
        ord <- rev(order(nl))
        bb  <- bb[ord]
        rhs <- rhs[ord]
        nl  <- nl[ord]
    }
    lhs <- lapply(bb, function(x)       # model matrices for r.e. terms
                  model.matrix(eval(substitute(~foo, list(foo=x[[2]]))), fr))
    x <- do.call(rbind, lapply(c(lhs, list(X)), t))
    dimnames(x) <- NULL
    nc  <- sapply(lhs, ncol)
    k   <- sum(nc)                      # number of non-zeros per row of Z
    q   <- sum(nc * nl)                 # total number of random effects
    if (any(dd <- duplicated(names(rhs)))) { # repeated grouping factors
        stop("Code not yet written")
    }
    theta <- unlist(lapply(nc, function(n) {
        if (n == 1L) return(1)
        id <- diag(nrow=n)
        id[lower.tri(id,diag=TRUE)]
    }))
    names(theta) <- NULL
    lower <- ifelse(theta == 0, -Inf, 0)
    uboff <- cumsum(c(0L, nl*nc, rep.int(1L, p))) # offsets into ubeta for rows of x
    ii <- do.call(rbind,c(lapply(rhs, as.integer),
                          list(matrix(1L,nrow=p,ncol=n)))) +
                              (uboff[seq_len(k+p)] - 1L)
    A <- tcrossprod(sparseMatrix(i=as.vector(ii),
                                 j=rep(0:(nrow(X)-1L), each=nrow(x)),
                                 x=as.vector(x),
                                 index1=FALSE)) +
                                     Diagonal(x=rep.int(c(1,0), c(q,p)))
    ## Cholesky results are not saved in this function but are cached in A
    ## FIXME: need to work out the permutation keeping Z and X parts distinct.
    ## For now set perm=FALSE
    Cholesky(A, perm=FALSE, LDL=FALSE) 
    new("RSC", x=x, i=ii[1:k,,drop=FALSE], theta=theta, lower=lower, A=A,
        ubeta=numeric(q+p), uboff=uboff)
}

# Update for the penalized least squares problem
#' @importFrom stats update
#' @S3method update RSC
update.RSC <- function(object, ...) {
    dots <- list(...)
    weights <- dots$weights
    resid <- dots$resid
    if (is.null(resid)) stop("resid must be specified in RSC update")
    resid <- as.vector(as.numeric(resid))
    if (is.null(weights)) weights <- rep.int(1, length(resid))
    .Call(lme4_RSCupdate, object, resid, weights)
}

##' @importFrom stats fitted
##' @S3method fitted RSC
fitted.RSC <- function(pred) .Call(lme4_RSCfitted, pred)

#' @export
RSCdevfun <- function(object, resp) {
    function(theta) {
        object@theta[] <- theta
        rr <- update(object, resid = resp)
        wrss <- sum((rr$linpred - resp)^2)
        ubeta <- rr$del_ubeta           # increment is coefficient vector for LMMs
        qpp <- length(ubeta)
        n <- length(resp)
        p <- nrow(object@x) - nrow(object@i)
        ussq <- sum(ubeta[seq_len(length(ubeta) - p)]^2)
        rr[["ldL2"]] + n * (1 + log(2 * pi * (wrss + ussq)/n))
    }
}

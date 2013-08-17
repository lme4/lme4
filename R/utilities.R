if(getRversion() < "2.15")
    paste0 <- function(...) paste(..., sep = '')

### Utilities for parsing and manipulating mixed-model formulas

##' From the result of \code{\link{findbars}} applied to a model formula and
##' and the evaluation frame, create the model matrix, etc. associated with
##' random-effects terms.  See the description of the returned value for a
##' detailed list.
##'
##' @title Create Z, Lambda, Lind, etc.
##' @param bars a list of parsed random-effects terms
##' @param fr a model frame in which to evaluate these terms
##' @return a list with components
##' \item{Zt}{transpose of the sparse model matrix for the random effects}
##' \item{Lambdat}{transpose of the sparse relative covariance factor}
##' \item{Lind}{an integer vector of indices determining the mapping of the
##'     elements of the \code{theta} to the \code{"x"} slot of \code{Lambdat}}
##' \item{theta}{initial values of the covariance parameters}
##' \item{lower}{lower bounds on the covariance parameters}
##' \item{flist}{list of grouping factors used in the random-effects terms}
##' \item{cnms}{a list of column names of the random effects according to
##'     the grouping factors}
##' @importFrom Matrix sparseMatrix rBind drop0
##' @importMethodsFrom Matrix coerce
##' @family utilities
##' @export
mkReTrms <- function(bars, fr) {
    if (!length(bars))
	stop("No random effects terms specified in formula")
    stopifnot(is.list(bars), all(sapply(bars, is.language)),
	      inherits(fr, "data.frame"))
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))

    ## auxiliary {named, for easier inspection}:
    mkBlist <- function(x) {
        frloc <- fr
        ## convert grouping variables to factors as necessary
        for (i in all.vars(x[[3]])) {
            frloc[[i]] <- factor(frloc[[i]])
        }
	ff <- eval(substitute(factor(fac), list(fac = x[[3]])), frloc)
	if (all(is.na(ff)))
	    stop("Invalid grouping factor specification, ",
		 deparse(x[[3]]))
	nl <- length(levels(ff))
	mm <- model.matrix(eval(substitute( ~ foo,
					   list(foo = x[[2]]))), fr)
	nc <- ncol(mm)
	nseq <- seq_len(nc)
	sm <- as(ff, "sparseMatrix")
	if (nc	> 1)
	    sm <- do.call(rBind, lapply(nseq, function(i) sm))
        ## hack for NA values contained in factor (FIXME: test elsewhere for consistency?)
	sm@x[] <- t(mm[!is.na(ff),])
	## When nc > 1 switch the order of the rows of sm
	## so the random effects for the same level of the
	## grouping factor are adjacent.
	if (nc > 1)
	    sm <- sm[as.vector(matrix(seq_len(nc * nl),
				      ncol = nl, byrow = TRUE)),]
	list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
    }
    blist <- lapply(bars, mkBlist)
    nl <- unlist(lapply(blist, "[[", "nl")) # no. of levels per term

    ## order terms stably by decreasing number of levels in the factor
    if (any(diff(nl)) > 0) {
	ord <- rev(order(nl))
	blist <- blist[ord]
	nl <- nl[ord]
    }
    Zt <- do.call(rBind, lapply(blist, "[[", "sm"))
    q <- nrow(Zt)

    ## Create and install Lambdat, Lind, etc.  This must be done after
    ## any potential reordering of the terms.
    cnms <- lapply(blist, "[[", "cnms")
    nc <- vapply(cnms, length, 1L)    # no. of columns per term
    nth <- as.integer((nc * (nc+1))/2)	# no. of parameters per term
    nb <- nc * nl                     # no. of random effects per term

    stopifnot(sum(nb) == q)
    boff <- cumsum(c(0L, nb))		# offsets into b
    thoff <- cumsum(c(0L, nth))		# offsets into theta
### FIXME: should this be done with cBind and avoid the transpose
### operator?  In other words should Lambdat be generated directly
### instead of generating Lambda first then transposing?
    Lambdat <-
        t(do.call(sparseMatrix,
                  do.call(rBind,
                          lapply(seq_along(blist), function(i)
                             {
                                 mm <- matrix(seq_len(nb[i]), ncol = nc[i],
                                              byrow = TRUE)
                                 dd <- diag(nc[i])
                                 ltri <- lower.tri(dd, diag = TRUE)
                                 ii <- row(dd)[ltri]
                                 jj <- col(dd)[ltri]
                                 dd[cbind(ii, jj)] <- seq_along(ii)
                                 data.frame(i = as.vector(mm[, ii]) + boff[i],
                                            j = as.vector(mm[, jj]) + boff[i],
                                            x = as.double(rep.int(seq_along(ii),
                                            rep.int(nl[i], length(ii))) +
                                            thoff[i]))
                             }))))
    thet <- numeric(sum(nth))
    ll <- list(Zt=Matrix::drop0(Zt), theta=thet, Lind=as.integer(Lambdat@x),
               Gp=unname(c(0L, cumsum(nb))))
    ## lower bounds on theta elements are 0 if on diagonal, else -Inf
    ll$lower <- -Inf * (thet + 1)
    ll$lower[unique(diag(Lambdat))] <- 0
    ll$theta[] <- is.finite(ll$lower) # initial values of theta are 0 off-diagonal, 1 on
    Lambdat@x[] <- ll$theta[ll$Lind]  # initialize elements of Lambdat
    ll$Lambdat <- Lambdat
                                        # massage the factor list
    fl <- lapply(blist, "[[", "ff")
                                        # check for repeated factors
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        fl <- fl[match(ufn, fnms)]
        asgn <- match(fnms, ufn)
    } else asgn <- seq_along(fl)
    names(fl) <- ufn
    fl <- do.call(data.frame, c(fl, check.names = FALSE))
    attr(fl, "assign") <- asgn
    ll$flist <- fl
    ll$cnms <- cnms
    ll
} ## {mkReTrms}

##' Create an lmerResp, glmResp or nlsResp instance
##'
##' @title Create an lmerResp, glmResp or nlsResp instance
##' @param fr a model frame
##' @param REML logical scalar, value of REML for an lmerResp instance
##' @param family the optional glm family (glmResp only)
##' @param nlenv the nonlinear model evaluation environment (nlsResp only)
##' @param nlmod the nonlinear model function (nlsResp only)
##' @param ... where to look for response information if \code{fr} is missing.
##'   Can contain a model response, \code{y}, offset, \code{offset}, and weights,
##'   \code{weights}.
##' @return an lmerResp or glmResp or nlsResp instance
##' @family utilities
##' @export
mkRespMod <- function(fr, REML=NULL, family = NULL, nlenv = NULL, nlmod = NULL, ...) {

  if(!missing(fr)){
    y <- model.response(fr)
    offset <- model.offset(fr)
    weights <- model.weights(fr)
    N <- n <- nrow(fr)
    etastart_update <- model.extract(fr, "etastart")
  }
  else{
    fr <- list(...)
    y <- fr$y
    N <- n <- if(is.matrix(y)) nrow(y) else length(y)
    offset <- fr$offset
    weights <- fr$weights
    etastart_update <- fr$etastart
  }

  ## FIXME: may need to add X, or pass it somehow, if we want to use glm.fit
  #y <- model.response(fr)
    if(length(dim(y)) == 1) {
	## avoid problems with 1D arrays, but keep names
	nm <- rownames(y)
	dim(y) <- NULL
	if(!is.null(nm)) names(y) <- nm
    }
    rho <- new.env()
    rho$y <- if (is.null(y)) numeric(0) else y
    if (!is.null(REML)) rho$REML <- REML
    rho$etastart <- fr$etastart
    rho$mustart <- fr$mustart
    #N <- n <- nrow(fr)
    if (!is.null(nlenv)) {
        stopifnot(is.language(nlmod),
                  is.environment(nlenv),
                  is.numeric(val <- eval(nlmod, nlenv)),
                  length(val) == n,
                  is.matrix(gr <- attr(val, "gradient")),
                  mode(gr) == "numeric",
                  nrow(gr) == n,
                  !is.null(pnames <- colnames(gr)))
        N <- length(gr)
        rho$mu <- as.vector(val)
        rho$sqrtXwt <- as.vector(gr)
        rho$gam <-
            unname(unlist(lapply(pnames,
                                 function(nm) get(nm, envir=nlenv))))
    }
    if (!is.null(offset)) {
        if (length(offset) == 1L) offset <- rep.int(offset, N)
        stopifnot(length(offset) == N)
        rho$offset <- unname(offset)
    } else rho$offset <- rep.int(0, N)
    if (!is.null(weights)) {
        stopifnot(length(weights) == n, all(weights >= 0))
        rho$weights <- unname(weights)
    } else rho$weights <- rep.int(1, n)
    if (is.null(family)) {
        if (is.null(nlenv)) return(do.call(lmerResp$new, as.list(rho)))
        return(do.call(nlsResp$new,
                       c(list(nlenv=nlenv,
                              nlmod=substitute(~foo, list(foo=nlmod)),
                              pnames=pnames), as.list(rho))))
    }
    stopifnot(inherits(family, "family"))
                              # need weights for initialize evaluation
    rho$nobs <- n
    eval(family$initialize, rho)
    family$initialize <- NULL     # remove clutter from str output
    ll <- as.list(rho)
    ans <- do.call(new, c(list(Class="glmResp", family=family),
                          ll[setdiff(names(ll), c("m", "nobs", "mustart"))]))
  ans$updateMu(if (!is.null(es <- etastart_update)) es else
                 family$linkfun(get("mustart", rho)))
    ans
}

##' From the right hand side of a formula for a mixed-effects model,
##' determine the pairs of expressions that are separated by the
##' vertical bar operator.  Also expand the slash operator in grouping
##' factor expressions.
##'
##' @title Determine random-effects expressions from a formula
##' @seealso \code{\link{formula}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
##' @param term a mixed-model formula
##' @return pairs of expressions that were separated by vertical bars
##' @section Note: This function is called recursively on individual
##' terms in the model, which is why the argument is called \code{term} and not
##' a name like \code{form}, indicating a formula.
##' @examples
##' findbars(f1 <- Reaction ~ Days + (Days|Subject))
##' ## => list( Days | Subject )
##' findbars(y ~ Days + (1|Subject) + (0+Days|Subject))
##' ## => list of length 2:  list ( 1 | Subject ,  0+Days|Subject)
##' findbars(~ 1 + (1|batch/cask))
##' ## => list of length 2:  list ( 1 | cask:batch ,  1 | batch)
##' identical(findbars(~ 1 + (Days || Subject)),
##'     findbars(~ 1 + (1|Subject) + (0+Days|Subject)))
##' \dontshow{
##' stopifnot(identical(findbars(f1),
##'                     list(expression(Days | Subject)[[1]])))
##' }
##' @family utilities
##' @keywords models utilities
##' @export
findbars <- function(term)
{
    ## Recursive function applied to individual terms
    fb <- function(term)
    {
        if (is.name(term) || !is.language(term)) return(NULL)
        if (term[[1]] == as.name("(")) return(fb(term[[2]]))
        stopifnot(is.call(term))
        if (term[[1]] == as.name('|')) return(term)
        if (length(term) == 2) return(fb(term[[2]]))
        c(fb(term[[2]]), fb(term[[3]]))
    }
    ## Expand any slashes in the grouping factors returned by fb
    expandSlash <- function(bb)
    {
        ## Create the interaction terms for nested effects
        makeInteraction <- function(x)
        {
            if (length(x) < 2) return(x)
            trm1 <- makeInteraction(x[[1]])
            trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
            list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
        }
        ## Return the list of '/'-separated terms
        slashTerms <- function(x)
        {
            if (!("/" %in% all.names(x))) return(x)
            if (x[[1]] != as.name("/"))
                stop("unparseable formula for grouping factor")
            list(slashTerms(x[[2]]), slashTerms(x[[3]]))
        }

        if (!is.list(bb)) return(expandSlash(list(bb)))
        ## lapply(unlist(... - unlist returns a flattened list
        unlist(lapply(bb, function(x) {
            if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
                return(lapply(unlist(makeInteraction(trms)),
                              function(trm) substitute(foo|bar,
                                                       list(foo = x[[2]],
                                                            bar = trm))))
            x
        }))
    }
    expandSlash(fb(term))
}

##' Remove the random-effects terms from a mixed-effects formula,
##' thereby producing the fixed-effects formula.
##'
##' @title Omit terms separated by vertical bars in a formula
##' @param term the right-hand side of a mixed-model formula
##' @return the fixed-effects part of the formula
##' @section Note: This function is called recursively on individual
##' terms in the model, which is why the argument is called \code{term} and not
##' a name like \code{form}, indicating a formula.
##' @examples
##' nobars(Reaction ~ Days + (Days|Subject)) ## => Reaction ~ Days
##' @seealso \code{\link{formula}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
##' @family utilities
##' @keywords models utilities
##' @export
nobars <- function(term)
{
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
	nb <- nobars(term[[2]])
	if (is.null(nb)) return(NULL)
	term[[2]] <- nb
	return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

##' Substitute the '+' function for the '|' function in a mixed-model
##' formula.  This provides a formula suitable for the current
##' model.frame function.
##'
##' @title "Sub[stitute] Bars"
##' @param term a mixed-model formula
##' @return the formula with all | operators replaced by +
##' @section Note: This function is called recursively on individual
##' terms in the model, which is why the argument is called \code{term} and not
##' a name like \code{form}, indicating a formula.
##' @examples
##' subbars(Reaction ~ Days + (Days|Subject)) ## => Reaction ~ Days + (Days + Subject)
##' @seealso \code{\link{formula}}, \code{\link{model.frame}}, \code{\link{model.matrix}}.
##' @family utilities
##' @keywords models utilities
##' @export
subbars <- function(term)
{
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
	term[[2]] <- subbars(term[[2]])
	return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
	term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

##' Does every level of f1 occur in conjunction with exactly one level
##' of f2? The function is based on converting a triplet sparse matrix
##' to a compressed column-oriented form in which the nesting can be
##' quickly evaluated.
##'
##' @title Is f1 nested within f2?
##'
##' @param f1 factor 1
##' @param f2 factor 2
##'
##' @return TRUE if factor 1 is nested within factor 2
##' @examples
##' with(Pastes, isNested(cask, batch))   ## => FALSE
##' with(Pastes, isNested(sample, batch))  ## => TRUE
##' @export
isNested <- function(f1, f2)
{
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    stopifnot(length(f1) == length(f2))
    k <- length(levels(f1))
    sm <- as(new("ngTMatrix",
		 i = as.integer(f2) - 1L,
		 j = as.integer(f1) - 1L,
		 Dim = c(length(levels(f2)), k)),
             "CsparseMatrix")
    all(sm@p[2:(k+1L)] - sm@p[1:k] <= 1L)
}

subnms <- function(form, nms) {
    ## Recursive function applied to individual terms
    sbnm <- function(term)
    {
        if (is.name(term))
            if (any(term == nms)) return(0) else return(term)
        switch(length(term),
               return(term),
           {
               term[[2]] <- sbnm(term[[2]])
               return(term)
           },
           {
               term[[2]] <- sbnm(term[[2]])
               term[[3]] <- sbnm(term[[3]])
               return(term)
           })
        NULL
    }
    sbnm(form)
}

## Check for a constant term (a literal 1) in an expression
##
## In the mixed-effects part of a nonlinear model formula, a constant
## term is not meaningful because every term must be relative to a
## nonlinear model parameter.  This function recursively checks the
## expressions in the formula for a a constant, calling stop() if
## such a term is encountered.
## @title Check for constant terms.
## @param expr an expression
## @return NULL.  The function is executed for its side effect.
chck1 <- function(expr) {
    if ((le <- length(expr)) == 1) {
        if (is.numeric(expr) && expr == 1)
            stop("1 is not meaningful in a nonlinear model formula")
        return()
    } else
    for (j in seq_len(le)[-1]) Recall(expr[[j]])
}

##' Check and manipulate the formula for a nonlinear model.
##'
##' The model formula for a nonlinear mixed-effects model is of the form
##' \code{resp ~ nlmod ~ mixed} where \code{"resp"} is an expression
##' (usually just a name) for the response, \code{nlmod} is the
##' call to the nonlinear model function, and \code{mixed} is the
##' mixed-effects formula defining the linear predictor for the
##' parameter matrix.  If the formula is to be used for optimizing
##' designs, the \code{"resp"} part can be omitted.
##'
##'
##' @title Manipulate a nonlinear model formula.
##' @param mc matched call from the calling function.  Should have arguments named
##' \describe{
##'     \item{formula}{a formula of the form \code{resp ~ nlmod ~ meform}
##'                    where \code{resp} is an expression for the response,
##'                    \code{nlmod} is the nonlinear model expression and
##'                    \code{meform} is the mixed-effects model formula. \code{resp}
##'                    can be omitted when, e.g., optimizing a design.}
##'     \item{data}{a data frame in which to evaluate the model function}
##'     \item{start}{either a numeric vector containing initial estimates for the
##'                  nonlinear model parameters or a list with components
##'          \describe{
##'              \item{nlpars}{the initial estimates of the nonlinear model parameters}
##'              \item{theta}{the initial estimates of the variance component parameters}
##'          }
##'     }
##' }
##' @return a list with components
##'  \item{"respMod"}{a response module of class \code{"\linkS4class{nlsResp}"}}
##'  \item{"frame"}{the model frame, including a terms attribute}
##'  \item{"X"}{the fixed-effects model matrix}
##'  \item{"reTrms"}{the random-effects terms object}
##' @export
##' @family utilities
nlformula <- function(mc) {
    start <- eval(mc$start, parent.frame(2L))
    if (is.numeric(start)) start <- list(nlpars = start)
    stopifnot(is.numeric(nlpars <- start$nlpars),
              all(sapply(nlpars, length) == 1L),
              length(pnames <- names(nlpars)) == length(nlpars),
              length(form <- as.formula(mc$formula)) == 3L,
              is(nlform <- eval(form[[2]]), "formula"),
              all(pnames %in%
                  (av <- all.vars(nlmod <- as.call(nlform[[lnl <- length(nlform)]])))))
    nlform[[lnl]] <- parse(text= paste(setdiff(all.vars(form), pnames), collapse=' + '))[[1]]
    nlform <- eval(nlform)
    environment(nlform) <- environment(form)
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
	       names(mc), 0)
    mc <- mc[c(1, m)]
    mc$drop.unused.levels <- TRUE
    mc[[1]] <- as.name("model.frame")
    mc$formula <- nlform
    fr <- eval(mc, parent.frame(2L))
    n <- nrow(fr)
    nlenv <- list2env(fr, parent=parent.frame(2L))
    lapply(pnames, function(nm) nlenv[[nm]] <- rep.int(nlpars[[nm]], n))
    respMod <- mkRespMod(fr, nlenv=nlenv, nlmod=nlmod)

    chck1(meform <- form[[3L]])
    pnameexpr <- parse(text=paste(pnames, collapse='+'))[[1]]
    nb <- nobars(meform)
    fe <- eval(substitute(~ 0 + nb + pnameexpr))
    environment(fe) <- environment(form)
    frE <- do.call(rbind, lapply(seq_along(nlpars), function(i) fr)) # rbind s copies of the frame
    for (nm in pnames) # convert these variables in fr to indicators
	frE[[nm]] <- as.numeric(rep(nm == pnames, each = n))
    X <- model.matrix(fe, frE)
    rownames(X) <- NULL

    reTrms <- mkReTrms(lapply(findbars(meform),
                              function(expr) {
                                  expr[[2]] <- substitute(0+foo, list(foo=expr[[2]]))
                                  expr
                              }), frE)
    list(respMod=respMod, frame=fr, X=X, reTrms=reTrms, pnames=pnames)
}

##' Create an object in a subclass of \code{\linkS4class{merMod}}
##' from the environment of the objective function and the value
##' returned by the optimizer.
##'
##' @title Create a merMod object
##' @param rho the environment of the objective function
##' @param opt the value returned by the optimizer
##' @param reTrms reTrms list from the calling function
##' @param fr model frame
##' @param mc matched call from the calling function
##' @return an object from a class that inherits from \code{\linkS4class{merMod}}
##' @export
mkMerMod <- function(rho, opt, reTrms, fr, mc) {
    if(missing(mc)) mc <- match.call()
    stopifnot(is.environment(rho),
              is(pp <- rho$pp, "merPredD"),
              is(resp <- rho$resp, "lmResp"),
              is.list(opt),
              "par" %in% names(opt),
              all(c("conv","fval") %in% substr(names(opt),1,4)), ## "conv[ergence]", "fval[ues]"
              is.list(reTrms),
              all(c("flist", "cnms", "Gp", "lower") %in%
                  names(reTrms)))
    rcl  <- class(resp)
    n    <- nrow(pp$V)
    p    <- ncol(pp$V)
    dims <- c(N=nrow(pp$X), n=n, p=p, nmp=n-p,
              nth=length(pp$theta), q=nrow(pp$Zt),
              nAGQ=rho$nAGQ,
              compDev=rho$compDev,
              ## 'use scale' in the sense of whether dispersion parameter should
              ##  be reported/used (*not* whether theta should be scaled by sigma)
              useSc=(rcl != "glmResp" ||
                     !resp$family$family %in% c("poisson","binomial")),
              reTrms=length(reTrms$cnms),
              spFe=0L,
              REML=if (rcl=="lmerResp") resp$REML else 0L,
              GLMM=(rcl=="glmResp"),
              NLMM=(rcl=="nlsResp"))
    storage.mode(dims) <- "integer"
    fac     <- as.numeric(rcl != "nlsResp")
    sqrLenU <- pp$sqrL(fac)
    wrss    <- resp$wrss()
    pwrss   <- wrss + sqrLenU
    weights <- resp$weights
    beta    <- pp$beta(fac)
    sigmaML <- pwrss/sum(weights)
    if (rcl != "lmerResp") {
        pars <- opt$par
        if (length(pars) > length(pp$theta)) beta <- pars[-(seq_along(pp$theta))]
    }
    cmp <- c(ldL2=pp$ldL2(), ldRX2=pp$ldRX2(), wrss=wrss,
             ussq=sqrLenU, pwrss=pwrss,
             drsum=if (rcl=="glmResp") resp$resDev() else NA,
             REML=if (rcl=="lmerResp" && resp$REML != 0L) opt$fval else NA,
             ## FIXME: construct 'REML deviance' here?
             dev=if (rcl=="lmerResp" && resp$REML != 0L) NA else opt$fval,
             sigmaML=sqrt(unname(if (!dims["useSc"]) NA else sigmaML)),
             sigmaREML=sqrt(unname(if (rcl!="lmerResp") NA else sigmaML*(dims['n']/dims['nmp']))),
             tolPwrss=rho$tolPwrss)
    # TODO:  improve this hack to get something in frame slot (maybe need weights, etc...)
    if(missing(fr)) fr <- data.frame(resp$y)
    new(switch(rcl, lmerResp="lmerMod", glmResp="glmerMod", nlsResp="nlmerMod"),
        call=mc, frame=fr, flist=reTrms$flist, cnms=reTrms$cnms,
        Gp=reTrms$Gp, theta=pp$theta, beta=beta, u=pp$u(fac),
        lower=reTrms$lower, devcomp=list(cmp=cmp, dims=dims), pp=pp, resp=resp,
        optinfo=list(optimizer=attr(opt,"optimizer"),
                     control=attr(opt,"control"),
                     conv=opt$conv,
                     warnings=attr(opt,"warnings"))
        )
}

## generic argument checking
## 'type': name of calling function ("glmer", "lmer", "nlmer")
##
checkArgs <- function(type,...) {
    l... <- list(...)
    if (isTRUE(l...[["sparseX"]])) warning("sparseX = TRUE has no effect at present")
    ## '...' handling up front, safe-guarding against typos ("familiy") :
    if(length(l... <- list(...))) {
        if (!is.null(l...[["family"]])) {  # call glmer if family specified
            ## we will only get here if 'family' is *not* in the arg list
            warning("calling lmer with family() is deprecated: please use glmer() instead")
            type <- "glmer"
        }
        ## Check for method argument which is no longer used
        ## (different meanings/hints depending on glmer vs lmer)
        if (!is.null(method <- l...[["method"]])) {
            msg <- paste("Argument", sQuote("method"), "is deprecated.")
            if (type=="lmer") msg <- paste(msg,"Use the REML argument to specify ML or REML estimation.")
            if (type=="glmer") msg <- paste(msg,"Use the nAGQ argument to specify Laplace (nAGQ=1) or adaptive",
                "Gauss-Hermite quadrature (nAGQ>1).  PQL is no longer available.")
            warning(msg)
            l... <- l...[names(l...) != "method"]
        }
        if(length(l...)) {
            warning("extra argument(s) ",
                    paste(sQuote(names(l...)), collapse=", "),
                    " disregarded")
        }
    }
}

## check formula and data: return an environment suitable for evaluating
##  the formula.
## (1) if data is specified, return it
## (2) otherwise, if formula has an environment, use it
## (3) otherwise [e.g. if formula was passed as a string], try to use parent.frame(2)

## if #3 is true *and* the user is doing something tricky with nested functions,
## this may fail ...

checkFormulaData <- function(formula,data,debug=FALSE) {
    dataName <- deparse(substitute(data))
    missingData <- inherits(tryCatch(eval(data), error=function(e)e), "error")
    ## data not found (this *should* only happen with garbage input,
    ## OR when strings used as formulae -> drop1/update/etc.)
    ##
    ## alternate attempt (fails)
    ##
    ## ff <- sys.frames()
    ## ex <- substitute(data)
    ## ii <- rev(seq_along(ff))
    ## for(i in ii) {
    ##     ex <- eval(substitute(substitute(x, env=sys.frames()[[n]]),
    ##                           env = list(x = ex, n=i)))
    ## }
    ## origName <- deparse(ex)
    ## missingData <- !exists(origName)
    ## (!dataName=="NULL" && !exists(dataName))
    if (missingData) {
        varex <- function(v,env) exists(v,envir=env,inherits=FALSE)
        allvars <- all.vars(as.formula(formula))
        allvarex <- function(vvec=allvars,...) { all(sapply(vvec,varex,...)) }
        if (allvarex(env=(ee <- environment(formula)))) {
            stop("'data' not found, but variables found in environment of formula: ",
                    "try specifying 'formula' as a formula rather ",
                    "than a string in the original model")
        } else stop("'data' not found, and some variables missing from formula environment")
    } else {
        if (is.null(data)) {
            if (!is.null(ee <- environment(formula))) {
                ## use environment of formula
                denv <- ee
            } else {
                ## e.g. no environment, e.g. because formula is a character vector
                ## parent.frame(2L) works because [g]lFormula (our calling environment)
                ## has been called within [g]lmer with env=parent.frame(1L)
                ## If you call checkFormulaData in some other bizarre way such that
                ## parent.frame(2L) is *not* OK, you deserve what you get
                ## calling checkFormulaData directly from the global
                ## environment should be OK, since trying to go up beyond the global
                ## environment keeps bringing you back to the global environment ...
                denv <- parent.frame(2L)
            }
        } else {
            ## data specified
            denv <- list2env(data)
        }
    }
    ## FIXME: set enclosing environment of denv to environment(formula), or parent.frame(2L) ?
    if (debug) {
        cat("Debugging parent frames in checkFormulaData:\n")
        ## find global environment -- could do this with sys.nframe() ?
        glEnv <- 1
        while (!identical(parent.frame(glEnv),.GlobalEnv)) {
            glEnv <- glEnv+1
        }
        ## where are vars?
        for (i in 1:glEnv) {
            OK <- allvarex(env=parent.frame(i))
            cat("vars exist in parent frame ",i)
            if (i==glEnv) cat(" (global)")
            cat(" ",OK,"\n")
        }
        cat("vars exist in env of formula ",allvarex(env=denv),"\n")
    } ## if (debug)

    stopifnot(length(as.formula(formula,env=denv)) == 3)  ## check for two-sided formula
    return(denv)
}

## checkFormulaData <- function(formula,data) {
##     ee <- environment(formula)
##     if (is.null(ee)) {
##         ee <- parent.frame(2)
##     }
##     if (missing(data)) data <- ee
##     stopifnot(length(as.formula(formula,env=as.environment(data))) == 3)
##     return(data)
## }



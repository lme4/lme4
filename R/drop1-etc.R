extractAIC.mer <- function(fit, scale = 0, k = 2, ...) {
    L <- logLik(fit)
    edf <- attr(L,"df")
    c(edf,-2*L + k*edf)
}

## doesn't install properly in current form
## should import nobs from stats?
## setMethod("nobs","mer",
##           function(object,...) {
##             nrow(object@frame)
##           })

## "Horribly" this is needed for stats::drop.scope() which does not see the S4 terms() method:
terms.mer <- function(x,...) attr(x@frame,"terms")


## hacked stats:::drop1.default
## FIXME: add F test (with specified denom df)?
drop1.mer <- function(object, scope, scale = 0, test = c("none", "Chisq"),
                      k = 2, trace = FALSE, ...)
{
    tl <- attr(terms(object), "term.labels")
    if(missing(scope)) scope <- drop.scope(object)
    else {
	if(!is.character(scope))
	    scope <- attr(terms(update.formula(object, scope)), "term.labels")
	if(!all(match(scope, tl, 0L) > 0L))
	    stop("scope is not a subset of term labels")
    }
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1L, ncol = 2L,
                  dimnames =  list(c("<none>", scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
    ## BMB: avoid nobs, to avoid dependence on 2.13
    ## n0 <- nobs(object, use.fallback = TRUE)
    n0 <- nrow(object@frame)
    env <- environment(formula(object))
    for(i in seq(ns)) {
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying -", tt, "\n", sep='')
	    utils::flush.console()
        }
        nfit <- update(object, as.formula(paste("~ . -", tt)),
                       evaluate = FALSE)
	nfit <- eval(nfit, envir = env) # was  eval.parent(nfit)
	ans[i+1, ] <- extractAIC(nfit, scale, k = k, ...)
        ## BMB: avoid nobs, to avoid dependence on 2.13
        ## nnew <- nobs(nfit, use.fallback = TRUE)
        nnew <- nrow(nfit@frame)
        if(all(is.finite(c(n0, nnew))) && nnew != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[1L , 1L] - ans[, 1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[,2])
    test <- match.arg(test)
    if(test == "Chisq") {
        dev <- ans[, 2L] - k*ans[, 1L]
        dev <- dev - dev[1L] ; dev[1L] <- NA
        nas <- !is.na(dev)
        P <- dev
        ## BMB: hack to extract safe_pchisq
        P[nas] <- stats:::safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
        aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    } else if (test == "F") {
        stop("F test STUB -- unfinished maybe forever")
        dev <- ans[, 2L] - k*ans[, 1L]
        dev <- dev - dev[1L] ; dev[1L] <- NA
        nas <- !is.na(dev)
        P <- dev
        ## BMB: hack to extract safe_pchisq
        P[nas] <- stats:::safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
        aod[, c("LRT", "Pr(F)")] <- list(dev, P)
    }
    head <- c("Single term deletions", "\nModel:", deparse(formula(object)),
	      if(scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

### Adapted from Rune Haubo's ordinal code
### extended convergence checking
### http://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions

##' @param derivs typically the "derivs" attribute of optimizeLmer(); with
##' "gradients" and possibly "Hessian" component
##' @param coefs estimated function value
##' @param ctrl list of lists, each with \code{action} character strings specifying
##' what should happen when a check triggers, and \code{tol} numerical tolerances,
##' as is the result of \code{\link{lmerControl}()$checkConv}.
##' @param lbound vector of lower bounds \emph{for random-effects parameters only}
##'   (length is taken to determine number of RE parameters)
##' @param debug useful if some checks are on "ignore", but would "trigger"
checkConv <- function(derivs, coefs, ctrl, lbound, debug = FALSE)
{

    if (is.null(derivs)) return(NULL)  ## bail out
    if (anyNA(derivs$gradient))
        return(list(code = -5L,
                    messages = gettextf("Gradient contains NAs")))
    ntheta <- length(lbound)
    res <- list()

    ## check singularity first, and unconditionally
    ## (ignore "ignore")
    ccl <- ctrl[[cstr <- "check.conv.singular"]] ; checkCtrlLevels(cstr, cc <- ccl[["action"]])
    ## similar logic to isSingular, but we don't have the fitted object to test
    bcoefs <- seq(ntheta)[lbound==0]
    is.singular <- any(coefs[bcoefs] < ccl$tol)
    
    if (doCheck(cc)) {
        ## singular fit
        ## are there other circumstances where we can get a singular fit?
        if (is.singular) {
            wstr <- "boundary (singular) fit: see ?isSingular"
            res$messages <- c(res$messages,wstr)
            switch(cc,
                   "message" = message(wstr),
                   "warning" = warning(wstr),
                   "stop" = stop(wstr),
                   stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
        }
    }

    ## DON'T check remaining gradient issues
    if (is.singular) return(res)
    
    ## gradients:
    ## check absolute gradient (default)
    ccl <- ctrl[[cstr <- "check.conv.grad"]] ; checkCtrlLevels(cstr, cc <- ccl[["action"]])
    wstr <- NULL
    if (doCheck(cc)) {
        scgrad <- tryCatch(with(derivs,solve(chol(Hessian),gradient)),
                           error=function(e)e)
        if (inherits(scgrad, "error")) {
            wstr <- "unable to evaluate scaled gradient"
            res$code <- -1L
        } else {
            ## find parallel *minimum* of scaled and absolute gradient
            ## the logic here is that we can sometimes get large
            ##  *scaled* gradients even when the *absolute* gradient
            ##  is small because the curvature is very flat as well ...
            mingrad <- pmin(abs(scgrad),abs(derivs$gradient))
            maxmingrad <- max(mingrad)
            if (maxmingrad > ccl$tol) {
                w <- which.max(maxmingrad)
                res$code <- -1L
                wstr <- gettextf("Model failed to converge with max|grad| = %g (tol = %g, component %d)",
                                 maxmingrad, ccl$tol,w)
            }
        }
        if (!is.null(wstr)) {
            res$messages <- wstr
            switch(cc,
                   "warning" = warning(wstr),
                   "stop" = stop(wstr),
                   stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
        }
        ## note: kktc package uses gmax > kkttol * (1 + abs(fval))
        ##  where kkttol defaults to 1e-3 and fval is the objective f'n value
        ## check relative gradient (only if enabled)
        if (!is.null(ccl$relTol) &&
            (max.rel.grad <- max(abs(derivs$gradient/coefs))) > ccl$relTol) {
                res$code <- -2L
                wstr <-
                    gettextf("Model failed to converge with max|relative grad| = %g (tol = %g)",
                             max.rel.grad, ccl$relTol)
                res$messages <- wstr
                switch(cc,
                       "warning" = warning(wstr),
                       "stop" = stop(wstr),
                       stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
        }
    }

    ccl <- ctrl[[cstr <- "check.conv.hess"]] ; checkCtrlLevels(cstr, cc <- ccl[["action"]])
    if (doCheck(cc)) {
        if (length(coefs) > ntheta) {
            ## GLMM, check for issues with beta parameters
            H.beta <- derivs$Hessian[-seq(ntheta),-seq(ntheta)]
            resHess <- checkHess(H.beta, ccl$tol, "fixed-effect")
            if (any(resHess$code!=0)) {
                res$code <- resHess$code
                res$messages <- c(res$messages,resHess$messages)
                wstr <- paste(resHess$messages,collapse=";")
                switch(cc,
                       "warning" = warning(wstr),
                       "stop" = stop(wstr),
                       stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
            }
        }
        resHess <- checkHess(derivs$Hessian, ccl$tol)
        if (any(resHess$code != 0)) {
            res$code <- resHess$code
            res$messages <- c(res$messages,resHess$messages)
            wstr <- paste(resHess$messages,collapse=";")
            switch(cc,
                   "warning" = warning(wstr),
                   "stop" = stop(wstr),
                   stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
        }
    }
    if (debug && length(res$messages) > 0) {
        print(res$messages)
    }
    res
}

checkHess <- function(H, tol, hesstype="") {
    ## FIXME: not sure why we decided to save messages as a list
    ## rather than as a character vector??
    res <- list(code=numeric(0),messages=list())
    evd <- tryCatch(eigen(H, symmetric=TRUE, only.values=TRUE)$values,
                    error=function(e)e)
    if (inherits(evd,"error")) {
        res$code <- -6L
        res$messages <- gettextf("Problem with Hessian check (infinite or missing values?)")
    } else {
        negative <- sum(evd < -tol)
        if(negative) {
            res$code <- -3L
            res$messages <-
                gettextf(paste("Model failed to converge:",
                               "degenerate",hesstype,"Hessian with %d negative eigenvalues"),
                         negative)
        } else {
            zero <- sum(abs(evd) < tol)
            if(zero || inherits(tryCatch(chol(H), error=function(e)e), "error")) {
                res$code <- -4L
                res$messages <-
                    paste(hesstype,"Hessian is numerically singular: parameters are not uniquely determined")
            } else {
                res$cond.H <- max(evd) / min(evd)
                if(max(evd) * tol > 1) {
                    res$code <- c(res$code, 2L)
                    res$messages <-
                        c(res$messages,
                          paste("Model is nearly unidentifiable: ",
                                "very large eigenvalue",
                                "\n - Rescale variables?", sep=""))
                }
                if ((min(evd) / max(evd)) < tol) {
                    res$code <- c(res$code, 3L)
                    ## consider skipping warning message if we've
                    ## already hit the previous flag?
                    if(!5L %in% res$code) {
                        res$messages <-
                            c(res$messages,
                              paste("Model is nearly unidentifiable: ",
                                    "large eigenvalue ratio",
                                    "\n - Rescale variables?", sep=""))
                    }
                }
            }
        }
    }
    if (length(res$code)==0) res$code <- 0
    res
}

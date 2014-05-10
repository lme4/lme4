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
    if (any(is.na(derivs$gradient)))
	return(list(code = -5L,
		    messages = gettextf("Gradient contains NAs")))
    ntheta <- length(lbound)
    res <- list()
    ## gradients:
    ## check absolute gradient (default)
    ccl <- ctrl[[cstr <- "check.conv.grad"]] ; checkCtrlLevels(cstr, cc <- ccl[["action"]])
    wstr <- NULL
    if (doCheck(cc)) {
        scgrad <- try(with(derivs,solve(Hessian,gradient)),silent=TRUE)
        if (inherits(scgrad,"try-error")) {
            wstr <- "unable to evaluate scaled gradient"
            res$code <- -1L
        } else {
            if ((max.grad <- max(abs(scgrad))) > ccl$tol) {
                w <- which.max(abs(scgrad))
                res$code <- -1L
                wstr <- gettextf("Model failed to converge with max|grad| = %g (tol = %g, component %d)",
                                 max.grad, ccl$tol,w)
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

    ccl <- ctrl[[cstr <- "check.conv.singular"]] ; checkCtrlLevels(cstr, cc <- ccl[["action"]])
    if (doCheck(cc)) {
        bcoefs <- seq(ntheta)[lbound==0]
        if (any(coefs[bcoefs] < ccl$tol)) {
            ## singular fit
            ## are there other circumstances where we can get a singular fit?
            wstr <- "singular fit"
            res$messages <- c(res$messages,wstr)
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
    res <- list(code=numeric(0),messages=list())
    evd <- eigen(H, symmetric=TRUE, only.values=TRUE)$values
    negative <- sum(evd < -tol)
    if(negative) {
        res$code <- -3L
        res$messages <-
            gettextf(paste("Model failed to converge:",
                           "degenerate",hesstype,"Hessian with %d negative eigenvalues"),
                     negative)
    } else {
        zero <- sum(abs(evd) < tol)
        ch <- try(chol(H), silent=TRUE)
        if(zero || inherits(ch, "try-error")) {
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
    if (length(res$code)==0) res$code <- 0
    res
}


## original ordinal-package code
##convCheck_orig##  <-
##     function(fit, control=NULL, Theta.ok=NULL, tol=sqrt(.Machine$double.eps), ...)
## ### check convergence along the
## ### way.
## ### fit: clm-object or the result of clm.fit.NR()
## ###
## ### Return: list with elements
## ### vcov, conv, cond.H, messages and
## {
##     if(missing(control))
##         control <- fit$control
##     if(is.null(control))
##         stop("'control' not supplied - cannot check convergence")

##     if(!is.null(control$tol))
##         tol <- control$tol
##     if(tol < 0)
##         stop(gettextf("numerical tolerance is %g, expecting non-negative value",
##                       tol), call.=FALSE)
## ### FIXME: test this.
##     H <- fit$Hessian
##     g <- fit$gradient
##     max.grad <- max(abs(g))
##     cov <- array(NA_real_, dim=dim(H), dimnames=dimnames(H))
##     cond.H <- NA_real_
##     res <- list(vcov=vcov, code=integer(0L), cond.H=cond.H,
##                 messages=character(0L))
##     class(res) <- "conv.check"
##     if(is.list(code <- fit$convergence))
##         code <- code[[1L]]
##     mess <-
##         switch(as.character(code),
##                "0" = "Absolute and relative convergence criteria were met",
##                "1" = "Absolute convergence criterion was met, but relative criterion was not met",
##                "2" = "iteration limit reached",
##                "3" = "step factor reduced below minimum",
##                "4" = "maximum number of consecutive Newton modifications reached")
##     res <- c(res, alg.message=mess)
##     ## }
##     if(max.grad > control$gradTol) {
##         res$code <- -1L
##         res$messages <-
##             gettextf("Model failed to converge with max|grad| = %g (tol = %g)",
##                      max.grad, control$gradTol)
##         return(res)
##     }
##     evd <- eigen(H, symmetric=TRUE, only.values=TRUE)$values
##     negative <- sum(evd < -tol)
##     if(negative) {
##         res$code <- -2L
##         res$messages <-
##             gettextf(paste("Model failed to converge:",
##                            "degenerate Hessian with %d negative eigenvalues"),
##                      negative)
##         return(res)
##     }
##     if(!is.null(Theta.ok) && !Theta.ok) {
##         res$code <- -3L
##         res$messages <-
##             "not all thresholds are increasing: fit is invalid"
##         return(res)
##     }
##     zero <- sum(abs(evd) < tol)
##     ch <- try(chol(H), silent=TRUE)
##     if(zero || inherits(ch, "try-error")) {
##         res$code <- 1L
##         res$messages <-
##             "Hessian is numerically singular: parameters are not uniquely determined"
##         return(res)
##     }
##     ## Add condition number to res:
##     res$cond.H <- max(evd) / min(evd)
##     ## Compute Newton step:
##     step <- c(backsolve(ch, backsolve(ch, g, transpose=TRUE)))
##     ## Compute var-cov:
##     res$vcov[] <- chol2inv(ch)
##     if(max(abs(step)) > control$relTol) {
##         res$code <- c(res$code, 1L)
##         corDec <- as.integer(min(cor.dec(step)))
##         res$messages <-
##             c(res$messages,
##               gettextf("some parameters may have only %d correct decimals",
##                        corDec))
##     }
##     if(max(evd) * tol > 1) {
##         res$code <- c(res$code, 2L)
##         res$messages <-
##             c(res$messages,
##               paste("Model is nearly unidentifiable: ",
##                     "very large eigenvalue",
##                     "\n - Rescale variables?", sep=""))
##     }
##     if((min(evd) / max(evd)) < tol) {
##         res$code <- c(res$code, 3L)
##         if(!5L %in% res$code) {
##             res$messages <-
##                 c(res$messages,
##                   paste("Model is nearly unidentifiable: ",
##                         "large eigenvalue ratio",
##                         "\n - Rescale variables?", sep=""))
##         }
##     }
##     if(!length(res$code)) {
##         res$code <- 0L
##         res$messages <- "successful convergence"
##     }
##     res
## }

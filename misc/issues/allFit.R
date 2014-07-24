##
library(lme4)
require(optimx)
require(nloptr)

namedList <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)),deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
    setNames(L,nm)
}

## originally from https://github.com/lme4/lme4/issues/98 :

nloptWrap <- function(fn, par, lower, upper, control=list(), ...) {
    defaultControl <- list(xtol_rel = 1e-6, maxeval = 1e5)
    for (n in names(defaultControl))
      if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
    res <- nloptr(x0=par, eval_f=fn, lb=lower, ub=upper, opts=control, ...)
    ##     ------
    with(res,list(par=solution,
                  fval=objective,
                  feval=iterations,
                  conv=if (status>0) 0 else status,
                  message=message))
}

##' Attempt to re-fit a [g]lmer model with a range of optimizers.
##' The default is to use all known optimizers for R that satisfy the
##' requirements (do not require explicit gradients, allow
##' box constraints), in three categories; (i) built-in
##' (minqa::bobyqa, lme4::Nelder_Mead), (ii) wrapped via optimx
##' (most of optimx's optimizers that allow box constraints require
##' an explicit gradient function to be specified; the two provided
##' here are really base R functions that can be accessed via optimx,
##' (iii) wrapped via nloptr.
##'
##' @param m a fitted model
##' @param meth.tab a matrix (or data.frame) with columns
##' - method  the name of a specific optimization method to pass to the optimizer
##'           (leave blank for built-in optimizers)
##' - optimizer  the \code{optimizer} function to use
##' @param verbose print progress messages?
##' @return a list of fitted \code{merMod} objects
##' @seealso slice, slice2D in the bbmle package
##' @examples
##' library(lme4)
##' gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'                  data = cbpp, family = binomial)
##' gm_all <- allFit(gm1)
##' t(sapply(gm_all,fixef))              ## extract fixed effects
##' sapply(gm_all,logLik)                ## log-likelihoods
##' sapply(gm_all,getME,"theta")         ## theta parameters
##' !sapply(gm_all,inherits,"try-error") ## was fit OK?
allFit <- function(m, meth.tab = cbind(optimizer=
                      rep(c("bobyqa","Nelder_Mead", "optimx",  "nloptWrap"),
                          c(    1,         1,           2,         2)),
                      method= c("",        "",  "nlminb","L-BFGS-B",
                               "NLOPT_LN_NELDERMEAD", "NLOPT_LN_BOBYQA")),
                   verbose=TRUE,
                   
                   maxfun=1e5)
{
    stopifnot(length(dm <- dim(meth.tab)) == 2, dm[1] >= 1, dm[2] >= 2,
	      is.character(optimizer <- meth.tab[,"optimizer"]),
	      is.character(method    <- meth.tab[,"method"]))
    fit.names <- paste(optimizer, method, sep=".")
    res <- setNames(as.list(fit.names), fit.names)
    for (i in seq_along(fit.names)) {
        if (verbose) cat(fit.names[i],": ")
        ctrl <- list(optimizer=optimizer[i])
        ctrl$optCtrl <- switch(optimizer[i],
                               optimx    = list(method   = method[i]),
                               nloptWrap = list(algorithm= method[i]),
                               list(maxfun=maxfun))
        ctrl <- do.call(if(isGLMM(m)) glmerControl else lmerControl, ctrl)
        tt <- system.time(rr <- tryCatch(update(m, control = ctrl), error = function(e) e))
        attr(rr, "optCtrl") <- ctrl$optCtrl # contains crucial info here
        attr(rr, "time") <- tt  # store timing info
        res[[i]] <- rr
        if (verbose) cat("[OK]\n")
    }
    ## 
    res
}

summary.allfit <- function(object, ...) {
    which.OK <- !sapply(object,is,"error")
    msgs <- lapply(object[which.OK],function(x) x@optinfo$conv$lme4$messages)
    fixef <- t(sapply(object[which.OK],fixef))
    llik <- sapply(object[which.OK],logLik)
    times <- t(sapply(object[which.OK],attr,"time"))
    feval <- sapply(object[which.OK],function(x) x@optinfo$feval)
    sdcor <- t(sapply(object[which.OK],function(x) {
        aa <- as.data.frame(VarCorr(x))
        setNames(aa[,"sdcor"],lme4:::tnames(object[which.OK][[1]]))
    }))
    namedList(which.OK,msgs,fixef,llik,sdcor,times,feval)
}
print.summary.allfit <- function(object,...) {
    if (!which.OK==seq(length(object))) {
        cat("some optimizers failed: ",
            paste(names(object)[!which.OK],collapse=","),"\n")
    }

}

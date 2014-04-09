##
library(lme4)
require(optimx)
require(nloptr)

## copied from https://github.com/lme4/lme4/issues/98:
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-6,maxeval=1e5)
nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
    for (n in names(defaultControl)) 
      if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
    res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
    with(res,list(par=solution,
                  fval=objective,
                  feval=iterations,
                  conv=if (status>0) 0 else status,
                  message=message))
}

mtab <- data.frame(optimizer=rep(c("bobyqa","Nelder_Mead",
                   "optimx","nloptwrap2"),c(1,1,2,2)),
                   method=c("","","nlminb","L-BFGS-B",
                   "NLOPT_LN_NELDERMEAD","NLOPT_LN_BOBYQA"),
                   stringsAsFactors=FALSE)

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
##' @param method the name of a specific optimization method to pass to the optimizer (leave blank for built-in optimizers)
##' @param optimizer the \code{optimizer} function to use
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
allFit <- function(m,   
                   method=mtab$method,
                   optimizer=mtab$optimizer,
                   verbose=TRUE) {
    res <- list()
    fit.names <- paste(optimizer,method,sep=".")
    for (i in seq_along(fit.names)) {
        if (verbose) cat(fit.names[i],"\n")
        c.arglist <- list()
        c.arglist <- list(optimizer=optimizer[i])
        c.arglist$optCtrl <- switch(optimizer[i],
                                    optimx=list(method=method[i]),
                                    nloptwrap2=list(algorithm=method[i]),
                                    NULL)
        cfun <- if (isGLMM(m)) glmerControl else lmerControl
        res[[i]] <- try(update(m,control=do.call(cfun,c.arglist)),silent=TRUE)
    }
    names(res) <- fit.names
    return(res)
}


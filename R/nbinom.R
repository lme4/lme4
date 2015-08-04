##' @importFrom MASS negative.binomial
##' @importFrom MASS theta.ml

## ==> user should use  getME(object, "glmer.nb.theta")
getNBdisp <- function(object) environment(object@resp$family$aic)[[".Theta"]]

## Hidden, originally used at least once, well tested (!), but
if(FALSE) # not needed anymore currently
getNBdisp.fam <- function(familyString)
    as.numeric(sub(".*([-+]?\\d+(\\.\\d*)?([Ee][-+]?\\d+)?){1}.*",
                   "\\1", familyString))

## Package "constants" {only on depending the glmResp definition in ./AllClass.R}:
glmResp.f.nms <- names(glmResp$fields())
glmNB.to.change <- setdiff(glmResp.f.nms, c("Ptr", "family"))

##' setNBdisp(object,theta) :=
##' NB-object with changed [DISP]ersion parameter 'theta' (and all that entails)
setNBdisp <- function(object,theta) {
  ## assign(".Theta",theta,envir=environment(object@resp$family$aic))
  rr <- object@resp
  newresp <- do.call(glmResp$new,
                     c(lapply(setNames(nm=glmNB.to.change), rr$field),
                       list(family = negative.binomial(theta=theta))))
  newresp$setOffset(rr$offset)
  newresp$updateMu(rr$eta - rr$offset)
  object@resp <- newresp
  object
}

refitNB <- function(object, theta, control = NULL) {
  refit(setNBdisp(object, theta), control = control)
}

##' @title Optimize over the Negative Binomial Parameter Theta
##' @param object a "glmerMod" object, updated from poisson to negative.binomial()
##' @param interval and
##' @param tol      are both passed to \code{\link{optimize}()}.
##' @param verbose  ## show our own progress
##' @param control passed to \code{\link{refit}()}
##' @return the last fit, an object like 'object'
optTheta <- function(object,
                     interval = c(-5,5),
                     tol = .Machine$double.eps^0.25,
                     verbose = FALSE,
                     control = NULL)
{
  lastfit <- object
  it <- 0L
  NBfun <- function(t) {
      ## Kluge to retain last value and evaluation count {good enough for ..}
      dev <- -2*logLik(lastfit <<- refitNB(lastfit,
                                           theta = exp(t),
                                           control = control))
      it <<- it+1L
      if (verbose)
          cat(sprintf("%2d: th=%#15.10g, dev=%#14.8f, beta[1]=%#14.8f\n",
                      it, exp(t), dev, lastfit@beta[1]))
      dev
  }
  optval <- optimize(NBfun, interval=interval, tol=tol)
  stopifnot(all.equal(optval$minimum, log(getNBdisp(lastfit))))
  ## FIXME: return eval count info somewhere else? MM: new slot there, why not?
  attr(lastfit,"nevals") <- it
  ## fix up the 'th' expression, replacing it by the real number,
  ## so effects:::mer.to.glm() can eval() it:
  lastfit@call$family[["theta"]] <- exp(optval$minimum)
  lastfit
}

## use MASS machinery to estimate theta from residuals
est_theta <- function(object, limit = 20,
                      eps = .Machine$double.eps^0.25, trace = 0)
{
  Y <- model.response(model.frame(object))
  mu <- fitted(object)
  theta.ml(Y, mu, weights = object@resp$weights,
           limit = limit, eps = eps, trace = trace)
}

## -------> ../man/glmer.Rd
##' glmer() for Negative Binomial
##' @param ... formula, data, etc: the arguments for
##' \code{\link{glmer}(..)} (apart from \code{family}!).
glmer.nb <- function(..., interval = log(th) + c(-3,3),
                     tol = 5e-5, verbose = FALSE, nb.control = NULL,
                     initCtrl = list(limit = 20, eps = 2*tol, trace = verbose))
{
    dotE <- as.list(substitute(E(...))[-1])
    ## nE <- names(dotE <- as.list(substitute(E(...))[-1]))
    ## i <- match(c("formula",""), nE)
    ## i.frml <- i[!is.na(i)][[1]] # the index of "formula" in '...'
    ## dots <- list(...)
    g0 <- glmer(..., family = poisson, verbose = verbose >= 2)
    th <- est_theta(g0, limit = initCtrl$limit,
		    eps = initCtrl$eps, trace = initCtrl$trace)
    ## using format() on purpose, influenced by options(digits = *) :
    if(verbose) cat("th := est_theta(glmer(..)) =", format(th))
    g1 <- update(g0, family = negative.binomial(theta=th))
    if(verbose) cat(" --> dev.= -2*logLik(.) =", format(-2*logLik(g1)),"\n")
    ## fix the 'data' part (only now!)
    if("data" %in% names(g1@call))
	g1@call[["data"]] <- dotE[["data"]]
    else
        warning("no 'data = *' in glmer.nb() call.. Not much is guaranteed")
    if ("verbose" %in% names(g1@call)) {
        g1@call[["verbose"]] <- dotE[["verbose"]]
    }
    optTheta(g1, interval=interval, tol=tol, verbose=verbose,
	     control = nb.control)
}

## do we want to facilitate profiling on theta??
## save evaluations used in optimize() fit?
## ('memoise'?)
## Again, I think that a reference class object would be a better approach.

##' @importFrom MASS negative.binomial
##' @importFrom MASS theta.ml

## should be getME(object,"NBdisp") ?
## MM: should the *user* use it?  if yes, consider method sigma() or AIC() ?
getNBdisp <- function(object) {
  get(".Theta",envir=environment(object@resp$family$aic))
}


glmResp.f.nms <- names(glmResp$fields())

## should be setME(object,"NBdisp") ?
setNBdisp <- function(object,theta) {
  ## assign(".Theta",theta,envir=environment(object@resp$family$aic))
  ff <- setdiff(names(getRefClass("glmResp")$fields()),c("Ptr","family"))
  rr <- object@resp
  arg1 <- lapply(ff,rr$field)
  names(arg1) <- ff
  newresp <- do.call(glmResp$new,
                     c(arg1, list(family=negative.binomial(theta=theta))))
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
  optval <- optimize(function(t) {
      ## Kluge to retain last value and evaluation count {good enough for ..}
      dev <- -2*logLik(lastfit <<- refitNB(lastfit,
                                           theta = exp(t),
                                           control = control))
      it <<- it+1L
      if (verbose) cat(sprintf("%2d: th=%#15.10g, dev=%#14.8f\n", it, exp(t), dev))
      dev
  }, interval=interval, tol=tol)
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
    if(verbose) cat("th := est_theta(glmer(..)) =", format(th),"\n")
    g1 <- update(g0, family = negative.binomial(theta=th))
    ## fix the 'data' part (only now!)
    if("data" %in% names(g1@call))
	g1@call[["data"]] <- dotE[["data"]]
    else
        warning("no 'data = *' in glmer.nb() call.. Not much is guaranteed")
    optTheta(g1, interval=interval, tol=tol, verbose=verbose,
	     control = nb.control)
}

## do we want to facilitate profiling on theta??
## save evaluations used in optimize() fit?
## ('memoise'?)
## Again, I think that a reference class object would be a better approach.

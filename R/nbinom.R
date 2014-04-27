##' @importFrom MASS negative.binomial
##' @importFrom MASS theta.ml

## should be getME(object,"NBdisp") ?
## MM: should the *user* use it?  if yes, consider method sigma() or AIC() ?
getNBdisp <- function(object) {
  get(".Theta",envir=environment(object@resp$family$aic))
}


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

optTheta <- function(object,
                     interval = c(-5,5),
                     maxit = 20,
                     verbose = FALSE,
                     control = NULL) {
  lastfit <- object
  evalcnt <- 0
  optval <- optimize(function(t) {
      ## FIXME: kluge to retain last value and evaluation count
      ## Perhaps use a reference class object to keep track of this
      ## auxilliary information?  DB
      dev <- deviance(lastfit <<- refitNB(lastfit,
                                          theta = exp(t),
                                          control = control))
      evalcnt <<- evalcnt+1
      if (verbose) cat(evalcnt,exp(t),dev,"\n")
      dev
  }, interval=interval)
  stopifnot(all.equal(optval$minimum,log(getNBdisp(lastfit))))
  ## FIXME: return eval count info somewhere else? MM: new slot there, why not?
  attr(lastfit,"nevals") <- evalcnt
  lastfit
}

## use MASS machinery to estimate theta from residuals
est_theta <- function(object) {
  Y <- model.response(model.frame(object))
  mu <- fitted(object)
  w <- object@resp$weights
  control <- list(maxit=20,trace=0)
  th <- theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                 trace = control$trace > 2)
}

## FIXME: really should document glmer.nb() on the same help page as glmer()
## I (MM) don't know how to use roxygen for that..

##' glmer() for Negative Binomial
##' @param ... formula, data, etc: the arguments for
##' \code{\link{glmer}(..)} (apart from \code{family}!).
##' @param interval interval in which to start the optimization
##' @param verbose logical indicating how much progress information
##' should be printed.
##' @export
glmer.nb <- function(..., interval = log(th)+c(-3,3), verbose=FALSE)
{
    control <- list(...)$control
    g0 <- glmer(..., family = poisson)
    th <- est_theta(g0)
    if(verbose) cat("th := est_theta(glmer(..)) =", format(th),"\n")
    g1 <- update(g0, family = negative.binomial(theta=th))
    ## if (is.null(interval)) interval <- log(th)+c(-3,3)
    optTheta(g1, interval = interval, verbose = verbose, control = control)
}

## do we want to facilitate profiling on theta??
## save evaluations used in optimize() fit?
## ('memoise'?)
## Again, I think that a reference class object would be a better approach.

##' @importFrom MASS negative.binomial
##' @importFrom MASS theta.ml
##require(MASS)

## should be getME(object,"NBdisp") ?
getNBdisp <- function(object) { 
  get(".Theta",envir=environment(object@resp$family$aic))
}


## should be setME(object,"NBdisp") ?
setNBdisp <- function(object,theta) {
  ## assign(".Theta",theta,envir=environment(object@resp$family$aic))
  ff <- setdiff(names(getRefClass("glmResp")$fields()),c("Ptr","family"))
  arg1 <- lapply(ff,object@resp$field)
  names(arg1) <- ff
  newresp <- do.call(glmResp$new,c(arg1,
                                   list(family=negative.binomial(theta=theta))))
  object@resp <- newresp
  object
}

refitNB <- function(object,theta) {
  orig_theta <- getNBdisp(object)
  object <- setNBdisp(object,theta)  ## new/copied object
  refit(object,newresp=model.response(model.frame(object)))
  ## FIXME: should refit() take this response as a default??
  ## Yes, I think that is a good idea.  DB 2012-02-23
}

optTheta <- function(object,
                     interval=c(-5,5),
                     maxit=20,
                     debug=FALSE) {
  lastfit <- object
  evalcnt <- 0
  optval <- optimize(function(t) {
    ## FIXME: kluge to retain last value and evaluation count
      ## Perhaps use a reference class object to keep track of this
      ## auxilliary information?  DB
    L <- -logLik(lastfit <<- refitNB(lastfit,theta=exp(t)))
    evalcnt <<- evalcnt+1
    if (debug) {
         cat(evalcnt,exp(t),L,"\n")
      }
      L
  },interval=interval)
  stopifnot(all.equal(optval$minimum,log(getNBdisp(lastfit))))
  ## FIXME: return eval count info somewhere else?
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

## wrapper for glmer stuff
glmer.nb <- function(...,
                     interval=NULL,
                     debug=FALSE) {
  g0 <- glmer(...,family=poisson)
  th <- est_theta(g0)
  g1 <- update(g0,family=negative.binomial(theta=th))
  if (is.null(interval)) interval <- log(th)+c(-3,3)
  optTheta(g1,interval=interval,debug=debug)
}

## do we want to facilitate profiling on theta??
## save evaluations used in optimize() fit?
## ('memoise'?)
## Again, I think that a reference class object would be a better approach.

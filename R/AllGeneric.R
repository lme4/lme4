## utilities, these *exported*:
##' @export getL
setGeneric("getL", function(x) standardGeneric("getL"))

fixed.effects <- function(object, ...) {
    ## fixed.effects was an alternative name for fixef
    .Deprecated("fixef")
    mCall = match.call()
    mCall[[1]] = as.name("fixef")
    eval(mCall, parent.frame())
}

random.effects <- function(object, ...) {
    ## random.effects was an alternative name for ranef
    .Deprecated("ranef")
    mCall = match.call()
    mCall[[1]] = as.name("ranef")
    eval(mCall, parent.frame())
}

## Create a Markov chain Monte Carlo sample from the posterior
## distribution of the parameters
##
##
## @title Create an MCMC sample
## @param object a fitted model object
## @param n number of samples to generate.  Defaults to 1; for real use values of 200-1000 are more typical
## @param verbose should verbose output be given?
## @param ... additional, optional arguments (not used)
## @return a Markov chain Monte Carlo sample as a matrix
mcmcsamp <- function(object, n = 1L, verbose = FALSE, ...) UseMethod("mcmcsamp")

if(getRversion() < "3.3") {
    sigma <- function(object, ...) UseMethod("sigma")
}

isREML <- function(x, ...) UseMethod("isREML")

isLMM <- function(x, ...) UseMethod("isLMM")
isGLMM <- function(x, ...) UseMethod("isGLMM")
isNLMM <- function(x, ...) UseMethod("isNLMM")


##' Refit a model using the maximum likelihood criterion
##'
##' This function is primarily used to get a maximum likelihood fit of
##' a linear mixed-effects model for an \code{\link{anova}} comparison.
##'
##' @title Refit a model by maximum likelihood criterion
##' @param x a fitted model, usually of class \code{"\linkS4class{lmerMod}"},
##'     to be refit according to the maximum likelihood criterion
##' @param ... optional additional parameters.  None are used at present.
##' @return an object like \code{x} but fit by maximum likelihood
##' @export
refitML <- function(x, ...) UseMethod("refitML")

refit <- function(object, newresp, ...) UseMethod("refit")

if (FALSE) {
setGeneric("HPDinterval",
           function(object, prob = 0.95, ...) standardGeneric("HPDinterval"))
}

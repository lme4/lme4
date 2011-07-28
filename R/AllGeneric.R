setGeneric("lmList",
           function(formula, data, family, subset, weights,
                    na.action, offset, pool, ...)
           standardGeneric("lmList"))

## if (getRversion() < '2.14.0') {
##     getCall <- function(x, ...) UseMethod("getCall")
## }

## utilities, these *exported*:
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

setGeneric("sigma", function(object, ...) standardGeneric("sigma"))

if (FALSE) {
setGeneric("HPDinterval",
           function(object, prob = 0.95, ...) standardGeneric("HPDinterval"))
}

if (FALSE) {
setGeneric("mcmcsamp",
           function(object, n = 1, verbose = FALSE, ...)
           standardGeneric("mcmcsamp"))
}

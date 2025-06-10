# 1. Virtual Class: Virtual Covariance

##' Virtual Class for Covariance Structures
##'
##' This abstract S4 class defines the common interface and required methods
##' for all specific covariance structure implementations within the package.
##' It serves as a base class for inheritance, ensuring consistent behavior
##' across different covariance types.
##'
##' @slot dimension An integer indicating the size of the square covariance matrix (N x N).
##' @slot parameters A numeric vector storing the unconstrained parameters used by the optimizer.
##' @export
setClass("VirtualCovariance", contains = "VIRTUAL",
    slots = c(
        dimension = "integer",
        parameters = "numeric"
    ),
    prototype = list(
        dimension = 0L,
        parameters = numeric(0)
    )
)

# Define ALL Generic Functions for Covariance Structures

##' Get Number of Parameters
##'
##' Retrieves the number of unconstrained parameters for a covariance object. 
##' @param object An S4 object inheriting from `VirtualCovariance`
##' @return An integer representing the number of unconstrained parameter
##' export
setGeneric("n_parameters", function(object) standardGeneric("n_parameters"))

##' Get Unconstrained Parameters
##' 
##' Returns the current unconstrained parameter vector used for optimization.
##' @param object An S4 object inheriting from `VirtualCovariance`.
##' @return A numeric vector of unconstrained parameters.
##' @export
setGeneric("get_parameters", function(object), standardGeneric("get_parameters"))

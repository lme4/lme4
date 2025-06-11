
# S4 Class and Method Definition for Covariance Structures

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

# Generic Functions for Covariance Structures 

##' Get Number of Parameters
##'
##' Generic function to retrieve the number of unconstrained parameters for a covariance object. 
##' @param object An S4 object inheriting from `VirtualCovariance`,
##' @return An integer representing the number of unconstrained parameters.
##' @export
setGeneric("n_parameters", function(object) standardGeneric("n_parameters"))

##' Get Unconstrained Parameters
##' 
##' Generic function to return the current unconstrained parameter vector used for optimization.
##' @param object An S4 object inheriting from `VirtualCovariance`.
##' @return A numeric vector of unconstrained parameters.
##' @export
setGeneric("get_parameters", function(object) standardGeneric("get_parameters"))

##' Set Unconstrained Parameters
##'
##' Generic function to update the internal parameters of a covariance object from a provided unconstrained vector. 
##' @param object An S4 object inheriting from `VirtualCovariance`.
##' @param value A numeric vector of unconstrained parameters. 
##' @return The modified covariance object. 
##' @export
setGeneric("set_parameters", function(object,value) standardGeneric("set_parameters"))

##' Compute Covariance Matrix
##' 
##' Generic function to compute the full covariance function matrix (Sigma) based on current parameters and auxillary data. 
##' @param object  An S4 object inheriting from `Virtual Covariance`.
##' @param data_context Optional. Auxillary data needed for computation (e.g., grouping factors, time points, coordinates).
##' @return A numeric matrix representing the covariance matrix. 
##' @export
setGeneric("compute_inverse_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_inverse_covariance_matrix"))

##' Compute Log-Determinant of Covariance Matrix
##'  
##' Generic function to compute the log-determinant of the covariance matrix.
##' @param object An S4 object inheriting from `VirtualCovariance`. 
##' @param data_context Optional. Auxillary data needed for computation. 
##' @return A numeric value representing the log-determinant. 
##' @export
setGeneric("compute_log_det_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_log_det_covariance_matrix"))

##' Compute Inverse Covariance Matrix
##'
##' Generic function to compute the inverse of the covariance matrix (Sigma_inverse).
##' @param object An S4 object inheriting from `VirtualCovariance`.
##' @param data_context Optional. Auxillary data needed for computation. 
##' @return A numeric matrix representing the Cholesky factor. 
##' @export 
setGeneric("get_cholesky_factor", function(object, data_context = NULL) standardGeneric("compute_inverse_covariance_matrix"))

##' Get Cholesky Factor
##'
##' Generic function to compute the Cholesky factor (e.g., L) such that Sigma = L %*% t(L).
##' @param object An S4 object inheriting from `VirtualCovariance`.
##' @param data_context Optional. Auxillary data needed for computation. 
##' @return A numeric matrix representing the Cholesky factor. 
##' @export
setGeneric("get_cholesky_factor", function(object, data_context = NULL) standardGeneric("get_cholesky_factor"))

##' Get Start Values for Parameters
##' 
##' Generic function to provide sensible initial values for the unconstrained parameters
##' for optimization.
##' @param object An S4 object inheriting from `VirtualCovariance`. 
##' @return A numeric vector of initial unconstrained parameter values. 
##' @export
setGeneric("get_start_values", function(object) standardGeneric("get_start_values"))

##' Validate Parameters
##' 
##' Generic function to validate if the current parameters lead to a valid
##' covariance matrix (e.g., positive definite). 
##' @param object An S4 object inheriting from `VirtualCovariance`. 
##' @return 'TRUE' if parameters are valid, otherwise throws an error. 
##' @export
setGeneric("validate_parameters", function(object) standardGeneric("validate_parameters"))

##' Get Interpretable Parameters
##' 
##' Generic function to return interpretable parameters (e.g., variance, correlation,
##' scale) from the unconstrained optimization parameters. 
##' @param object An S4 object inheriting from `VirtualCovariance`.
##' @return A named list of interpretable parameters. 
##' @export
setGeneric("get_interpretable_parameters", function(object) standardGeneric("get_interpretable_parameters"))

##' Check if Covariance Matrix is Diagonal 
##' 
##' Generic function to check if the covariance matrix represented by the object is diagonal. 
##' @param object. An S4 object inheriting from `VirtualCovariance`.
##' @return 'TRUE' if the matrix is diagonal, 'FALSE' otherwise. 
##' @export
setGeneric("is_diagonal", function(object) standardGeneric("is_diagonal"))

# Virtual Classes for Parameterization Types
# These serve as markers for method dispatch based on parameterization strategy.

##' Virtual Class for Log-Cholesky Parameterization
##'
##' This virtual class indicates that the covariance structure's parameters
##' are handled using a log-Cholesky decomposition for unconstrained optimization.
##'
##' @export
setClass("VirtualParameterizationLogChol", contains = "VIRTUAL")

##' Virtual Class for Log-Scale and Bounded Correlation Parameterization
##'
##' This virtual class indicates that the covariance structure's parameters
##' involve log-transformed scales/variances and transformed (e.g., atanh)
##' correlations or shape parameters for unconstrained optimization.
##'
##' @export
setClass("VirtualParameterizationLogScaleBoundedCor", contains = "VIRTUAL")

# Concrete Covariance Structure Classes 

##' Diagonal Covariance Structure 
##'
##' Represents a diagonal covariance matrix, where random effects are
##' independent and potentially have heterogeneous variances.
##'
##' @slot internal_diag_values A numeric vector storing the actual positive
##' diagonal variance values. 
##' @export
setClass("DiagonalCov", 
	 contains = c("VirtualCovariance", "VirtualParameterizationLogScaleBoundedCor"),
	 slots = c(
	     internal_diag_values = "numeric"
	 ),
	 prototype = list(
	     dimension = 0L,
	     parameters = numeric(0),
	     internal_diag_values = numeric(0)
	)
)

# Method for DiagonalCov 

##' @rdname n_parameters
setMethod("n_parameters", "DiagonalCov", function(object) {
	object@dimension
})

##' @rdname get_parameters
setMethod("get_parameters", "DiagonalCov", function(object) {
	log(object@internal_diag_values)
})

##' @rdname set_parameters
setMethod("set_parameters", "DiagonalCov", function(object, value) {
	if (length(value) != object@dimension) {
		stop("Incorrect number of parameters for DiagonalCov. Expected ", object@dimension)
	}
	object@internal_diag_values <- exp(value)
	object@parameters <- value 
	object
})

##' @rdname compute_covariance_matrix
setMethod("compute_covariance_matrix", "DiagonalCov", function (object, data_context = NULL) {
	  diag(object@internal_diag_values)
})
					   

##' @rdname compute_log_det_covariance_matrix
setMethod("compute_log_det_covariance_matrix", "DiagonalCov", function(object, data_context = NULL) {
    sum(log(object@internal_diag_values))
})

##' @rdname compute_inverse_covariance_matrix
setMethod("compute_inverse_covariance_matrix", "DiagonalCov", function(object, data_context = NULL) {
    diag(1 / object@internal_diag_values)
})

##' @rdname get_cholesky_factor
setMethod("get_cholesky_factor", "DiagonalCov", function(object, data_context = NULL) {
    diag(sqrt(object@internal_diag_values))
})

##' @rdname get_start_values
setMethod("get_start_values", "DiagonalCov", function(object) {
    rep(log(1), object@dimension) # Start with log(1) for all variances (i.e., variance = 1)
})

##' @rdname validate_parameters
setMethod("validate_parameters", "DiagonalCov", function(object) {
    if (any(object@internal_diag_values <= 0)) {
        stop("DiagonalCov: All variances must be positive.")
    }
    TRUE
})

##' @rdname get_interpretable_parameters
setMethod("get_interpretable_parameters", "DiagonalCov", function(object) {
    list(variances = object@internal_diag_values)
})

##' @rdname is_diagonal
setMethod("is_diagonal", "DiagonalCov", function(object) {
    TRUE
})

##' Show Method for DiagonalCov
##' 
##' Prints a summary of the DiagonalCov Object.
##' @param object A `DiagonalCov` object. 
##' @export
setMethod("show", "DiagonalCov" function(object) { 
    	cat(class(object), "object(dimension:", object@dimension, ")\n")
	cat(" Variances:", round(object@internal_diag_values, 4), "\n")
	cat(" Unconstrained parameters:", round(object@parameters, 4), "\n")
})


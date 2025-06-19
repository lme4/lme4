
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
setGeneric("compute_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_covariance_matrix"))

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
setGeneric("compute_inverse_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_inverse_covariance_matrix"))

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

##' Get Lower Bounds for Parameters
##'
##' Generic function to retrieve the vector of lower bounds for the
##' unconstrained optimization parameters.
##' @param object An S4 object inheriting from `VirtualCovariance`.
##' @return A numeric vector of lower bounds for each parameter.
##' @export
setGeneric("get_lower_bounds", function(object) standardGeneric("get_lower_bounds"))

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
	contains = c("VirtualCovariance", "VirtualParameterizationLogScaleBoundedCor")
)

setValidity("DiagonalCov", function(object) {
	if (length(object@parameters) != object@dimension) {
        return(sprintf("Expected %d parameters, but slot has %d.", object@dimension, length(object@parameters)))
	}
	TRUE
})

# Method for DiagonalCov 

##' @rdname n_parameters
setMethod("n_parameters", "DiagonalCov", function(object) {
	object@dimension
})

##' @rdname get_parameters
setMethod("get_parameters", "DiagonalCov", function(object) {
	object@parameters
})

##' @rdname set_parameters
setMethod("set_parameters", "DiagonalCov", function(object, value) {
	if (length(value) != object@dimension) {
		stop("Incorrect number of parameters for DiagonalCov. Expected ", object@dimension)
	}
	object@parameters <- value 
	object
})

##' @rdname compute_covariance_matrix
setMethod("compute_covariance_matrix", "DiagonalCov", function (object, data_context = NULL) {
	variances <- exp(object@parameters)
	Matrix::Diagonal(x = variances)
})
					   

##' @rdname compute_log_det_covariance_matrix
setMethod("compute_log_det_covariance_matrix", "DiagonalCov", function(object, data_context = NULL) {
	sum(object@parameters)
})

##' @rdname compute_inverse_covariance_matrix
setMethod("compute_inverse_covariance_matrix", "DiagonalCov", function(object, data_context = NULL) {
	inv_variances <- exp(-object@parameters)
	Matrix::Diagonal(x = inv_variances)
})

##' @rdname get_cholesky_factor
setMethod("get_cholesky_factor", "DiagonalCov", function(object, data_context = NULL) {
	std_devs <- exp(0.5 * object@parameters)
	Matrix::Diagonal(x = std_devs)
})

##' @rdname get_start_values
setMethod("get_start_values", "DiagonalCov", function(object) {
	rep(log(1), object@dimension) # Start with log(1) for all variances
})

##' @rdname get_interpretable_parameters
setMethod("get_interpretable_parameters", "DiagonalCov", function(object) {
	list(variances = exp(object@parameters))
})

##' @rdname is_diagonal
setMethod("is_diagonal", "DiagonalCov", function(object) {
	TRUE
})

##' @rdname get_lower_bounds
setMethod("get_lower_bounds", "DiagonalCov", function(object) {
	rep(-Inf, n_parameters(object))
})

##' Show Method for DiagonalCov
##' 
##' Prints a summary of the DiagonalCov Object.
##' @param object A `DiagonalCov` object. 
##' @export
setMethod("show", "DiagonalCov", function(object) {
	cat(class(object), " object (dimension:", object@dimension, ")\n")
	try({
        	variances <- exp(object@parameters)
        	cat("  Variances:", round(variances, 4), "\n")
	}, silent = TRUE)
})

##' Unstructured Covariance Structure
##'
##' Represents a general, unstructured positive semi-definite covariance matrix,
##' typically parameterized using a log-Cholesky decomposition.
##'
##' @slot internal_L A dtrMatrix storing the lower-triangular Cholesky factor (L)
##' such that Sigma = L %*% t(L).
##' @export
setClass("UnstructuredCov",
	contains = c("VirtualCovariance", "VirtualParameterizationLogChol")
)

setValidity("UnstructuredCov", function(object) {
	expected_len <- (object@dimension * (object@dimension + 1)) / 2
	if (length(object@parameters) != expected_len) {
        	return(sprintf("Expected %d parameters, but slot has %d.",
			       expected_len,
			       length(object@parameters)
    			)
		)
	} 
	TRUE
})
    

# Method for UnstructuredCov

##' @rdname n_parameters
setMethod("n_parameters", "UnstructuredCov", function(object) {
	d <- object@dimension 
	d * (d + 1) / 2
})

##' @rdname get_parameters
setMethod("get_parameters", "UnstructuredCov", function(object) {
	object@parameters
})

##' @rdname set_parameters
setMethod("set_parameters",
	signature(object = "UnstructuredCov", value = "numeric"),
	function(object, value) {
        	if (length(value) != n_parameters(object)) {
            stop("Incorrect number of parameters for UnstructuredCov.")
        }
        object@parameters <- value
        object
	}
)

##' @rdname set_parameters
setMethod("set_parameters",
    signature(object = "UnstructuredCov", value = "matrix"),
    function(object, value) {
        d <- object@dimension
        if (!all(dim(value) == d)) {
          stop(sprintf("Input matrix must have dimension %d x %d.", d, d))
        }
        L_vec <- value[lower.tri(value, diag = TRUE)]
        diag_indices <- vech_diag_indices(d)
        L_vec[diag_indices] <- log(L_vec[diag_indices])

        object@parameters <- L_vec
        object
    }
)

##' @rdname compute_covariance_matrix
setMethod("compute_covariance_matrix", "UnstructuredCov", function(object, data_context = NULL) {
	L <- get_chol_from_params(object)
	tcrossprod(L)
})

##' @rdname compute_log_det_covariance_matrix
setMethod("compute_log_det_covariance_matrix", "UnstructuredCov", function(object, data_context = NULL) {
	diag_indices <- vech_diag_indices(object@dimension)
	2 * sum(object@parameters[diag_indices])
})

#' @rdname compute_inverse_covariance_matrix
setMethod("compute_inverse_covariance_matrix", "UnstructuredCov", function(object, data_context = NULL) {
	L <- get_chol_from_params(object)
	inv_L <- solve(L)
	crossprod(inv_L)
})

##' @rdname get_cholesky_factor
setMethod("get_cholesky_factor", "UnstructuredCov", function(object, data_context = NULL) {
	get_chol_from_params(object)
})

##' @rdname get_start_values

setMethod("get_start_values", "UnstructuredCov", function(object) {
	rep(0, n_parameters(object)) # Assuming 0 for off-diagonals and 0 for log(1) diagonals
})

##' @rdname get_interpretable_parameters
setMethod("get_interpretable_parameters", "UnstructuredCov", function(object) {
	Sigma <- compute_covariance_matrix(object)
	st_devs <- sqrt(diag(Sigma))

	list(
		st_devs = st_devs,
		covariance_matrix = Sigma	
	)
})

##' @rdname is_diagonal
setMethod("is_diagonal", "UnstructuredCov", function(object) {

	FALSE # By definition
})

##' @rdname get_lower_bounds
setMethod("get_lower_bounds", "UnstructuredCov", function(object) {
    	# The log-Cholesky parameterization is unconstrained.	
    	rep(-Inf, n_parameters(object))
})

##' Show Method for UnstructuredCov
##'
##' Prints a summary of the UnstructuredCov object.
##' @param object A `UnstructuredCov` object.
##' @export
    setMethod("show", "UnstructuredCov", function(object) {
    cat(sprintf("%s object (dimension: %d)\n", class(object), object@dimension))
    try({
        params <- get_interpretable_parameters(object)
        cat("  Standard Deviations:", round(params$st_devs, 4), "\n")
        if (object@dimension > 1) {
            cat("  Covariance Matrix:\n")
            print(round(params$covariance_matrix, 4))
        }
    }, silent = TRUE)
})  


# Compound Symmetry Covariance Structure 

##' Compound Symmetry Covariance Structure
##'
##' Represents a compound symmetric (or exchangeable) covariance matrix,
##' characterized by a single common variance and a single common correlation
##' for all pairs of elements. This structure is parameterized by two
##' unconstrained values: the log of the variance and the inverse hyperbolic
##' tangent (`atanh`) of the correlation.
##'
##' @export
setClass("CompoundSymmetryCov",
    contains = c("VirtualCovariance", "VirtualParameterizationLogScaleBoundedCor")
)
setValidity("CompoundSymmetryCov", function(object) {
    if(length(object@parameters) != 2) {
        return("Expected 2 parameters (log-variance, atanh-correlation).")
    }

    d <- object@dimension
    if (d > 1) {
        rho <- tanh(object@parameters[2])
        if (rho <= -1 / (d - 1)) {
            return(sprintf("Implied correlation (rho) of %f is not valid for dimension %d. Must be >%f.", 
                           rho, d, -1 / (d - 1)))
        }
    }
    TRUE
})


##' @rdname n_parameters 
setMethod("n_parameters", "CompoundSymmetryCov", function(object) {
    2L 
})

##' @rdname get_parameters
setMethod("get_parameters", "CompoundSymmetryCov", function(object) {
    object@parameters
})

##' @rdname set_parameters
setMethod("set_parameters", "CompoundSymmetryCov", function(object,value) {
    if (length(value) != 2) {
        stop("Incorrect number of parameters for CompoundSymmetryCov. Expected 2.")
    }
    object@parameters <- value 
    validObject(object)
    object
})

##' @rdname get_interpretable_parameters
setMethod("get_interpretable_parameters", "CompoundSymmetryCov", function(object, data_context = NULL) {
    get_cs_interpretable_params(object)
})

##' @rdname compute_covariance_matrix
setMethod("compute_covariance_matrix", "CompoundSymmetryCov", function(object, data_context = NULL) {    
    d <- object@dimension
    if (d == 0) return (new("dsyMatrix", Dim = c(0L, 0L)))

    params <- get_cs_interpretable_params(object)
    sigma_sq <- params$variance
    rho <- params$correlation 

    I <- Diagonal(d)
    J <- Matrix(1, nrow = d, ncol = d)
    
    Sigma <- sigma_sq * ((1 - rho) * I + rho * J)
    Sigma 
  })

##' @rdname get_cholesky_factor
setMethod("get_cholesky_factor", "CompoundSymmetryCov", function(object, data_context = NULL) {
    chol(compute_covariance_matrix(object))
})


##' @rdname compute_log_det_covariance_matrix 
setMethod("compute_log_det_covariance_matrix", "CompoundSymmetryCov", function(object, data_context = NULL) {
    d <- object@dimension
    if (d == 0) return(0)

    params <- get_cs_interpretable_params(object)

    # See: Searle, S. R. (1982). Matrix Algebra Useful for Statistics. Wiley. (Specifically, in sections on patterned matrices).
    # det(Sigma) = [sigma^2(1 - rho)]^(d-1) * [sigma^2(1 + (d - 1)rho)]
    # Formula: log((sigma_sq*(1-rho))^(d-1) * (sigma_sq*(1+(d-1)*rho)))

    (d - 1) * log(params$variance * (1 - params$correlation)) + log(params$variance * (1 + (d - 1) * params$correlation))
})

##' @rdname compute_inverse_covariance_matrix 



# S4 Class and Method Definition for Covariance Structures

# SECTION 1. TOP-LEVEL VIRTUAL CLASSES AND GENERICS

##' @title Virtual Class for All Covariance Structures
##' @description This abstract S4 class defines the common interface for all covariance models.
##' @slot dimension An integer indicating the size of the square covariance matrix.
##' @slot vparameters A numeric vector for variance-related parameters.
##' @slot cparameters A numeric vector for covariance-related parameters.
#3' @keywords internal 
setClass("VirtualCovariance",
    contains = "VIRTUAL",
    slots = c(
        dimension = "integer",
        vparameters = "numeric",
        cparameters = "numeric"
    )
)
##' @rdname VirtualComponentClasses
setValidity("VirtualCovariance", function(object) {
    if (object@dimension < 1L) {
        return("Covariance dimension must be at least 1 (cannot have zero random effects)")
    }
    TRUE
})

##' @name CovarianceMethods
##' @title Public S4 Methods for Covariance Objects
##'
##' @description A collection of S4 generic functions to interact with covariance
##'   structure objects.
##'
##' @param object A covariance structure object.
##' @param value A numeric vector of parameters to set.
##' @param data_context A context object, currently unused (for future expansion).
##'
##' @details
##' The following methods form the main public API for these objects:
##' - `n_parameters()`: Get the total number of parameters.
##' - `get_parameters()`: Get the complete vector of current parameters.
##' - `set_parameters()`: Set the parameters from a numeric vector.
##' - `get_interpretable_parameters()`: Get a list of human-readable parameters
##'   (e.g., standard deviations and correlations).
##' - `compute_covariance_matrix()`: Compute the final covariance matrix.
##' - `compute_cholesky_factor()`: Compute the Cholesky factor of the covariance matrix.
##' - `compute_inverse_covariance_matrix()`: Compute the inverse of the
##'   covariance matrix.
##' - `compute_log_det_covariance_matrix()`: Compute the log-determinant of the
##'   covariance matrix.
NULL

##' @rdname CovarianceMethods
##' @export
setGeneric("n_parameters", function(object) standardGeneric("n_parameters"))
##' @rdname CovarianceMethods
##' @export
setGeneric("get_parameters", function(object) standardGeneric("get_parameters"))
##' @rdname CovarianceMethods
##' @export
setGeneric("set_parameters", function(object, value) standardGeneric("set_parameters"))
##' @rdname CovarianceMethods
##' @export
setGeneric("get_interpretable_parameters", function(object) standardGeneric("get_interpretable_parameters"))
##' @rdname CovarianceMethods
##' @export
setGeneric("compute_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_covariance_matrix"))
##' @rdname CovarianceMethods
##' @export
setGeneric("compute_cholesky_factor", function(object, data_context = NULL) standardGeneric("compute_cholesky_factor"))
##' @rdname CovarianceMethods
##' @export
setGeneric("compute_inverse_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_inverse_covariance_matrix"))
##' @rdname CovarianceMethods
##' @export
setGeneric("compute_log_det_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_log_det_covariance_matrix"))
##' @name InternalCovarianceMethods
##' @title Internal S4 Generics for Covariance Objects
##' @description These generic functions are intended for internal package use
##'  and are not part of the public API.
##' @param object A covariance structure object.
##' @keywords internal
NULL

##' @rdname InternalCovarianceMethods
setGeneric("get_start_values", function(object) standardGeneric("get_start_values"))
##' @rdname InternalCovarianceMethods
setGeneric("get_lower_bounds", function(object) standardGeneric("get_lower_bounds"))
##' @rdname InternalCovarianceMethods
setGeneric("is_diagonal", function(object) standardGeneric("is_diagonal"))
##' @rdname InternalCovarianceMethods
setGeneric("get_lambda", function(object) standardGeneric("get_lambda"))
##' @rdname InternalCovarianceMethods
setGeneric("get_lind", function(object) standardGeneric("get_lind"))
##' @rdname InternalCovarianceMethods
setGeneric("compute_correlation_matrix", function(object) standardGeneric("compute_correlation_matrix"))
##' @rdname InternalCovarianceMethods
setGeneric("compute_lambdat_x", function(object, theta) standardGeneric("compute_lambdat_x"))


# SECTION 2: COMPONENT CLASS HIERARCHY & VALIDITY

##' @title Virtual Classes for Covariance Structure Components
##' @description These virtual classes define the building blocks for the
##'   concrete covariance models and are intended for internal use only.
##'                  Unstructured     Diagonal     CS ...
##' Homogeneous      HomUns...        HomDiag.    HomCS..
##' Heterogeneous    HetUns...        HetDiag.    HetCS..
##' VirtualCovariance is a virtual superclass of everything else.
##' Row and column classes are virtual subclasses of VirtualCovariance.
##' Matrix elements are nonvirtual subclasses of the corresponding row
##' row and column classes.
##' @keywords internal
##' @name VirtualComponentClasses
NULL

##' @rdname VirtualComponentClasses
setClass("HomogeneousVariance", contains = c("VIRTUAL", "VirtualCovariance"))
##' @rdname VirtualComponentClasses
setClass("HeterogeneousVariance", contains = c("VIRTUAL", "VirtualCovariance"))

##' @rdname VirtualComponentClasses
setClass("DiagonalCovariance", contains = c("VIRTUAL", "VirtualCovariance"))
setValidity("DiagonalCovariance", function(object) {
    if (length(object@cparameters) == 0L) TRUE else "Diagonal models must have 0 cparameters."
})

##' @title Concrete Class for Unstructured Covariance
##' @description Defines the Unstructured structure, which uses all elements
##'   of the lower triangle as free parameters.
##' @keywords internal
setClass("UnstructuredCovariance", contains = "VirtualCovariance")
setValidity("UnstructuredCovariance", function(object) {
    d <- object@dimension
    expected_len <- d * (d + 1) / 2
    if (length(object@cparameters) == expected_len) TRUE else "Incorrect number of cparameters for Unstructured."
})


##' @title Virtual Class for Compound Symmetry (CS) Covariance
##' @description Defines the CS structure, which uses 2 `cparameters` for the
##'   diagonal and off-diagonal correlation elements.
##' @keywords internal
setClass("CSCovariance", contains = c("VIRTUAL", "VirtualCovariance"))
setValidity("CSCovariance", function(object) {
    d <- object@dimension
    expected_len <- if (d > 1) 1L else 0L
    if (length(object@cparameters) == expected_len) TRUE else "CS models must have 1 cparameters (diag and off-diag)."
})

##' @title Virtual Class for Autoregressive Order 1 (AR1) Covariance
##' @description Defines the AR1 structure, which uses `d` `cparameters`,
##'   one for each correlation lag.
##' @keywords internal
setClass("AR1Covariance", contains = c("VIRTUAL", "VirtualCovariance"))
setValidity("AR1Covariance", function(object) {
    d <- object@dimension 
    expected_len <- if (d > 1) 1L else 0L
    if (length(object@cparameters) == expected_len) TRUE else "AR1 models must have 1 cparameter."
})

##' @title Concrete Covariance Structure Classes
##' @description These are the user-facing S4 classes for specifying covariance
##' structures. They can be initialized with `new("ClassName", dimension = d)`.
##'
##' @name ConcreteCovarianceClasses
##' @export HomogeneousDiagonalCovariance
##' @export HeterogeneousDiagonalCovariance
##' @export HomogeneousCSCovariance
##' @export HeterogeneousCSCovariance
##' @export HomogeneousAR1Covariance
##' @export HeterogeneousAR1Covariance
##' @export UnstructuredCovariance
NULL

# Create the Concrete Classes
for (variance_type in c("Homogeneous", "Heterogeneous")) {
    for (structure_type in c("Diagonal", "CS", "AR1")) {
        setClass(
            paste0(variance_type, structure_type, "Covariance"),
            contains = c(
                paste0(variance_type, "Variance"),
                paste0(structure_type, "Covariance")
            )
        )
    }
}

# SECTION 3: PARAMETER MANAGEMENT METHODS

##' @title Initialize a Covariance Object
##' @description This method ensures that when a new object is created, it is
##' immediately populated with valid default starting parameters. Not for direct use.
##' @keywords internal
setMethod("initialize", "VirtualCovariance",    
  function(.Object, ...) { 
    args <- list(...)

    if (is.null(args$vparameters) && is.null(args$cparameters) && !is.null(args$dimension)) {
      
      temp_obj <- .Object 
      temp_obj@dimension <- as.integer(args$dimension)

      if (is(temp_obj, "UnstructuredCovariance")) {  
          args$vparameters <- numeric(0)
          args$cparameters <- get_start_values(temp_obj)
      } else {
          start_params <- get_start_values(temp_obj)
          n_v_params <- if (is(temp_obj, "HomogeneousVariance")) 1L else temp_obj@dimension
          
          args$vparameters <- start_params[seq_len(n_v_params)]
          args$cparameters <- start_params[-(seq_len(n_v_params))]
      }
    } 

    do.call(callNextMethod, c(list(.Object), args))
})

# The number of parameters for a given covariance structure is defined as the number
# of unique, estimated values needed to define its cholesky factor Lambda. 

##' @rdname CovarianceMethods
setMethod("n_parameters", "UnstructuredCovariance", function(object) {
    d <- object@dimension 
    as.integer(d * (d + 1) / 2)
})

##' @rdname CovarianceMethods
setMethod("n_parameters", "DiagonalCovariance", function(object) {
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else object@dimension
    as.integer(n_v_params) # No correlation parameters
})

##' @rdname CovarianceMethods
setMethod("n_parameters", "CSCovariance", function(object) {
    d <- object@dimension
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    n_c_params <- if (d > 1) 1L else 0L
    as.integer(n_v_params + n_c_params)
})

##' @rdname CovarianceMethods
setMethod("n_parameters", "AR1Covariance", function(object) {
    d <- object@dimension
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    n_c_params <- if (d > 1) 1L else 0L
    as.integer(n_v_params + n_c_params)
})
##' @rdname CovarianceMethods
setMethod("set_parameters", "VirtualCovariance", function(object, value) {
    n_total <- n_parameters(object)
    if (length(value) != n_total) {
        stop(sprintf("Incorrect number of parameters. Expected %d, got %d.", n_total, length(value)))
    }
    
    if (is(object, "UnstructuredCovariance")) {
        object@vparameters <- numeric(0)
        object@cparameters <- value
    } else {
        d <- object@dimension
        if (is(object, "HomogeneousVariance")) {
            n_v_params <- 1L
        } else {
            n_v_params <- max(d, 1L)
        }
        
        if (n_v_params <= length(value)) {
            object@vparameters <- value[seq_len(n_v_params)]
            if (length(value) > n_v_params) {
                object@cparameters <- value[-seq_len(n_v_params)]
            } else {
                object@cparameters <- numeric(0)
            }
        } else {
            object@vparameters <- value
            object@cparameters <- numeric(0)
        }
    }

    validObject(object)
    return(object)
})

##' @rdname CovarianceMethods
setMethod("get_start_values", "UnstructuredCovariance", function(object) {
    d <- object@dimension
    num_params <- d * (d + 1) / 2
    start_vals <- numeric(num_params)
    
    diag_indices <- vech_diag_indices(d)
    start_vals[diag_indices] <- 1.0
    
    return(start_vals)
})

##' @rdname CovarianceMethods
setMethod("get_start_values", "DiagonalCovariance", function(object) {
    num_params <- n_parameters(object)
    rep(1.0, num_params) 
})

##' @rdname CovarianceMethods
setMethod("get_start_values", "CSCovariance", function(object) {
    get_structured_start_values(object)
})

##' @rdname CovarianceMethods
setMethod("get_start_values", "AR1Covariance", function(object) {
    get_structured_start_values(object)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lower_bounds", "UnstructuredCovariance", function(object) {
    d <- object@dimension
    
    num_params <- d * (d + 1) / 2
    lower_bounds <- rep(-Inf, num_params)
    
    diag_indices <- vech_diag_indices(d)
    lower_bounds[diag_indices] <- 0
    
    return(lower_bounds)
})

setMethod("get_lower_bounds", "DiagonalCovariance", function(object) {
    d <- object@dimension
    if (is(object, "HomogeneousVariance")) {
        return(0) 
    } else {
        return(rep(0, d))
    }
})

##' @rdname InternalCovarianceMethods
setMethod("get_lower_bounds", "CSCovariance", function(object) {
    get_structured_lower_bounds(object)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lower_bounds", "AR1Covariance", function(object) {
    get_structured_lower_bounds(object)
})


##' @rdname CovarianceMethods
setMethod("get_parameters", "VirtualCovariance", function(object) c(object@vparameters, object@cparameters))
##' @rdname InternalCovarianceMethods
setMethod("is_diagonal", "DiagonalCovariance", function(object) TRUE)
##' @rdname InternalCovarianceMethods
setMethod("is_diagonal", "VirtualCovariance", function(object) FALSE)

# SECTION 4: MODEL FITTING INTEGRATION METHODS (GET_LAMBDA, GET_LIND, COMPUTE_LAMBDAT_X)

##' @rdname InternalCovarianceMethods
setMethod("get_lambda", "UnstructuredCovariance", function(object) {
    d <- object@dimension
    if (d == 0) return(new("dgTMatrix"))
    indices <- get_lower_tri_indices(d)
    sparseMatrix(i = indices$i, j = indices$j, x = 1.0, dims = c(d, d))
})

##' @rdname InternalCovarianceMethods
setMethod("get_lambda", "DiagonalCovariance", function(object) {
    d <- object@dimension
    if (d == 0) return(new("dgTMatrix"))
    sparseMatrix(i = 1:d, j = 1:d, x = 1.0, dims = c(d, d))
})

##' @rdname InternalCovarianceMethods
setMethod("get_lambda", "CSCovariance", function(object) {
    d <- object@dimension
    if (d == 0) return(new("dgTMatrix"))
    indices <- get_lower_tri_indices(d)
    sparseMatrix(i = indices$i, j = indices$j, x = 1.0, dims = c(d, d))
})

##' @rdname InternalCovarianceMethods
setMethod("get_lambda", "AR1Covariance", function(object) {
    d <- object@dimension
    if (d == 0) return(new("dgTMatrix"))
    indices <- get_lower_tri_indices(d)
    sparseMatrix(i = indices$i, j = indices$j, x = 1.0, dims = c(d, d))
})

###' @rdname InternalCovarianceMethods
setMethod("get_lind", "UnstructuredCovariance", function(object) {
    num_params <- n_parameters(object)
    if (num_params == 0) return(integer(0))
    
    seq_len(num_params)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lind", "HomogeneousDiagonalCovariance", function(object) {
    rep(1L, object@dimension)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lind", "HeterogeneousDiagonalCovariance", function(object) {
    seq_len(object@dimension)
})

setMethod("get_lind", "HomogeneousCSCovariance", function(object) {
    d <- object@dimension
    if (d < 2) return(integer(0))
    
    n_ltri <- d * (d + 1) / 2
    res <- integer(n_ltri)  
    diag_indices <- vech_diag_indices(d) 
    
    res[diag_indices] <- 1L
    res[-diag_indices] <- 2L
    
    return(res)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lind", "HeterogeneousCSCovariance", function(object) {
    d <- object@dimension
    if (d < 2) return(integer(0))
    n_v_params <- d 
    n_ltri <- d * (d + 1) / 2
    res <- integer(n_ltri)
    diag_indices <- vech_diag_indices(d) 
    res[diag_indices] <- seq_len(n_v_params)
    res[-diag_indices] <- n_v_params + 1L
    
    return(res)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lind", "HomogeneousAR1Covariance", function(object) {
    d <- object@dimension
    if (d < 2) return(integer(0))
    n_ltri <- d * (d + 1) / 2
    res <- integer(n_ltri)
    diag_indices <- vech_diag_indices(d)
    res[diag_indices] <- 1L  
    res[-diag_indices] <- 2L 
    return(res)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lind", "HeterogeneousAR1Covariance", function(object) {
    d <- object@dimension
    if (d < 2) return(integer(0))
    n_v_params <- d
    n_ltri <- d * (d + 1) / 2
    res <- integer(n_ltri)
    diag_indices <- vech_diag_indices(d)
    res[diag_indices] <- seq_len(n_v_params) 
    res[-diag_indices] <- n_v_params + 1L    
    return(res)
})

##' @rdname InternalCovarianceMethods
setMethod("compute_lambdat_x", "DiagonalCovariance", function(object, theta) {
    return(exp(theta))
})

##' @rdname InternalCovarianceMethods
setMethod("compute_lambdat_x", "UnstructuredCovariance", function(object, theta) {
    d <- object@dimension
    if (d == 0) return(numeric(0))
    
    L <- get_chol_from_params(theta, d)
    return(L[lower.tri(L, diag = TRUE)])
})

##' @rdname InternalCovarianceMethods
setMethod("compute_lambdat_x", "HomogeneousCSCovariance", function(object, theta) {
    d <- object@dimension
    if (d == 0) return(numeric(0))

    sigma_re <- exp(0.5 * theta[1])
    rho <- tanh(theta[2])

    R <- matrix(rho, d, d)
    diag(R) <- 1.0
    L_corr <- t(chol(R))
    L <- sigma_re * L_corr

    return(L[lower.tri(L, diag = TRUE)])
})

##' @rdname InternalCovarianceMethods
setMethod("compute_lambdat_x", "HeterogeneousCSCovariance", function(object, theta) {
    d <- object@dimension
    if (d == 0) return(numeric(0))
    if (d == 1) return(exp(0.5 * theta[1]))

    st_devs <- exp(0.5 * theta[1:d])
    rho <- tanh(theta[d + 1])

    R <- matrix(rho, d, d)
    diag(R) <- 1.0
    L_corr <- t(chol(R)) 
    
    D <- Diagonal(x = st_devs)
    L <- D %*% L_corr

    return(L[lower.tri(L, diag = TRUE)])
})


##' @rdname InternalCovarianceMethods
setMethod("compute_lambdat_x", "HomogeneousAR1Covariance", function(object, theta) {
    d <- object@dimension
    if (d == 0) return(numeric(0))

    sigma_re <- exp(0.5 * theta[1])
    rho <- tanh(theta[2])

    L <- matrix(0, d, d)
    
    first_col <- sigma_re * (rho^(0:(d - 1)))
    L[, 1] <- first_col
    
    if (d > 1) {
        s <- sqrt(1 - rho^2)
        for (j in 2:d) {
            L[j:d, j] <- s * first_col[1:(d - j + 1)]
        }
    }
    
    return(L[lower.tri(L, diag = TRUE)])
})

##' @rdname InternalCovarianceMethods
setMethod("compute_lambdat_x", "HeterogeneousAR1Covariance", function(object, theta) {
    d <- object@dimension
    if (d == 0) return(numeric(0))
    if (d == 1) return(exp(0.5 * theta[1]))

    st_devs <- exp(0.5 * theta[1:d])
    rho <- tanh(theta[d + 1])

    L_corr <- matrix(0, d, d)
    first_col_corr <- rho^(0:(d - 1))
    L_corr[, 1] <- first_col_corr
    if (d > 1) {
        s <- sqrt(1 - rho^2)
        for (j in 2:d) {
            L_corr[j:d, j] <- s * first_col_corr[1:(d - j + 1)]
        }
    }
    
    D <- Diagonal(x = st_devs)
    L <- D %*% L_corr
    
    return(L[lower.tri(L, diag = TRUE)])
})
          
# SECTION 5: POST-HOC COMPUTATION & INTERPRETATION METHODS

##' @rdname VirtualComponentClasses
setClass("VirtualStructuredCovariance", contains = c("VIRTUAL", "VirtualCovariance"))
setClass("DiagonalCovariance", contains = c("VIRTUAL", "VirtualCovariance", "VirtualStructuredCovariance"))
setClass("CSCovariance", contains = c("VIRTUAL", "VirtualCovariance", "VirtualStructuredCovariance"))
setClass("AR1Covariance", contains = c("VIRTUAL", "VirtualCovariance", "VirtualStructuredCovariance"))

##' @rdname InternalCovarianceMethods
setMethod("compute_correlation_matrix", "DiagonalCovariance", function(object) {
    d <- object@dimension 
    R <- Diagonal(d)

    force_dsyMatrix(R)
})
##' @rdname InternalCovarianceMethods
setMethod("compute_correlation_matrix", "CSCovariance", function(object) {
    d <- object@dimension
    if (d == 1) return(Diagonal(1))
    rho <- tanh(object@cparameters[1])
    R <- Matrix(rho, nrow = d, ncol = d)
    diag(R) <- 1
    
    force_dsyMatrix(R)
})
##' @rdname InternalCovarianceMethods
setMethod("compute_correlation_matrix", "AR1Covariance", function(object) {
    d <- object@dimension
    if (d == 1) return(Diagonal(d))
    rho <- tanh(object@cparameters[1])
    time_diffs <- abs(outer(1:d, 1:d, "-"))
    R <- rho^time_diffs

    force_dsyMatrix(R)
})
##' @rdname InternalCovarianceMethods
setMethod("compute_correlation_matrix", "UnstructuredCovariance", function(object) {
    Sigma <- compute_covariance_matrix(object)
    R <- cov2cor(Sigma)

    force_dsyMatrix(R)
})

##' @rdname CovarianceMethods
setMethod("compute_covariance_matrix", "DiagonalCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    
    if (is(object, "HomogeneousVariance")) {
        sigma_sq <- exp(object@vparameters[1])
        Sigma <- Diagonal(d, x = sigma_sq)
    } else {
        variances <- exp(object@vparameters)
        Sigma <- Diagonal(d, x = variances)
    }

    force_dsyMatrix(Sigma)
})

##' @rdname CovarianceMethods
setMethod("compute_covariance_matrix", "CSCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    
    R <- compute_correlation_matrix(object)
    
    if (is(object, "HomogeneousVariance")) {
        sigma_sq <- exp(object@vparameters[1])
        Sigma <- sigma_sq * R
    } else {
        st_devs <- exp(0.5 * object@vparameters)
        D <- Diagonal(d, x = st_devs)
        Sigma <- D %*% R %*% D
    }
    
    force_dsyMatrix(Sigma)
})

##' @rdname CovarianceMethods
setMethod("compute_covariance_matrix", "AR1Covariance", function(object, data_context = NULL) {
    d <- object@dimension
    
    R <- compute_correlation_matrix(object)

    if (is(object, "HomogeneousVariance")) {
        sigma_sq <- exp(object@vparameters[1])
        Sigma <- sigma_sq * R
    } else {
        st_devs <- exp(0.5 * object@vparameters)
        D <- Diagonal(d, x = st_devs)
        Sigma <- D %*% R %*% D
    }
    
    force_dsyMatrix(Sigma)
})

##' @rdname CovarianceMethods
setMethod("compute_covariance_matrix", "UnstructuredCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    
    L <- get_chol_from_params(object@cparameters, d)
    
    Sigma <- tcrossprod(L)
    force_dsyMatrix(Sigma)
})

##' @rdname CovarianceMethods
setMethod("compute_cholesky_factor", "VirtualCovariance", function(object, data_context = NULL) {
    Sigma <- compute_covariance_matrix(object, data_context)
    t(chol(Sigma))
})
##' @rdname CovarianceMethods
setMethod("compute_cholesky_factor", "UnstructuredCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    L <- (get_chol_from_params(object@cparameters, d))
    as.matrix(L)
})

##' @rdname CovarianceMethods
setMethod("compute_log_det_covariance_matrix", "DiagonalCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    if (d == 0) return(0)
    
    if (is(object, "HomogeneousVariance")) {
        d * object@vparameters[1]
    } else {
        sum(object@vparameters)
    }
})

##' @rdname CovarianceMethods
setMethod("compute_log_det_covariance_matrix", "CSCovariance", function(object, data_context = NULL) {
    compute_log_det_structured(object)
})

##' @rdname CovarianceMethods
setMethod("compute_log_det_covariance_matrix", "AR1Covariance", function(object, data_context = NULL) {
    compute_log_det_structured(object)
})

##' @rdname CovarianceMethods
setMethod("compute_log_det_covariance_matrix", "UnstructuredCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    if (d == 0) return(0)

    L <- compute_cholesky_factor(object)
    2 * sum(log(base::diag(L)))
})

##' @rdname CovarianceMethods
setMethod("compute_inverse_covariance_matrix", "DiagonalCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
    
    if (is(object, "HomogeneousVariance")) {
        inv_sigma_sq <- exp(-object@vparameters[1])
        inv_Sigma <- Diagonal(d, x = inv_sigma_sq)
    } else {
        inv_variances <- exp(-object@vparameters)
        inv_Sigma <- Diagonal(d, x = inv_variances)
    }
    
    # Ensure consistent matrix type
    force_dsyMatrix(inv_Sigma)
})

##' @rdname CovarianceMethods
setMethod("compute_inverse_covariance_matrix", "CSCovariance", function(object, data_context = NULL) {
    compute_inverse_structured(object)
})

##' @rdname CovarianceMethods
setMethod("compute_inverse_covariance_matrix", "AR1Covariance", function(object, data_context = NULL) {
    compute_inverse_structured(object)
})

##' @rdname CovarianceMethods
setMethod("compute_inverse_covariance_matrix", "UnstructuredCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
    
    # inv(Sigma) = inv(L * L') = crossprod(solve(L))
    L <- compute_cholesky_factor(object)
    inv_L <- solve(L)
    inv_Sigma <- crossprod(inv_L)
    
    force_dsyMatrix(inv_Sigma)
})



##' @rdname CovarianceMethods
setMethod("get_interpretable_parameters", "UnstructuredCovariance", function(object) {
    Sigma <- compute_covariance_matrix(object)
    st_devs <- sqrt(diag(as.matrix(Sigma)))
    cor_matrix <- cov2cor(Sigma)
    
    params <- list(st_devs = st_devs, correlation = cor_matrix)
    return(params)
})

##' @rdname CovarianceMethods
setMethod("get_interpretable_parameters", "DiagonalCovariance", function(object) {
    params <- list()
    if (is(object, "HomogeneousVariance")) {
        params$st_dev <- exp(0.5 * object@vparameters[1])
    } else {
        params$st_devs <- exp(0.5 * object@vparameters)
    }
    return(params)
})

##' @rdname CovarianceMethods
setMethod("get_interpretable_parameters", "CSCovariance", function(object) {
    params <- list()
    if (is(object, "HomogeneousVariance")) {
        params$st_dev <- exp(0.5 * object@vparameters[1])
    } else {
        params$st_devs <- exp(0.5 * object@vparameters)
    }
    if (length(object@cparameters) > 0) {
        params$correlation <- tanh(object@cparameters[1])
    }
    
    return(params)
})
##' @rdname CovarianceMethods
setMethod("get_interpretable_parameters", "AR1Covariance", function(object) {
    params <- list()
    if (is(object, "HomogeneousVariance")) {
        params$st_dev <- exp(0.5 * object@vparameters[1])
    } else {
        params$st_devs <- exp(0.5 * object@vparameters)
    }    
    if (length(object@cparameters) > 0) {
        params$correlation <- tanh(object@cparameters[1])
    }
    
    return(params)
})
# SECTION 6: SHOW METHOD
setMethod("show", "VirtualCovariance", function(object) {
    model_name <- class(object)
    n_params <- n_parameters(object)
    cat(sprintf("%s object (dimension: %d, parameters: %d)\n",
        model_name, object@dimension, n_params))
    try({
        params <- get_interpretable_parameters(object)
        if (length(params) > 0) {
            cat("  Interpretable Parameters:\n")
            if (!is.null(params$st_devs)) {
                st_devs_formatted <- paste(sprintf("%.4f", params$st_devs), collapse = " ")
                cat("    Standard Deviations:", st_devs_formatted, "\n")
            } else if (!is.null(params$st_dev)) {
                cat(sprintf("    Standard Deviation (sigma): %.4f\n", params$st_dev))
            }
            if (!is.null(params$correlation)) {
                cat(sprintf("    Correlation (rho): %.4f\n", params$correlation))
            } else if (!is.null(params$correlation)) {
                cat(sprintf("    Lag-1 Correlation (rho): %.4f\n", params$correlation))
            }
        }
    }, silent = TRUE)
})


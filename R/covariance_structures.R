# S4 Class and Method Definition for Covariance Structures

library(Matrix)


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
##' - `get_cholesky_factor()`: Compute the Cholesky factor of the covariance matrix.
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
setGeneric("get_cholesky_factor", function(object, data_context = NULL) standardGeneric("get_cholesky_factor"))
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
setGeneric("build_correlation_matrix", function(object) standardGeneric("build_correlation_matrix"))

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
setClass("UnstructuredCovariance", contains = c("VIRTUAL", "VirtualCovariance"))
##' @rdname VirtualComponentClasses
setClass("DiagonalCovariance", contains = c("VIRTUAL", "VirtualCovariance"))

setValidity("UnstructuredCovariance", function(object) {
    d <- object@dimension
    expected_len <- d * (d + 1) / 2
    if (length(object@cparameters) == expected_len) TRUE else "Incorrect number of cparameters for Unstructured."
})

setValidity("DiagonalCovariance", function(object) {
    if (length(object@cparameters) == 0L) TRUE else "Diagonal models must have 0 cparameters."
})

##' @title Virtual Class for Compound Symmetry (CS) Covariance
##' @description Defines the CS structure, which uses 2 `cparameters` for the
##'   diagonal and off-diagonal correlation elements.
##' @keywords internal
setClass("CSCovariance", contains = c("VIRTUAL", "VirtualCovariance"))
setValidity("CSCovariance", function(object) {
    expected_len <- if (object@dimension > 0) 2L else 0L
    if (length(object@cparameters) == expected_len) TRUE else "CS models must have 2 cparameters (diag and off-diag)."
})

##' @title Virtual Class for Autoregressive Order 1 (AR1) Covariance
##' @description Defines the AR1 structure, which uses `d` `cparameters`,
##'   one for each correlation lag.
##' @keywords internal
setClass("AR1Covariance", contains = c("VIRTUAL", "VirtualCovariance"))
setValidity("AR1Covariance", function(object) {
    if (length(object@cparameters) == object@dimension) TRUE else "AR1 models must have d cparameters (one for each lag)."
})

##' @title Concrete Covariance Structure Classes
##' @description These are the user-facing S4 classes for specifying covariance
##'   structures. They can be initialized with `new("ClassName", dimension = d)`.
##'
##' @name ConcreteCovarianceClasses
##' @export HomogeneousDiagonalCovariance
##' @export HeterogeneousDiagonalCovariance
##' @export HomogeneousCSCovariance
##' @export HeterogeneousCSCovariance
##' @export HomogeneousAR1Covariance
##' @export HeterogeneousAR1Covariance
##' @export HomogeneousUnstructuredCovariance
##' @export HeterogeneousUnstructuredCovariance
NULL

# Create the Concrete Classes
for (variance_type in c("Homogeneous", "Heterogeneous")) {
    for (structure_type in c("Unstructured", "Diagonal", "CS", "AR1")) {
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

      start_params <- get_start_values(temp_obj)
      n_v_params <- if (is(temp_obj, "HomogeneousVariance")) 1L else temp_obj@dimension
      
      args$vparameters <- start_params[seq_len(n_v_params)]
      args$cparameters <- start_params[-seq_len(n_v_params)]
    }
        do.call(callNextMethod, c(list(.Object), args))
})

##' @rdname CovarianceMethods
setMethod("n_parameters", "UnstructuredCovariance", function(object) {
    d <- object@dimension
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    n_c_params <- d * (d + 1) / 2
    n_v_params + n_c_params
})

##' @rdname CovarianceMethods
setMethod("n_parameters", "DiagonalCovariance", function(object) {
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else object@dimension
    n_v_params # No correlation parameters
})

##' @rdname CovarianceMethods
setMethod("n_parameters", "CSCovariance", function(object) {
    d <- object@dimension
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    n_c_params <- if (d > 0) 2L else 0L
    n_v_params + n_c_params
})

##' @rdname CovarianceMethods
setMethod("n_parameters", "AR1Covariance", function(object) {
    d <- object@dimension
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    n_c_params <- d
    n_v_params + n_c_params
})

#' @rdname CovarianceMethods
setMethod("set_parameters", "VirtualCovariance", function(object, value) {
    n_total <- n_parameters(object)
    if (length(value) != n_total) {
        stop(sprintf("Incorrect number of parameters. Expected %d, got %d.", n_total, length(value)))
    }
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else object@dimension
    object@vparameters <- value[seq_len(n_v_params)]
    object@cparameters <- value[-seq_len(n_v_params)]
    validObject(object)
    return(object)
})

##' @rdname InternalCovarianceMethods
setMethod("get_start_values", "VirtualCovariance", function(object) {
    d <- object@dimension
    v_starts <- if (is(object, "HomogeneousVariance")) 0 else rep(0, d)
    c_starts <- if (is(object, "UnstructuredCovariance")) {
        rep(0, d * (d + 1) / 2)
    } else if (is(object, "DiagonalCovariance")) {
        numeric(0)
    } else if (is(object, "CSCovariance")) {
        if (d > 0) c(0, 0.1) else numeric(0)
    } else if (is(object, "AR1Covariance")) {
        rep(0.1, d)
    }
    c(v_starts, c_starts)
})

##' @rdname CovarianceMethods
setMethod("get_parameters", "VirtualCovariance", function(object) c(object@vparameters, object@cparameters))
##' @rdname InternalCovarianceMethods
setMethod("get_lower_bounds", "VirtualCovariance", function(object) rep(-Inf, n_parameters(object)))
##' @rdname InternalCovarianceMethods
setMethod("is_diagonal", "DiagonalCovariance", function(object) TRUE)
##' @rdname InternalCovarianceMethods
setMethod("is_diagonal", "VirtualCovariance", function(object) FALSE)

# SECTION 4: MODEL FITTING INTEGRATION METHODS (GET_LAMBDA, GET_LIND)

##' @rdname InternalCovarianceMethods
setMethod("get_lambda", "UnstructuredCovariance", function(object) {
    d <- object@dimension
    if (d == 0) return(new("dgTMatrix"))
    indices <- get_lower_tri_indices(d)
    x <- as.numeric(indices$i == indices$j)
    sparseMatrix(i = indices$i, j = indices$j, x = x, dims = c(d, d))
})

##' @rdname InternalCovarianceMethods
setMethod("get_lambda", "DiagonalCovariance", function(object) {
    Matrix::Diagonal(n = object@dimension)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lambda", "CSCovariance", function(object) {
    d <- object@dimension
    if (d == 0) return(new("dgTMatrix"))
    indices <- get_lower_tri_indices(d)
    x <- as.numeric(indices$i == indices$j)
    sparseMatrix(i = indices$i, j = indices$j, x = x, dims = c(d, d))
})

##' @rdname InternalCovarianceMethods
setMethod("get_lambda", "AR1Covariance", function(object) {
    d <- object@dimension
    if (d == 0) return(new("dgTMatrix"))
    indices <- get_lower_tri_indices(d)
    x <- as.numeric(indices$i == indices$j)
    sparseMatrix(i = indices$i, j = indices$j, x = x, dims = c(d, d))
})

##' @rdname InternalCovarianceMethods
setMethod("get_lind", "UnstructuredCovariance", function(object) {
    d <- object@dimension
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    n_c_params <- d * (d + 1) / 2
    seq.int(from = n_v_params + 1, length.out = n_c_params)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lind", "DiagonalCovariance", function(object) {
    integer(0)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lind", "CSCovariance", function(object) {
    d <- object@dimension
    if (d == 0) return(integer(0))
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    indices <- get_lower_tri_indices(d)
    ifelse(indices$i == indices$j, n_v_params + 1, n_v_params + 2)
})

##' @rdname InternalCovarianceMethods
setMethod("get_lind", "AR1Covariance", function(object) {
    d <- object@dimension
    if (d == 0) return(integer(0))
    n_v_params <- if (is(object, "HomogeneousVariance")) 1L else d
    indices <- get_lower_tri_indices(d)
    lag <- indices$i - indices$j
    n_v_params + lag + 1
})

# SECTION 5: POST-HOC COMPUTATION & INTERPRETATION METHODS

##' @rdname VirtualComponentClasses
setClass("VirtualStructuredCovariance", contains = c("VIRTUAL", "VirtualCovariance"))
setIs("DiagonalCovariance", "VirtualStructuredCovariance")
setIs("CSCovariance", "VirtualStructuredCovariance")
setIs("AR1Covariance", "VirtualStructuredCovariance")

##' @rdname InternalCovarianceMethods
setMethod("build_correlation_matrix", "DiagonalCovariance", function(object) Diagonal(object@dimension))
##' @rdname InternalCovarianceMethods
setMethod("build_correlation_matrix", "CSCovariance", function(object) {
    d <- object@dimension
    if (d < 2) return(Diagonal(d))
    rho <- tanh(object@cparameters[2])
    R <- Matrix(rho, nrow = d, ncol = d); diag(R) <- 1
    return(R)
})
##' @rdname InternalCovarianceMethods
setMethod("build_correlation_matrix", "AR1Covariance", function(object) {
    d <- object@dimension
    if (d < 2) return(Diagonal(d))
    rho <- tanh(object@cparameters[2])
    time_diffs <- abs(outer(1:d, 1:d, "-"))
    rho^time_diffs
})
##' @rdname InternalCovarianceMethods
setMethod("build_correlation_matrix", "UnstructuredCovariance", function(object) {
    Sigma <- compute_covariance_matrix(object)
    cov2cor(Sigma)
})

##' @rdname CovarianceMethods
setMethod("compute_covariance_matrix", "VirtualStructuredCovariance", function(object, data_context = NULL) {
    if (object@dimension == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
    R <- build_correlation_matrix(object)
    if (is(object, "HomogeneousVariance")) {
        sigma_sq <- exp(object@vparameters[1])
        return(sigma_sq * R)
    } else {
        st_devs <- exp(0.5 * object@vparameters)
        D <- Diagonal(object@dimension, x = st_devs)
        return(D %*% R %*% D)
    }
})

##' @rdname CovarianceMethods
setMethod("compute_covariance_matrix", "UnstructuredCovariance", function(object, data_context = NULL) {
    L <- get_chol_from_params(object@cparameters, object@dimension)
    tcrossprod(L)
})

##' @rdname CovarianceMethods
setMethod("get_cholesky_factor", "VirtualCovariance", function(object, data_context = NULL) {
    chol(compute_covariance_matrix(object, data_context))
})
##' @rdname CovarianceMethods
setMethod("get_cholesky_factor", "UnstructuredCovariance", function(object, data_context = NULL) {
    get_chol_from_params(object@cparameters, object@dimension)
})

##' @rdname CovarianceMethods
setMethod("compute_log_det_covariance_matrix", "UnstructuredCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    if (d == 0) return(0)
    diag_indices <- vech_diag_indices(d)
    2 * sum(object@cparameters[diag_indices])
})
##' @rdname CovarianceMethods
setMethod("compute_log_det_covariance_matrix", "VirtualStructuredCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    if (d == 0) return(0)
    R <- build_correlation_matrix(object)
    log_det_R <- as.numeric(determinant(R, logarithm = TRUE)$modulus)
    log_det_V <- if (is(object, "HomogeneousVariance")) {
        d * object@vparameters[1]
    } else {
        sum(object@vparameters)
    }
    return(log_det_V + log_det_R)
})

##' @rdname CovarianceMethods
setMethod("compute_inverse_covariance_matrix", "UnstructuredCovariance", function(object, data_context = NULL) {
    L <- get_cholesky_factor(object)
    inv_L <- solve(L)
    crossprod(inv_L)
})
##' @rdname CovarianceMethods
setMethod("compute_inverse_covariance_matrix", "VirtualStructuredCovariance", function(object, data_context = NULL) {
    d <- object@dimension
    if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
    inv_R <- solve(build_correlation_matrix(object))
    if (is(object, "HomogeneousVariance")) {
        inv_sigma_sq <- exp(-object@vparameters[1])
        return(inv_sigma_sq * inv_R)
    } else {
        inv_st_devs <- exp(-0.5 * object@vparameters)
        inv_D <- Diagonal(d, x = inv_st_devs)
        return(inv_D %*% inv_R %*% inv_D)
    }
})

##' @rdname CovarianceMethods
setMethod("get_interpretable_parameters", "VirtualStructuredCovariance", function(object) {
    params <- list()
    if (is(object, "HomogeneousVariance")) {
        params$st_dev <- exp(0.5 * object@vparameters[1])
    } else {
        params$st_devs <- exp(0.5 * object@vparameters)
    }
    if (is(object, "CSCovariance") && object@dimension > 1) {
        params$correlation <- tanh(object@cparameters[2])
    } else if (is(object, "AR1Covariance") && object@dimension > 1) {
        params$correlation_lag1 <- tanh(object@cparameters[2])
    }
    return(params)
})

##' @rdname CovarianceMethods
setMethod("get_interpretable_parameters", "UnstructuredCovariance", function(object) {
    Sigma <- compute_covariance_matrix(object)
    st_devs <- sqrt(diag(Sigma))
    
    params <- list()
    if (is(object, "HomogeneousVariance")) {
        params$st_dev <- mean(st_devs) 
    } else {
        params$st_devs <- st_devs
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
            } else if (!is.null(params$correlation_lag1)) {
                cat(sprintf("    Lag-1 Correlation (rho): %.4f\n", params$correlation_lag1))
            }
        }
    }, silent = TRUE)
})


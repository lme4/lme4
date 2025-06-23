# S4 Class and Method Definition for Covariance Structures

library(Matrix)

# SECTION 1: TOP-LEVEL VIRTUAL CLASS AND generics

##' Virtual Class for All Covariance Structures
##'
##' This abstract S4 class defines the common interface for all covariance models.
##' It serves as a base class for inheritance.
##'
##' @slot dimension An integer indicating the size of the square covariance matrix.
##' @export
setClass("VirtualCovariance",
	contains = "VIRTUAL",
	slots = c(dimension = "integer")
)

# Define all generic functions 
##' @rdname n_parameters
setGeneric("n_parameters", function(object) standardGeneric("n_parameters"))
##' @rdname get_parameters
setGeneric("get_parameters", function(object) standardGeneric("get_parameters"))
##' @rdname set_parameters
setGeneric("set_parameters", function(object, value) standardGeneric("set_parameters"))
##' @rdname get_start_values
setGeneric("get_start_values", function(object) standardGeneric("get_start_values"))
##' @rdname get_lower_bounds
setGeneric("get_lower_bounds", function(object) standardGeneric("get_lower_bounds"))
##' @rdname get_interpretable_parameters
setGeneric("get_interpretable_parameters", function(object) standardGeneric("get_interpretable_parameters"))
##' @rdname is_diagonal
setGeneric("is_diagonal", function(object) standardGeneric("is_diagonal"))
##' @rdname compute_covariance_matrix
setGeneric("compute_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_covariance_matrix"))
##' @rdname compute_inverse_covariance_matrix
setGeneric("compute_inverse_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_inverse_covariance_matrix"))
##' @rdname compute_log_det_covariance_matrix
setGeneric("compute_log_det_covariance_matrix", function(object, data_context = NULL) standardGeneric("compute_log_det_covariance_matrix"))
##' @rdname get_cholesky_factor
setGeneric("get_cholesky_factor", function(object, data_context = NULL) standardGeneric("get_cholesky_factor"))

##' @rdname get_coordinates
setGeneric("get_coordinates", function(data_context) standardGeneric("get_coordinates"))

# SECTION 2: COMPONENT CLASS HIERARCHIES

## Variance Model Components

##' Virtual Class for Variance Models
##'
##' This abstract class is the parent for all variance model components.
##' @slot parameters A numeric vector storing the unconstrained variance parameters.
##' @export
setClass("VirtualVarianceModel", contains = "VIRTUAL",
	slots = c(parameters = "numeric")
)

##' Homogeneous Variance Model
##'
##' Represents a single, common variance for all dimensions.
##' @export
setClass("HomogeneousVarianceModel", contains = "VirtualVarianceModel")
setValidity("HomogeneousVarianceModel", function(object) {
	if (length(object@parameters) != 1) "Expected 1 parameter for HomogeneousVarianceModel." else TRUE
})

##' Heterogeneous Variance Model
##'
##' Represents a unique variance for each dimension.
##' @slot dimension An integer specifying the number of unique variances.
##' @export
setClass("HeterogeneousVarianceModel",
	contains = "VirtualVarianceModel",
	slots = c(dimension = "integer")
)
setValidity("HeterogeneousVarianceModel", function(object) {
	if (length(object@parameters) != object@dimension) "Number of parameters must equal dimension." else TRUE
})

## Correlation Model Components

##' Virtual Class for Correlation Models
##'
##' This abstract class is the parent for all correlation model components.
##' @slot parameters A numeric vector storing the unconstrained correlation parameters.
##' @export
setClass("VirtualCorrelationModel", contains = "VIRTUAL",
	slots = c(parameters = "numeric")
)

##' Identity Correlation Model
##'
##' Represents an Identity (no correlation) structure. Has 0 parameters.
##' @export
setClass("IdentityCorrelationModel", contains = "VirtualCorrelationModel")
setValidity("IdentityCorrelationModel", function(object) {
	if (length(object@parameters) != 0) "Expected 0 parameters for IdentityCorrelationModel." else TRUE
})

##' Compound Symmetry (CS) Correlation Model
##'
##' Represents a Compound Symmetry (CS) correlation structure.
##' @export
setClass("CSCorrelationModel", contains = "VirtualCorrelationModel")
setValidity("CSCorrelationModel", function(object) {
	if (length(object@parameters) != 1) return("Expected 1 parameter for CSCorrelationModel.")
	TRUE
})

##' Autoregressive Order 1 (AR1) Correlation Model
##'
##' Represents a first-order Autoregressive (AR1) correlation structure.
##' @export
setClass("AR1CorrelationModel", contains = "VirtualCorrelationModel")
setValidity("AR1CorrelationModel", function(object) {
	if (length(object@parameters) != 1) return("Expected 1 parameter for AR1CorrelationModel.")
	rho <- tanh(object@parameters[1])
	if (abs(rho) >= 1) return("Correlation rho must be in (-1, 1).")
	TRUE
})

##' Virtual Class for Kernel-Based Correlation Models
##' @export
setClass("VirtualKernelModel", contains = "VirtualCorrelationModel")

##' Squared Exponential (Gaussian) Kernel Model
##' @export
setClass("SquaredExpKernelModel", contains = "VirtualKernelModel")
setValidity("SquaredExpKernelModel", function(object) {
	if (length(object@parameters) != 1) "Expected 1 parameter (log-length-scale)." else TRUE
})

##' Exponential Kernel Model
##' @export
setClass("ExponentialKernelModel", contains = "VirtualKernelModel")
setValidity("ExponentialKernelModel", function(object) {
	if (length(object@parameters) != 1) "Expected 1 parameter (log-length-scale)." else TRUE
})

##' Matern Kernel Model
##' @export
setClass("MaternKernelModel", contains = "VirtualKernelModel") # Params: log(length_scale), log(nu)
setValidity("MaternKernelModel", function(object) {
	if (length(object@parameters) != 2) "Expected 2 parameters (log-length-scale, log-nu)." else TRUE
})

##' Spherical Kernel Model
##' @export
setClass("SphericalKernelModel", contains = "VirtualKernelModel")
setValidity("SphericalKernelModel", function(object) {
	if (length(object@parameters) != 1) "Expected 1 parameter (log-range)." else TRUE
})

##' Rational Quadratic Kernel Model
##' @export
setClass("RationalQuadraticKernelModel", contains = "VirtualKernelModel") # Params: log(alpha), log(length_scale)
setValidity("RationalQuadraticKernelModel", function(object) {
	if (length(object@parameters) != 2) "Expected 2 parameters (log-alpha, log-length-scale)." else TRUE
})

# TOP-LEVEL CONCRETE CLASSES

##' Flexible D-R-D Covariance Structure
##'
##' A flexible covariance structure for the D*R*D family of models, where
##' Sigma = D %*% R %*% D. D is a diagonal matrix of standard deviations,
##' and R is a correlation matrix.
##'
##' @slot dimension An integer indicating the size of the square covariance matrix.
##' @slot variance_model An object from a `VirtualVarianceModel` class.
##' @slot correlation_model An object from a `VirtualCorrelationModel` class.
##' @export
setClass("FlexibleCovariance",
	contains = "VirtualCovariance",
	slots = c(
		dimension = "integer",
		variance_model = "VirtualVarianceModel",
		correlation_model = "VirtualCorrelationModel"
	)
)

setValidity("FlexibleCovariance", function(object) {
	validObject(object@variance_model)
	validObject(object@correlation_model)
	d <- object@dimension
	if (is(object@correlation_model, "CSCorrelationModel") && d > 1) {
		rho <- tanh(object@correlation_model@parameters[1])
		if (rho <= -1 / (d - 1)) {
			return(sprintf("Correlation rho of %f is not valid for dimension %d.", rho, d))
		}
	}
	TRUE
})

##' Unstructured Covariance Structure
##'
##' Represents a general, unstructured positive semi-definite covariance matrix,
##' parameterized using a log-Cholesky decomposition.
##'
##' @slot parameters A numeric vector storing the unconstrained parameters
##' (the vectorized log-Cholesky factor).
##' @export
setClass("UnstructuredCovariance",
	contains = "VirtualCovariance",
	slots = c(parameters = "numeric")
)

setValidity("UnstructuredCovariance", function(object) {
	expected_len <- (object@dimension * (object@dimension + 1)) / 2
	if (length(object@parameters) != expected_len) {
		return(sprintf("Expected %d parameters, but slot has %d.", expected_len, length(object@parameters)))
	}
	TRUE
})

# SECTION 4: METHODS FOR THE D*R*D FAMILY (FLEXIBLE COVARIANCE)

##' Compute Correlation Matrix
##'
##' Generic function to compute the correlation matrix (R) from a correlation model component.
##' @param object An S4 object inheriting from `VirtualCorrelationModel`.
##' @param d The dimension of the matrix.
##' @param data_context Optional. Auxiliary data needed for computation.
##' @return A matrix representing the correlation matrix.
##' @export
setGeneric("compute_correlation_matrix", function(object, d, data_context = NULL) standardGeneric("compute_correlation_matrix"))

##' Compute Log-Determinant of Correlation Matrix
##'
##' Generic function to compute the log-determinant of the correlation matrix.
##' @param object An S4 object inheriting from `VirtualCorrelationModel`.
##' @param d The dimension of the matrix.
##' @param data_context Optional. Auxiliary data needed for computation.
##' @return A numeric value representing the log-determinant.
##' @export
setGeneric("compute_log_det_correlation_matrix", function(object, d, data_context = NULL) standardGeneric("compute_log_det_correlation_matrix"))

##' Compute Inverse of Correlation Matrix
##'
##' Generic function to compute the inverse of the correlation matrix (R_inverse).
##' @param object An S4 object inheriting from `VirtualCorrelationModel`.
##' @param d The dimension of the matrix.
##' @param data_context Optional. Auxiliary data needed for computation.
##' @return A matrix representing the inverse correlation matrix.
##' @export
setGeneric("compute_inverse_correlation_matrix", function(object, d, data_context = NULL) standardGeneric("compute_inverse_correlation_matrix"))

##' @rdname compute_correlation_matrix
setMethod("compute_correlation_matrix", "IdentityCorrelationModel", function(object, d, data_context) Diagonal(d))
##' @rdname compute_log_det_correlation_matrix
setMethod("compute_log_det_correlation_matrix", "IdentityCorrelationModel", function(object, d, data_context) 0)
##' @rdname compute_inverse_correlation_matrix
setMethod("compute_inverse_correlation_matrix", "IdentityCorrelationModel", function(object, d, data_context) Diagonal(d))
##' @rdname compute_correlation_matrix
setMethod("compute_correlation_matrix", "CSCorrelationModel", function(object, d, data_context) {
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	rho <- tanh(object@parameters[1])
	R <- Matrix(rho, nrow = d, ncol = d)
    diag(R) <- 1
	return(R)
})
##' @rdname compute_log_det_correlation_matrix
setMethod("compute_log_det_correlation_matrix", "CSCorrelationModel", function(object, d, data_context) {
	if (d == 0) return(0)
	rho <- tanh(object@parameters[1])
	(d - 1) * log(1 - rho) + log(1 + (d - 1) * rho)
})
##' @rdname compute_inverse_correlation_matrix
setMethod("compute_inverse_correlation_matrix", "CSCorrelationModel", function(object, d, data_context) {
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	rho <- tanh(object@parameters[1])
	common_divisor <- (1 - rho) * (1 + (d - 1) * rho)
	val_A <- (1 + (d - 2) * rho) / common_divisor
	val_B <- -rho / common_divisor
	Inv_R <- Matrix(val_B, nrow = d, ncol = d)
    diag(Inv_R) <- val_A
	return(Inv_R)
})

##' @rdname compute_correlation_matrix
setMethod("compute_correlation_matrix", "AR1CorrelationModel", function(object, d, data_context) {
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	rho <- tanh(object@parameters[1])
	time_diffs <- abs(outer(1:d, 1:d, "-"))
    rho^time_diffs
})
##' @rdname compute_log_det_correlation_matrix
setMethod("compute_log_det_correlation_matrix", "AR1CorrelationModel", function(object, d, data_context) {
	if (d == 0) return(0)
	rho <- tanh(object@parameters[1])
	(d - 1) * log(1 - rho^2)
})
##' @rdname compute_inverse_correlation_matrix
setMethod("compute_inverse_correlation_matrix", "AR1CorrelationModel", function(object, d, data_context) {
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	if (d == 1) return(Matrix(1, 1, 1))

	rho <- tanh(object@parameters[1])
	scaling_factor <- 1 / (1 - rho^2)

	diag_vals <- rep(1 + rho^2, d)
	diag_vals[1] <- 1
	diag_vals[d] <- 1

	Inv_R <- bandSparse(d, d, k = c(-1, 0, 1),
		diagonals = list(rep(-rho, d - 1),
			diag_vals,
			rep(-rho, d - 1)))

	return(scaling_factor * Inv_R)
})

## KERNEL METHODS
##' @rdname compute_correlation_matrix
setMethod("compute_correlation_matrix", "SquaredExpKernelModel", function(object, d, data_context) {
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	coords <- get_coordinates(data_context)
	if (is.null(coords)) stop("Coordinates must be provided in data_context for kernel models.")
	dist_sq <- as.matrix(dist(coords)^2)
	length_scale <- exp(object@parameters[1])
	R <- exp(-dist_sq / (2 * length_scale^2))
	return(R)
})

##' @rdname compute_correlation_matrix
setMethod("compute_correlation_matrix", "ExponentialKernelModel", function(object, d, data_context) {
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	coords <- get_coordinates(data_context)
	if (is.null(coords)) stop("Coor<LeftMouse>dinates must be provided in data_context for kernel models.")
	dists <- as.matrix(dist(coords))
	length_scale <- exp(object@parameters[1])
	R <- exp(-dists / length_scale)
	return(R)
})

##' @rdname compute_correlation_matrix
setMethod("compute_correlation_matrix", "MaternKernelModel", function(object, d, data_context) {
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	coords <- get_coordinates(data_context)
	if (is.null(coords)) stop("Coordinates must be provided in data_context for kernel models.")
	dists <- as.matrix(dist(coords))
	length_scale <- exp(object@parameters[1])
	nu <- exp(object@parameters[2])
	dists[dists == 0] <- 1e-9
	term1 <- (sqrt(2 * nu) * dists) / length_scale
	term2 <- besselK(term1, nu)
	R <- (2^(1 - nu) / gamma(nu)) * (term1^nu) * term2
	diag(R) <- 1
	return(R)
})

##' @rdname compute_correlation_matrix
setMethod("compute_correlation_matrix", "SphericalKernelModel", function(object, d, data_context) {
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	coords <- get_coordinates(data_context)
	if (is.null(coords)) stop("Coordinates must be provided in data_context for kernel models.")
	dists <- as.matrix(dist(coords))
	range_param <- exp(object@parameters[1])
	R <- 1 - 1.5 * (dists / range_param) + 0.5 * (dists / range_param)^3
	R[dists > range_param] <- 0
	return(R)
})

##' @rdname compute_correlation_matrix
setMethod("compute_correlation_matrix", "RationalQuadraticKernelModel", function(object, d, data_context) {
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	coords <- get_coordinates(data_context)
	if (is.null(coords)) stop("Coordinates must be provided in data_context for kernel models.")
	dist_sq <- as.matrix(dist(coords)^2)
	alpha <- exp(object@parameters[1])
	length_scale <- exp(object@parameters[2])
	R <- (1 + dist_sq / (2 * alpha * length_scale^2))^(-alpha)
	return(R)
})

## Methods for the Composer Class (`FlexibleCovariance`)
##' @rdname compute_covariance_matrix
setMethod("compute_covariance_matrix", "FlexibleCovariance", function(object, data_context = NULL) {
	if (object@dimension == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	R <- compute_correlation_matrix(object@correlation_model, d = object@dimension, data_context = data_context)
	if (is(object@variance_model, "HomogeneousVarianceModel")) {
		sigma_sq <- exp(object@variance_model@parameters[1])
		return(sigma_sq * R)
	} else {
		st_devs <- exp(0.5 * object@variance_model@parameters)
		D <- Diagonal(object@dimension, x = st_devs)
		return(D %*% R %*% D)
	}
})

##' @rdname compute_log_det_covariance_matrix
setMethod("compute_log_det_covariance_matrix", "FlexibleCovariance", function(object, data_context = NULL) {
	d <- object@dimension
	if (d == 0) return(0)
	log_det_R <- compute_log_det_correlation_matrix(object@correlation_model, d, data_context)
	if (is(object@variance_model, "HomogeneousVarianceModel")) {
		log_det_V <- d * object@variance_model@parameters[1]
	} else {
		log_det_V <- sum(object@variance_model@parameters)
	}
	return(log_det_V + log_det_R)
})

##' @rdname compute_inverse_covariance_matrix
setMethod("compute_inverse_covariance_matrix", "FlexibleCovariance", function(object, data_context = NULL) {
	d <- object@dimension
	if (d == 0) return(new("dsyMatrix", Dim = c(0L, 0L)))
	inv_R <- compute_inverse_correlation_matrix(object@correlation_model, d, data_context)
	if (is(object@variance_model, "HomogeneousVarianceModel")) {
		inv_sigma_sq <- exp(-object@variance_model@parameters[1])
		return(inv_sigma_sq * inv_R)
	} else {
		inv_st_devs <- exp(-0.5 * object@variance_model@parameters)
		inv_D <- Diagonal(d, x = inv_st_devs)
		return(inv_D %*% inv_R %*% inv_D)
	}
})

##' @rdname get_cholesky_factor
setMethod("get_cholesky_factor", "FlexibleCovariance", function(object, data_context = NULL) {
	chol(compute_covariance_matrix(object, data_context))
})

##' @rdname n_parameters
setMethod("n_parameters", "FlexibleCovariance", function(object) {
	length(object@variance_model@parameters) + length(object@correlation_model@parameters)
})

##' @rdname get_parameters
setMethod("get_parameters", "FlexibleCovariance", function(object) {
	c(object@variance_model@parameters, object@correlation_model@parameters)
})

##' @rdname set_parameters
setMethod("set_parameters", "FlexibleCovariance", function(object, value) {
	n_var_params <- length(object@variance_model@parameters)
	n_cor_params <- length(object@correlation_model@parameters)
	if (length(value) != n_var_params + n_cor_params) {
		stop("Incorrect number of parameters provided.")
	}
	object@variance_model@parameters <- value[1:n_var_params]
	if (n_cor_params > 0) {
		object@correlation_model@parameters <- value[(n_var_params + 1):(n_var_params + n_cor_params)]
	}
	validObject(object)
	return(object)
})

##' @rdname get_start_values
setMethod("get_start_values", "FlexibleCovariance", function(object) {
	d <- object@dimension
	var_starts <- if (is(object@variance_model, "HomogeneousVarianceModel")) 0 else rep(0, d)
	cor_model <- object@correlation_model
	cor_starts <- if (is(cor_model, "IdentityCorrelationModel")) {
		numeric(0)
	} else if (is(cor_model, "MaternKernelModel") || is(cor_model, "RationalQuadraticKernelModel")) {
		c(0.1, 0.1)
	} else {
		0.1
	}
	c(var_starts, cor_starts)
})

##' @rdname get_lower_bounds
setMethod("get_lower_bounds", "FlexibleCovariance", function(object) {
	rep(-Inf, n_parameters(object))
})

##' @rdname is_diagonal
setMethod("is_diagonal", "FlexibleCovariance", function(object) {
	is(object@correlation_model, "IdentityCorrelationModel") || all(object@correlation_model@parameters == 0)
})

##' @rdname get_interpretable_parameters
setMethod("get_interpretable_parameters", "FlexibleCovariance", function(object) {
	params <- list()
	if (is(object@variance_model, "HomogeneousVarianceModel")) {
		params$st_dev <- exp(0.5 * object@variance_model@parameters[1])
	} else {
		params$st_devs <- exp(0.5 * object@variance_model@parameters)
	}
	cor_model <- object@correlation_model
	if (is(cor_model, "SquaredExpKernelModel") || is(cor_model, "ExponentialKernelModel") ||
		is(cor_model, "SphericalKernelModel")) {
		params$length_scale <- exp(cor_model@parameters[1])
	} else if (is(cor_model, "MaternKernelModel")) {
		params$length_scale <- exp(cor_model@parameters[1])
		params$nu <- exp(cor_model@parameters[2])
	} else if (is(cor_model, "RationalQuadraticKernelModel")) {
		params$alpha <- exp(cor_model@parameters[1])
		params$length_scale <- exp(cor_model@parameters[2])
	} else if (!is(cor_model, "IdentityCorrelationModel")) {
		params$correlation <- tanh(cor_model@parameters[1])
	}
	return(params)
})

##' Show Method for FlexibleCovariance
##'
##' Prints a summary of the FlexibleCovariance Object.
##' @param object A `FlexibleCovariance` object.
##' @export
setMethod("show", "FlexibleCovariance", function(object) {
	var_name <- sub("VarianceModel", "", class(object@variance_model))
	cor_name <- sub("KernelModel", "", class(object@correlation_model))
	cor_name <- sub("CorrelationModel", "", cor_name)
	model_name <- sprintf("%s-%s", var_name, cor_name)

	cat(sprintf("%s Covariance object (dimension: %d)\n", model_name, object@dimension))
	try({
		params <- get_interpretable_parameters(object)
		if (!is.null(params$st_devs)) {
			st_devs_formatted <- paste(sprintf("%.4f", params$st_devs), collapse = " ")
			cat("  Standard Deviations:", st_devs_formatted, "\n")
		} else if (!is.null(params$st_dev)) {
			cat(sprintf("  Standard Deviation (sigma): %.4f\n", params$st_dev))
		}
		if (!is.null(params$correlation)) {
			cat(sprintf("  Correlation (rho): %.4f\n", params$correlation))
		}
		if (!is.null(params$length_scale)) {
			cat(sprintf("  length_scale: %.4f\n", params$length_scale))
		}
		if (!is.null(params$nu)) {
			cat(sprintf("  nu: %.4f\n", params$nu))
		}
		if (!is.null(params$alpha)) {
			cat(sprintf("  alpha: %.4f\n", params$alpha))
		}
		if (!is.null(params$c_intercept)) {
			cat(sprintf("  c_intercept: %.4f\n", params$c_intercept))
		}
	}, silent = TRUE)
})

# SECTION 5: METHODS FOR THE LOG-CHOLESKY FAMILY (UNSTRUCTURED COVARIANCE)

##' @rdname n_parameters
setMethod("n_parameters", "UnstructuredCovariance", function(object) {
	d <- object@dimension
	d * (d + 1) / 2
})

##' @rdname get_parameters
setMethod("get_parameters", "UnstructuredCovariance", function(object) {
	object@parameters
})

##' @rdname set_parameters
setMethod("set_parameters", "UnstructuredCovariance", function(object, value) {
	if (length(value) != n_parameters(object)) {
		stop("Incorrect number of parameters for UnstructuredCovariance.")
	}
	object@parameters <- value
	validObject(object)
	return(object)
})

##' @rdname get_start_values
setMethod("get_start_values", "UnstructuredCovariance", function(object) {
	rep(0, n_parameters(object))
})

##' @rdname get_lower_bounds
setMethod("get_lower_bounds", "UnstructuredCovariance", function(object) {
	rep(-Inf, n_parameters(object))
})

##' @rdname get_interpretable_parameters
setMethod("get_interpretable_parameters", "UnstructuredCovariance", function(object) {
	Sigma <- compute_covariance_matrix(object)
	list(st_devs = sqrt(diag(Sigma)))
})

##' @rdname is_diagonal
setMethod("is_diagonal", "UnstructuredCovariance", function(object) {
	FALSE
})

##' @rdname get_cholesky_factor
setMethod("get_cholesky_factor", "UnstructuredCovariance", function(object, data_context = NULL) {
	get_chol_from_params(object@parameters, object@dimension)
})

##' @rdname compute_covariance_matrix
setMethod("compute_covariance_matrix", "UnstructuredCovariance", function(object, data_context = NULL) {
	L <- get_cholesky_factor(object)
	tcrossprod(L)
})

##' @rdname compute_log_det_covariance_matrix
setMethod("compute_log_det_covariance_matrix", "UnstructuredCovariance", function(object, data_context = NULL) {
	diag_indices <- vech_diag_indices(object@dimension)
	2 * sum(object@parameters[diag_indices])
})

##' @rdname compute_inverse_covariance_matrix
setMethod("compute_inverse_covariance_matrix", "UnstructuredCovariance", function(object, data_context = NULL) {
	L <- get_cholesky_factor(object)
	inv_L <- solve(L)
	crossprod(inv_L)
})

##' Show Method for UnstructuredCovariance
##'
##' Prints a summary of the UnstructuredCovariance object.
##' @param object A `UnstructuredCovariance` object.
##' @export
setMethod("show", "UnstructuredCovariance", function(object) {
	cat(sprintf("UnstructuredCovariance object (dimension: %d)\n", object@dimension))
	try({
		params <- get_interpretable_parameters(object)
		st_devs_formatted <- paste(sprintf("%.4f", params$st_devs), collapse = " ")
		cat("  Standard Deviations:", st_devs_formatted, "\n")
	}, silent = TRUE)
})

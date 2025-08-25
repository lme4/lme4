#' S4 Object Factory for Covariance Structures
#'
#' An internal helper function that translates the parsed components of a single
#' random effect term into the corresponding S4 covariance object.
#'
#' @param type A character string for the covariance structure type (e.g., "ar1", "us").
#' @param add_args_call The language object of the wrapper function call (e.g., quote(ar1(d = 2))).
#' @param cnms A character vector of column names for a single random effect term; its length determines the dimension of the covariance structure.
#'
#' @return An S4 object inheriting from `VirtualCovariance`.
#' @keywords internal
create_covariance_object_from_term <- function(type, cnms, add_args_call) {

    # Define the lookup table 
    cov_struct_table <- list(
        cs = list(prefix = "CS"), 
        ar1 = list(prefix = "AR1"), 
        diag = list(prefix = "Diagonal"), 
        us = list(prefix = "Unstructured")
    )  
    
    spec <- cov_struct_table[[type]]
    if (is.null(spec)) {
        stop("Unknown covariance structure provided: ", type)
    }
    structure_prefix <- spec$prefix 

    args <- as.list(add_args_call)[-1]

    if (type == "us") {
        return(new("UnstructuredCovariance", dimension = length(cnms)))
    }

    is_hom <- args$hom %||% FALSE
    variance_prefix <- if(is_hom) "Homogeneous" else "Heterogeneous"

    full_class_name <- paste0(variance_prefix, structure_prefix, "Covariance")

    dimension <- length(cnms)


    return(new(full_class_name, dimension = dimension))
}

#' Parse Covariance Structures from a Model Formula
#'
#' An internal function that uses `reformulas::splitForm` to parse a model
#' formula containing special covariance structures (e.g., `ar1(...)`) and
#' generates a list of corresponding S4 covariance objects.
#'
#' @param formula The model formula.
#' @param data The model data frame.
#'
#' @return A list of S4 objects, where each element corresponds to a
#'   random effect term in the formula. Returns an empty list if no
#'   random effects are found.
#' @keywords internal
parse_model_formula <- function(formula, data) {
    specials_list <- c("ar1", "cs", "diag", "us") 
    split_formula <- reformulas::splitForm(formula, specials = specials_list)
    s4_object_list <- list() 
    
    if (!is.null(split_formula$reTrmFormulas)) {
        temp_reTrms <- reformulas::mkReTrms(split_formula$reTrmFormulas, data, calc.lambdat = FALSE)
        
        for (i in seq_along(split_formula$reTrmClasses)) {
            type <- split_formula$reTrmClasses[i]
            add_args_call <- split_formula$reTrmAddArgs[[i]]
            cnms <- temp_reTrms$cnms[[i]]
            
            s4_object_list[[i]] <- create_covariance_object_from_term(
                type = type,
                cnms = cnms,
                add_args_call = add_args_call
            )
        }
    }
    return(s4_object_list)
}

#' Build Covariance Parameter Structures
#'
#' This function builds the Lambdat, Lind, theta, and lower components for a
#' mixed model. It uses S4 dispatch on a list of covariance objects to
#' determine the parameter structure for each random-effect term. For each
#' term, it constructs a block-diagonal `Lambdat` matrix using `kronecker()`
#' and then combines these blocks using `Matrix::bdiag()`.
#'
#' @param reTrms A list object returned by `reformulas::mkReTrms` (with
#'   calc.lambdat=FALSE). It must contain `$cnms` and `$nl.
#' @param s4_object_list A list of S4 covariance objects, one for each RE term.
#' @return A list with components: `Lambdat`, `Lind`, `theta`, and `lower`.
#' @keywords internal
mkReLambdat <- function(reTrms, s4_object_list) {
    # First pass: expand parameters and calculate actual parameter counts
    theta_list <- list()
    lind_list <- list()
    lower_list <- list()
    actual_nth_per_term <- integer(length(s4_object_list))
    
    for (i in seq_along(s4_object_list)) {
        s4_obj <- s4_object_list[[i]]
        original_theta <- get_start_values(s4_obj)
        
        if (needs_parameter_expansion(s4_obj)) {
            expansion <- expand_parameters_for_optimization(s4_obj, original_theta)
            theta_list[[i]] <- expansion$expanded_theta
            expanded_lind <- expansion$expanded_lind
            # Extend lower bounds to match expanded parameters
            original_lower <- get_lower_bounds(s4_obj)
            lower_list[[i]] <- c(original_lower, rep(original_lower[2], length(expansion$expanded_theta) - 2))
        } else {
            theta_list[[i]] <- original_theta
            expanded_lind <- get_lind(s4_obj)
            lower_list[[i]] <- get_lower_bounds(s4_obj)
        }
        
        actual_nth_per_term[i] <- length(theta_list[[i]])  # Use actual count after expansion
    }
    
    # Recalculate thoff with actual parameter counts
    thoff <- cumsum(c(0L, actual_nth_per_term))
    
    # Second pass: build Lind with correct offsets and Lambdat blocks
    lambdat_blocks <- list()
    for (i in seq_along(s4_object_list)) {
        s4_obj <- s4_object_list[[i]]
        nl <- reTrms$nl[i]
        L_template <- get_lambda(s4_obj)
        
        # Build Lambdat block
        lambdat_blocks[[i]] <- kronecker(Diagonal(nl), L_template)
        
        # Get the expanded Lind from first pass
        if (needs_parameter_expansion(s4_obj)) {
            expansion <- expand_parameters_for_optimization(s4_obj, get_start_values(s4_obj))
            expanded_lind <- expansion$expanded_lind
        } else {
            expanded_lind <- get_lind(s4_obj)
        }
        
        if (length(expanded_lind) > 0) {
            lind_list[[i]] <- rep(expanded_lind, times = nl) + thoff[i]
        } else {
            lind_list[[i]] <- integer(0)
        }
    }
    
    # Combine results
    theta <- unlist(theta_list)
    lower <- unlist(lower_list)
    Lambdat <- do.call(Matrix::bdiag, lambdat_blocks)
    Lind <- unlist(lind_list)
    
    if (length(Lambdat@x) != length(Lind)) {
        stop("Internal error: Mismatch between Lambdat@x length (", 
             length(Lambdat@x), ") and Lind length (", length(Lind), ")")
    }
    
    Lambdat@x[] <- theta[Lind]
    
    list(
        Lambdat = Lambdat,
        Lind = as.integer(Lind),
        theta = theta,
        lower = lower
    )
}

#' Store Covariance Structure Information in reTrms
#' 
#' This function attaches S4 covariance structure object to reTrms as 
#' attribute, following glmmTMB pattern of storing structure metadata. 
#'
#' @param reTrms The reTrms object from mkReTrms
#' @param s4_object_list  of S4 covariance objects from parse_model_formula
#' @return Modifed reTrms object with structure information attached 
#' @keywords internal
store_structure_info <- function(reTrms, s4_object_list) {

     attr(reTrms, "cov_structures") <- s4_object_list

    structure_types <- sapply(s4_object_list, function(obj) {
        if (is(obj, "UnstructuredCovariance")) "us"
        else if (is(obj, "CSCovariance")) "cs" 
        else if (is(obj, "AR1Covariance")) "ar1"
        else if (is(obj, "DiagonalCovariance")) "diag"
        else "unknown"
    })

    attr(reTrms, "structure_types") <- structure_types

    param_sizes <- sapply(s4_object_list, function(obj) {
        n_parameters(obj)
    })
    attr(reTrms, "param_sizes") <- param_sizes

    return(reTrms)
}

#' Extract Structure Information from merMod Object
#'
#' Helper function to retrieve stored structure information from a fitted merMod object.
#'
#' @param merMod A fitted merMod object
#' @return List containing structure objects, types, and parameter information
#' @keywords internal
extract_structure_info <- function(merMod) {
    
    cov_structures <- attr(merMod, "cov_structures")
    structure_types <- attr(merMod, "structure_types") 
    param_sizes <- attr(merMod, "param_sizes")
    
    if (is.null(cov_structures)) {
         return(NULL)
    }
    
    return(list(
        structures = cov_structures,
        types = structure_types,
        param_sizes = param_sizes,
        n_terms = length(cov_structures)
    ))
}


##' Integrate Structure Information into lFormula
##'
##' This function should be called from lFormula to store structure information
##' when covariance structures are detected in the formula.
##'
##' @param formula Model formula
##' @param data Data frame
##' @param reTrms Existing reTrms object from mkReTrms
##' @return Modified reTrms object with structure information
##' @keywords internal
integrate_structures_into_reTrms <- function(formula, data, reTrms) {
    
    # Parse formula for covariance structures using existing S4 infrastructure

        s4_object_list <- parse_model_formula(formula, data)
        
        if (length(s4_object_list) > 0) {
            # Store structure information in reTrms
            reTrms <- store_structure_info(reTrms, s4_object_list)
        }
        
        return(reTrms)
         
}

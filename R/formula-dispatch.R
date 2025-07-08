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
        dcov = list(prefix = "Diagonal"), 
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

    is_hom <- (is.null(args$hom) || args$hom == TRUE)
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
    specials_list <- c("ar1", "cs", "us", "dcov") 

    split_formula <- reformulas::splitForm(formula, specials = specials_list)

    s4_object_list <- list() 

    if (!is.null(split_formula$reTrmFormulas)) {

        temp_reTrms <- reformulas::mkReTrms(split_formula$reTrmFormulas, data, calc.lambdat = FALSE)

        for (i in seq_along(split_formula$reTrmClasses)) {

            type <- split_formula$reTrmClasses[i]
            add_args_call <- split_formula$reTrmAddArgs[[i]]
            cnms <- temp_reTrms$cnms[[i]]

            s4_object_list[[length(s4_object_list) + 1]] <- create_covariance_object_from_term(
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
    nth_per_term <- vapply(s4_object_list, n_parameters, integer(1))
    thoff <- cumsum(c(0L, nth_per_term))
    
    theta <- unlist(lapply(s4_object_list, get_start_values))
    lower <- unlist(lapply(s4_object_list, get_lower_bounds))

    lambdat_blocks <- list()
    lind_list <- list()

    for (i in seq_along(s4_object_list)) {
        s4_obj <- s4_object_list[[i]]
        nl <- reTrms$nl[i] # no. of levels for this term
        L_template <- get_lambda(s4_obj)
        
        # Create a sparse identity matrix of size nl x nl.
        # The Kronecker product then tiles L_template along the diagonal,
        # building the full sparse matrix for this term.
        lambdat_blocks[[i]] <- kronecker(Diagonal(nl), L_template) 

        lind_term <- get_lind(s4_obj)  

        # Generate the full Lind mapping for this term by replicating the
        # per-level index vector (lind_term) for each of the nl levels,
        # then add the global theta offset (thoff).
        if (length(lind_term) > 0) {
            lind_list[[i]] <- rep(lind_term, times = nl) + thoff[i] 
        } else {
            lind_list[[i]] <- integer(0)
        }
    }
    # Combine all the individual matrix blocks from the lambdat_blocks list
    # into the single, large block-diagonal Lambdat matrix for the entire model.
    Lambdat <- do.call(Matrix::bdiag, lambdat_blocks)
    Lind <- unlist(lind_list)
    
   if (length(Lambdat@x) != length(Lind)) {
        stop("Internal error: Mismatch between Lambdat@x length (", 
             length(Lambdat@x), ") and Lind length (", length(Lind), ")")
    }

    # Use the Lind index vector to populate the non-zero elements of the
    # Lambdat matrix with the corresponding values from the theta vector.
    Lambdat@x[] <- theta[Lind]
    
    list(
        Lambdat = Lambdat,
        Lind = as.integer(Lind),
        theta = theta,
        lower = lower
    )
}

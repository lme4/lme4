#' S4 Object Factory for Covariance Structures
#'
#' An internal helper function that translates the parsed components of a single
#' random effect term into the corresponding S4 covariance object.
#'
#' @param type A character string for the covariance structure type (e.g., "ar1", "us").
#' @param bar_formula The core `term | group` formula for this random effect.
#' @param add_args_call The language object of the wrapper function call (e.g., quote(ar1(d = 2))).
#' @param data The model data frame, used to determine the term's dimension.
#'
#' @return An S4 object inheriting from `VirtualCovariance`.
#' @keywords internal
create_covariance_object_from_term <- function(type, bar_formula, add_args_call, data) {

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

    is_hom <- (is.null(args$hom) || args$hom == TRUE)
    variance_prefix <- if(is_hom) "Homogeneous" else "Heterogeneous"

    full_class_name <- paste0(variance_prefix, structure_prefix, "Covariance")

    term_formula <- as.formula(paste("~", deparse(bar_formula[[2]])))
    dimension <- ncol(model.matrix(term_formula, data))

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
    browser()
    specials_list <- c("ar1", "cs", "us", "rr", "diag")
    split_formula <- reformulas::splitForm(formula, specials = specials_list)

    s4_object_list <- list() 

    if (!is.null(split_formula$reTrmFormulas)) {
        for (i in seq_along(split_formula$reTrmClasses)) {
    
            s4_object_list[[i]] <- create_covariance_object_from_term(
            type = split_formula$reTrmClasses[i],
            bar_formula = split_formula$reTrmFormulas[[i]],
            add_args_call = split_formula$reTrmAddArgs[[i]],
            data = data # 
            )
        }
    }
    return(s4_object_list)
}



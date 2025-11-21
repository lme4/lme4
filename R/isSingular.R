
getSingTol <- function() 
  getOption("lme4.singular.tolerance", 1e-4)

# Check for singular fit
#' @param object A fitted model.
#' @param ... Passed to methods
#' @export
isSingular <- function(object, ...) {
  UseMethod("isSingular")
}

#' @export
isSingular.merMod <- function(object, tol = getSingTol(), 
                              method = c("cholesky", "det"), ...) {
  method <- match.arg(method)
  
  if (method == "cholesky") {
    cc <- getME(object, "ST")
    diagvals <- lapply(cc, diag)
    min_diag <- min(unlist(diagvals))
    return(min_diag < tol)
  } else if (method == "det"){
    vc_list <- VarCorr(object)
    for (i in seq_along(vc_list)) {
      Sigma <- as.matrix(vc_list[[i]])
      if (length(Sigma) > 1) {
        if(det(Sigma) < tol) return(TRUE)
      }
    }
  }
  return(FALSE)
}


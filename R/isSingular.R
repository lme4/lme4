# Check for singular fit
#' @param object A fitted model.
#' @param ... Passed to methods
#' @export
isSingular <- function(object, ...) {
  UseMethod("isSingular")
}

#' @export
isSingular.merMod <- function(object, tol = 1e-6, method = c("eigen", "cholesky"), ...) {
  method <- match.arg(method)
  
  if (method == "cholesky") {
    cc <- getME(object, "ST")
    diagvals <- lapply(cc, diag)
    min_diag <- min(unlist(diagvals))
    return(min_diag < tol)
  }
  
  vc_list <- VarCorr(object)
  for (i in seq_along(vc_list)) {
    Sigma <- as.matrix(vc_list[[i]])
    if (length(Sigma) > 1) {
      eigvals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
      if (any(eigvals < tol)) return(TRUE)
    }
  }
  return(FALSE)
}


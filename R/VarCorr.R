##' Make variance and correlation matrices from \code{theta}
##'
##' @param sc scale factor (residual standard deviation)
##' @param cnms component names
##' @param nc numeric vector: number of terms in each RE component
##' @param theta theta vector (lower-triangle of Cholesky factors)
##' @param nms component names (FIXME: nms/cnms redundant: nms=names(cnms)?)
##' @seealso \code{\link{VarCorr}}
##' @return A matrix
##' @export
mkVarCorr <- function(sc, cnms, nc, theta, nms, reCovs = NULL) {
  if (is.null(reCovs)) {
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
  }
  else
    ncseq <- seq_along(reCovs)
  if(!all(nms == names(cnms))) ## the above FIXME
    warning("nms != names(cnms)  -- whereas lme4-authors thought they were --\n",
            "Please report!", immediate. = TRUE)
  ans <- lapply(ncseq, function(i)
  {
    ## Li := \Lambda_i, the i-th block diagonal of \Lambda(\theta)
    if (is.null(reCovs)) {
      Li <- diag(nrow = nc[i])
      Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
    }
    else
      Li <- getLambda(reCovs[[i]])
    rownames(Li) <- cnms[[i]]
    ## val := \Sigma_i = \sigma^2 \Lambda_i \Lambda_i', the
    val <- tcrossprod(sc * Li) # variance-covariance
    stddev <- sqrt(diag(val))
    corr <- t(val / stddev)/stddev
    diag(corr) <- 1
    structure(val, stddev = stddev, correlation = corr)
  })
  for(j in seq_along(ans)){
    cls <- sub("Covariance\\.", "", class(reCovs[[j]]))
    class(ans[[j]]) <- c(paste0("vcmat_", cls), class(ans[[j]]))
  }
  if(is.character(nms)) {
    ## FIXME: do we want this?  Maybe not.
    ## Potential problem: the names of the elements of the VarCorr() list
    ##  are not necessarily unique (e.g. fm2 from example("lmer") has *two*
    ##  Subject terms, so the names are "Subject", "Subject".  The print method
    ##  for VarCorrs handles this just fine, but it's a little awkward if we
    ##  want to dig out elements of the VarCorr list ... ???
    if (anyDuplicated(nms))
      nms <- make.names(nms, unique = TRUE)
    names(ans) <- nms
  }
  structure(ans, sc = sc)
}

##' Extract variance and correlation components
##'
VarCorr.merMod <- function(x, sigma = 1, ...)
{
  ## TODO: now that we have '...', add  type=c("varcov","sdcorr","logs" ?
  if (is.null(cnms <- x@cnms))
    stop("VarCorr methods require reTrms, not just reModule")
  if(missing(sigma))
    sigma <- sigma(x)
  nc <- lengths(cnms) # no. of columns per term
  structure(mkVarCorr(sigma, cnms = cnms, nc = nc, theta = x@theta,
                      nms = { fl <- x@flist; names(fl)[attr(fl, "assign")]},
                      reCovs = getReCovs(x)),
            useSc = as.logical(x@devcomp$dims[["useSc"]]),
            class = "VarCorr.merMod")
}

##' @S3method print VarCorr.merMod
print.VarCorr.merMod <- function(x, digits = max(3, getOption("digits") - 2),
                                 comp = "Std.Dev.", formatter = format, ...) {
  print(formatVC(x, digits=digits, comp=comp, formatter=formatter),
        quote = FALSE)
  invisible(x)
}

##' format columns corresponding to std. dev. and/or variance
##' @param reStdDev a vector of standard deviations
##' @param use.c a character vector indicating which scales to include
##' @param digits number of significant digits
##' @param formatter formatting function
##' @param ... additional arguments to formatter
##' @examples
##' invisible(gt_load("test_data/models.rda"))
##' format_sdvar(reStdDev = 1:3, use.c = c("Variance", "Std.Dev."))
##' format_sdvar(attr(VarCorr(fm1)$cond$Subject, "stddev"))
##' @export
## FIXME: avoid repeating defaults
format_sdvar <- function(reStdDev, use.c = "Std.Dev.", formatter=format,
                         digits = max(3, getOption("digits") - 2), ...) {
  res <- list()
  if("Variance" %in% use.c)
    res <- c(res,
             list(Variance = formatter(unlist(reStdDev)^2, digits = digits, ...)))
  if("Std.Dev." %in% use.c)
    res <- c(res, list(`Std.Dev.`=formatter(unlist(reStdDev),   digits = digits, ...)))
  mat <- do.call(cbind, res)
  colnames(mat) <- names(res)
  rownm <- names(res[[1]])  %||% ""
  mat <- cbind(Name = rownm, mat)
  rownames(mat) <- NULL
  return(mat)
}


##' @rdname format_sdvar
##' @param x a square numeric matrix
##' @param maxdim maximum number of rows/columns to display
##' @param digits digits for format
##' @param maxlen maximum number of rows to display
##' @param ... additional parameters
## FIXME: avoid repeating defaults
##' @export
format_corr <- function(x, maxdim=Inf, digits=2, maxlen = 10, ...) {
  UseMethod("format_corr")
}

##' @rdname format_sdvar
##' @export
get_sd <- function(x, ...) {
  UseMethod("get_sd")
}

##' @export
get_sd.default <- function(x, ...) {
  attr(x, "stddev")
}

##' @export
get_sd.vcmat_ar1 <- function(x, ...) {
  attr(x, "stddev")[1]
}

##' @export
format_corr.default <- function(x, maxdim = Inf, digits=2, ...) {
  if (length(x)==0) return("")
  x <- attr(x, "correlation")
  x <- as(x, "matrix")
  extra_rows <- (nrow(x) > maxdim)
  newdim <- min(maxdim, nrow(x))
  if (identical(c(x), NaN)) {
    cc <- matrix("(not stored)")
  } else {
    cc <- format(round(x, digits), nsmall = digits)
    cc[upper.tri(cc, diag = TRUE)] <- ""  ## empty upper triangle
    if (extra_rows) cc <- rbind(cc, "...")
  }
  cc
}

#' @export
format_corr.vcmat_diag <- function(x, maxdim = Inf, digits=2, ...) {
  ## empty correlation
  return(matrix(""))
}

#' @export
format_corr.vcmat_ar1 <- function(x, maxdim = Inf, digits=2, ...) {
  x <- attr(x, "correlation")
  if (length(x)==1) {
    cc <- format(round(x, digits), nsmall = digits)
  } else {
    cc <- format(round(x[2,1], digits), nsmall = digits)
  }
  return(matrix(paste(cc, "(ar1)")))
}

#' @export
format_corr.vcmat_cs <- function(x, maxdim = Inf, digits=2, ...) {
  x <- attr(x, "correlation")
  cc <- format(round(x[2,1], digits), nsmall = digits)
  return(matrix(paste(cc, "(cs)")))
}

## FIXME: get specials for ou, compsymm, spatial matrices, etc..

## original from lme4, *heavily* modified
##' "format()" the 'VarCorr' matrix of the random effects -- for
##' print()ing and show()ing
##'
##' @title Format the 'VarCorr' Matrix of Random Effects
##' @param varcor a \code{\link{VarCorr}} (-like) matrix with attributes.
##' @param digits the number of significant digits for standard deviations and variances
##' @param corr_digits the number of significant digits for correlations
##' @param comp character vector of length one or two indicating which
##' columns out of "Variance" and "Std.Dev." should be shown in the
##' formatted output.
##' @param formatter the \code{\link{function}} to be used for
##' formatting the standard deviations and or variances (but
##' \emph{not} the correlations which (currently) are always formatted
##' as "0.nnn"
##' @param useScale whether to report a scale parameter (e.g. residual standard deviation)
##' @param maxdim maximum dimensions (numbers of standard deviations/variances and number of
##' rows of correlation matrices) to report per random effects term
##' @param ... optional arguments for \code{formatter(*)} in addition
##' to the first (numeric vector) and \code{digits}.
##' @return a character matrix of formatted VarCorr entries from \code{varcor}.
##' @importFrom methods as
formatVC <- function(varcor, digits = max(3, getOption("digits") - 2),
                     corr_digits = max(2, digits-2),
                     maxdim = 10,
                     comp = "Std.Dev.", formatter = format,
                     useScale = attr(varcor, "useSc"),
                     ...)
{
  comp_opts <- c("Variance", "Std.Dev.")
  if(anyNA(mcc <- pmatch(comp, comp_opts))) {
    stop("Illegal 'comp': ", comp[is.na(mcc)])
  }
  use.c <- comp_opts[mcc]
  if (length(use.c) == 0) {
    stop("Must report either standard deviations or variances")
  }
  
  termnames <- names(varcor)
  
  ## get std devs (wait until after processing useScale to create output matrices)
  ## ugh, want to restrict lengths of sd that get reported: do we need methods/special
  ##   cases for this as well?
  
  reStdDev <- lapply(varcor, get_sd)
  
  ## get corr outputs
  corr_out <- lapply(varcor, format_corr, digits = corr_digits, maxdim = maxdim)
  
  if(useScale) {
    reStdDev <- c(reStdDev,
                  list(Residual = unname(attr(varcor, "sc"))))
    termnames <- c(termnames, "Residual")
    ## dummy correlation for Residual
    corr_out <- c(corr_out, list(matrix("")))
  }
  
  ## in order to get everything formatted consistently we have to collapse the std devs to a single
  ## vector, format them all at once, then split them back up (e.g. to insert extra spaces where necessary)
  
  trunc_rows <- function(x) {
    if (nrow(x) > maxdim) {
      x <- rbind(x[1:maxdim,,drop = FALSE], rep("...", ncol(x)))
    }
    return(x)
  }
  
  formatted_sdvar <- format_sdvar(unlist(unname(reStdDev)), digits = digits, comp = comp_opts, formatter = formatter, use.c = use.c)
  ## split back into chunks
  sdvar_out <- split.data.frame(formatted_sdvar,
                                rep(seq(length(reStdDev)), lengths(reStdDev)))
  sdvar_out <- lapply(sdvar_out, trunc_rows)
  
  names(sdvar_out) <- names(reStdDev)
  
  ## stick it all back together, properly spaced
  assemble_sdcor(sdvar_out, corr_out, termnames)
}

pad_blank <- function(m, max_rows=0, max_cols=0) {
  m <- as.matrix(m) ## handle scalar case
  if ((xrows <- (max_rows - nrow(m))) > 0) {
    m <- rbind(m, matrix("", nrow = xrows, ncol = ncol(m)))
  }
  if ((xcols <- (max_cols - ncol(m))) > 0) {
    m <- cbind(m, matrix("", ncol = xcols, nrow = nrow(m)))
  }
  return(m)
}

## patch together sd/var info, correlation info, group names
assemble_sdcor <- function(sdvar_out, corr_out, termnames) {
  
  sdvar_rows <- vapply(sdvar_out, nrow, numeric(1))
  corr_rows <- vapply(corr_out, nrow, numeric(1))
  max_rows <- pmax(sdvar_rows, corr_rows)
  
  nt <- length(corr_out)
  corr_cols <- vapply(corr_out, ncol, numeric(1))
  max_cols <- rep(max(corr_cols), nt)
  
  termnames_out <- mapply(pad_blank, termnames, max_rows, SIMPLIFY = FALSE)
  termnames_out <- do.call(rbind, termnames_out)
  colnames(termnames_out) <- "Groups"
  
  sdvar_out <- mapply(pad_blank, sdvar_out, max_rows, max_cols, SIMPLIFY = FALSE)
  sdvar_out <- do.call(rbind, sdvar_out)
  
  corr_out <- mapply(pad_blank, corr_out, max_cols = max_cols, max_rows = max_rows, SIMPLIFY = FALSE)
  corr_out <- do.call(rbind, corr_out)
  if (all(corr_out == "")) {
    corr_out <- NULL
  } else {
    colnames(corr_out) <- c("Corr", rep("", ncol(corr_out)-1))
  }
  ## FIXME: should we enable column names here? spacing, abbrev, etc to worry about
  ##  (first, making sure that null correlation matrices are unnamed)
  
  res <- cbind(termnames_out, sdvar_out, corr_out)
  rownames(res) <- rep("", nrow(res))
  
  return(res)
  
}

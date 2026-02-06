##' Make variance and correlation matrices from \code{theta}
##'
##' @param sc scale factor (residual standard deviation)
##' @param cnms component names
##' @param nc numeric vector: number of terms in each RE component
##' @param theta theta vector (lower-triangle of Cholesky factors)
##' @param nms component names (FIXME: nms/cnms redundant: nms=names(cnms)?)
##' @param full_cor specifies whether the full correlation matrix should be produced
##' @param is_lmm specifies whether the information came from a LMM
##' @seealso \code{\link{VarCorr}}
##' @return A matrix
##' @export
mkVarCorr <- function(sc, cnms, nc, theta, nms, reCovs = NULL, 
                      full_cor = NULL, is_lmm = NULL) {
  #browser()
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
    # the scalar only applies for LMMs
    if(!is.null(is_lmm) && is_lmm == FALSE){sc = 1}
    # val := \Sigma_i = \sigma^2 \Lambda_i \Lambda_i', the
    val <- tcrossprod(sc * Li) # variance-covariance 
    stddev <- sqrt(diag(val))
    ## if null, then only print the covariance matrix if nc <= 20.
    ## 20 is an arbitrary threshold, but anything larger than that would likely
    ## be too tedious to work with.
    if((is.null(full_cor) && nc[[i]] <= 20)
        || inherits(reCovs[[i]], "Covariance.ar1") 
        || full_cor == TRUE){
      corr <- t(val / stddev)/stddev
      diag(corr) <- 1
    } else {
      corr <- matrix(NaN)
    }
    ## adding more information depending on the covariance type
    if(inherits(reCovs[[i]], "Covariance.ar1") 
       || inherits(reCovs[[i]], "Covariance.cs")){
      structure(val, stddev = stddev, correlation = corr,
                theta = getTheta(reCovs[[i]]), profpar = getProfPar(reCovs[[i]]),
                rho = getVC(reCovs[[i]])$ccomp)
    } else {
      structure(val, stddev = stddev, correlation = corr,
                theta = getTheta(reCovs[[i]]), profpar = getProfPar(reCovs[[i]]))
    }
  })
  for(j in seq_along(ans)){
    reCov <- reCovs[[j]]
    cls <- sub("Covariance\\.", "", class(reCov))
    ## There are separate vcmat_hetar1, vcmat_homcs classes
    ## Inconsistent naming because AR1 default is typically homogenous,
    ## CS default is heterogeneous (cf glmmTMB)
    hom_status <- ""
    if ("hom" %in% slotNames(reCov)) {
      if (!reCov@hom && inherits(reCov, "Covariance.ar1")) hom_status <- "het"
      if (reCov@hom && inherits(reCov, "Covariance.cs")) hom_status <- "hom"
    }
    class(ans[[j]]) <- c(paste0("vcmat_", hom_status, cls), class(ans[[j]]))
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

## FIXME: automate this from list of known Covariance.* classes ... 
for (varclass in c("us",
                   c(outer(c("hom", "het"), c("ar1", "cs", "diag"),
                           function(x, y) paste(x, y, sep = "_"))))) {
  setOldClass(c(paste0("vcmat_", varclass), "matrix", "array"))
}
                             
##' Extract variance and correlation components
##'
VarCorr.merMod <- function(x, sigma = 1, full_cor = NULL, ...)
{
  ## TODO: now that we have '...', add  type=c("varcov","sdcorr","logs" ?
  if (is.null(cnms <- x@cnms))
    stop("VarCorr methods require reTrms, not just reModule")
  if(missing(sigma))
    sigma <- sigma(x)
  nc <- lengths(cnms) # no. of columns per term
  structure(mkVarCorr(sigma, cnms = cnms, nc = nc, theta = x@theta,
                      nms = { fl <- x@flist; names(fl)[attr(fl, "assign")]},
                      reCovs = getReCovs(x),
                      full_cor = full_cor,
                      is_lmm = isLMM(x)),
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


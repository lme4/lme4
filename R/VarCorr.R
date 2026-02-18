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
mkVarCorr <-
function(sc, cnms, nc, theta, nms, reCovs = NULL,
         full_cor = NULL, is_lmm = NULL) {
    if (!identical(names(cnms), nms)) # FIXME? deprecate 'nms'
        warning("names(cnms) and 'nms' are not identical; please report",
                immediate. = TRUE)
    nc <- as.integer(nc)
    if (is.null(reCovs))
        reCovs <- .mapply(new,
                          list(nc = nc,
                               par = split(theta, rep(seq_along(nc), (nc * (nc + 1L)) %/% 2L))),
                          list(Class = "Covariance.us"))
    ans <- lapply(seq_along(reCovs), function(i) {
        ## Li := \Lambda_i, the i-th diagonal block of \Lambda(\theta)
        ## Si := \sigma^2 \Lambda_i \Lambda_i'
        object <- reCovs[[i]]
        nci <- nc[i]
        Li <- getLambda(object)
        rownames(Li) <- cnms[[i]]
        if (!(is.null(is_lmm) || is_lmm))
            sc <- 1 # nonunit only for LMMs
        jj <- seq.int(from = 1L, by = nci + 1L, length.out = nci)
        Si <- sc * sc * tcrossprod(Li)
        Si.sd <- sqrt(Si[jj])
        if (full_cor %||% nci <= 20L) { # FIXME? condition on structure
            Si.cor <- Si/Si.sd/rep(Si.sd, each = nci)
            Si.cor[jj] <- 1
        }
        else
            Si.cor <- matrix(NaN)
        typ0 <- typ1 <- sub("^Covariance[.]", "", class(object))
        switch(typ0, # respecting the conventions of, e.g., glmmTMB:
               "cs"  = if ( object@hom) typ1 <- "homcs",
               "ar1" = if (!object@hom) typ1 <- "hetar1")
        class(Si) <- c(paste0("vcmat_", typ1), "matrix", "array")
        attr(Si, "stddev") <- Si.sd
        attr(Si, "correlation") <- Si.cor
        attr(Si, "theta") <- getTheta(object)
        attr(Si, "profpar") <- getProfPar(object)
        if (typ0 %in% c("cs", "ar1"))
        attr(Si, "rho") <- if (nci > 1L) getVC(object)$ccomp else NaN
        Si
    })
    if (is.character(nms)) {
        ## FIXME:
        ## Do we want this?  Maybe not.  'nms' are not necessarily
        ## unique, e.g., 'fm2' from example("lmer") has *two* Subject
        ## terms, so the names are "Subject", "Subject".  The 'print'
        ## method for 'VarCorr.merMod' handles this just fine, but
        ## it's a little awkward if we want to dig out components of
        ## the list ...
        if (anyDuplicated(nms))
            nms <- make.names(nms, unique = TRUE)
        names(ans) <- nms
    }
    attr(ans, "sc") <- sc
    ans
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

as.data.frame.VarCorr.merMod <- function(x,row.names = NULL,
                                         optional = FALSE,
                                         order = c("cov.last", "lower.tri"),
                                         ...)  {
    order <- match.arg(order)
    tmpf <- function(v,grp) {
        vcov <- c(diag(v), v[lt.v <- lower.tri(v, diag = FALSE)])
        sdcor <- c(attr(v,"stddev"),
                   attr(v,"correlation")[lt.v])
        nm <- rownames(v)
        n <- nrow(v)
        dd <- data.frame(grp = grp,
                         var1 = nm[c(seq(n), col(v)[lt.v])],
                         var2 = c(rep(NA,n), nm[row(v)[lt.v]]),
                         vcov,
                         sdcor,
                         stringsAsFactors = FALSE)
        if (order=="lower.tri") {
            ## reorder *back* to lower.tri order
            m <- matrix(NA,n,n)
            diag(m) <- seq(n)
            m[lower.tri(m)] <- (n+1):(n*(n+1)/2)
            dd <- dd[m[lower.tri(m, diag=TRUE)],]
        }
        dd
    }
    r <- do.call(rbind,
                 c(mapply(tmpf, x,names(x), SIMPLIFY = FALSE),
                   deparse.level = 0))
    if (attr(x,"useSc")) {
        ss <- attr(x,"sc")
        r <- rbind(r,data.frame(grp = "Residual",var1 = NA,var2 = NA,
                                vcov = ss^2,
                                sdcor = ss),
                   deparse.level = 0)
    }
    rownames(r) <- NULL
    r
}

##' @S3method print VarCorr.merMod
print.VarCorr.merMod <- function(x, digits = max(3, getOption("digits") - 2),
                                 comp = "Std.Dev.", formatter = format, ...) {
  print(formatVC(x, digits=digits, comp=comp, formatter=formatter),
        quote = FALSE)
  invisible(x)
}


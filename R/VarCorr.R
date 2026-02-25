mkVarCorr <-
function(sc, cnms, nc = lengths(cnms, use.names = FALSE),
         theta, nms = names(cnms), reCovs = NULL) {
    if (!missing(nc)) {
        nc <- as.integer(nc)
        if (!missing(cnms) && !identical(nc, lengths(cnms, use.names = FALSE)))
            stop("'nc' and lengths(cnms) are inconsistent")
    }
    if (!missing(nms)) {
        if (!missing(cnms) && !identical(nms, names(cnms)))
            stop("'nms' and names(cnms) are inconsistent")
    }
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
        jj <- seq.int(from = 1L, by = nci + 1L, length.out = nci)
        Si <- if (is.null(sc)) tcrossprod(Li) else sc * sc * tcrossprod(Li)
        Si.sd <- sqrt(Si[jj])
        Si.cor <- Si/Si.sd/rep(Si.sd, each = nci)
        Si.cor[jj] <- 1
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
    if (!is.null(nms)) {
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

## MJ:
## Maintain this list manually as prefixing convention is inconsistent.
for (.nm in c("us", "diag", "cs", "homcs", "ar1", "hetar1"))
    setOldClass(c(paste0("vcmat_", .nm), "matrix", "array"))
rm(.nm)
                             
VarCorr.merMod <-
function(x, sigma = 1, ...) {
    ## TODO? add argument type=c("varcov", "sdcor", "logs")
    useSc <- isLMM(x)
    sc <- if (useSc) { if (missing(sigma)) sigma(x) else sigma }
    cnms <- x@cnms
    nc <- lengths(cnms, use.names = FALSE)
    theta <- x@theta
    nms <- names(x@flist)[attr(x@flist, "assign")]
    reCovs <- getReCovs(x)
    ans <- mkVarCorr(sc = sc, cnms = cnms, nc = nc,
                     theta = theta, nms = nms, reCovs = reCovs)
    class(ans) <- "VarCorr.merMod"
    attr(ans, "useSc") <- useSc # for reformulas::formatVC
    ans
}

as.data.frame.VarCorr.merMod <-
function(x, row.names = NULL, optional = FALSE,
         order = c("cov.last", "lower.tri"), ...) {
    order <- match.arg(order)
    dlist <- lapply(seq_along(x), function(i) {
        Si <- x[[i]]
        Si.sd <- attr(Si, "stddev")
        Si.cor <- attr(Si, "correlation")
        nci <- ncol(Si)
        cnmsi <- colnames(Si)
        jj <- seq.int(from = 1L, by = nci + 1L, length.out = nci)
        ij1 <- sequence.default(from = jj, nvec = nci:1L)
        ij0 <- sequence.default(from = jj + 1L, nvec = (nci - 1L):0L)
        if (order == "cov.last")
            data.frame(grp = names(x)[i],
                       var1 = cnmsi[c(1L:nci,
                                      1L + (ij0 - 1L) %/% nci)],
                       var2 = cnmsi[c(rep(NA_integer_, nci),
                                      1L + (ij0 - 1L) %% nci)],
                       vcov = Si[c(jj, ij0)],
                       sdcor = c(Si.sd, Si.cor[ij0]),
                       row.names = NULL,
                       stringsAsFactors = FALSE)
        else {
            jj. <- cumsum(c(1L, if (nci > 1L) nci:2L))
            data.frame(grp = names(x)[i],
                       var1 = cnmsi[1L + (ij1 - 1L) %/% nci],
                       var2 = `[<-`(cnmsi[1L + (ij1 - 1L) %% nci],
                                    jj., NA_character_),
                       vcov = Si[ij1],
                       sdcor = `[<-`(Si.cor[ij1], jj., Si.sd),
                       row.names = NULL,
                       stringsAsFactors = FALSE)
        }
    })
    if (!is.null(sc <- attr(x, "sc")) && (attr(x, "useSc") %||% TRUE)) {
        r <- data.frame(grp = "Residual",
                        var1 = NA_character_,
                        var2 = NA_character_,
                        vcov = sc * sc,
                        sdcor = sc,
                        row.names = NULL,
                        stringsAsFactors = FALSE)
        dlist <- c(dlist, list(r))
    }
    do.call(rbind, dlist)
}

print.VarCorr.merMod <-
function(x, digits = max(3L, getOption("digits") - 2L),
         comp = "Std.Dev.", formatter = format, ...) {
    y <- formatVC(x, digits = digits, comp = comp, formatter = formatter)
    print(y, quote = FALSE)
    invisible(x)
}

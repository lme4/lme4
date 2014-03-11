##' Flexible lmer
##'
##' @param formula Formula
##' @param data Data
##' @param specials Special functions in \code{formula}
##' @param verbose Verbose
##' @return \code{merMod} object
flexLmer <- function(formula, data, specials = c("cs", "d", "ar1d"), verbose = 0L) {
    # see example(modular)
    lmod <- flexFormula(formula, data, family = NULL, specials = specials, verbose = verbose)
    devfun <- do.call(mkLmerDevfun, lmod)
    opt <- optimizeLmer(devfun, verbose = verbose)
    # FIXME: this should not be necessary...
    environment(devfun)$pp$theta <- opt$par
    list(model = mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr), 
        opt = opt, devfun = devfun)
}

flexGlmer <- function(formula, data, family = gaussian(), specials = c("cs", "d", "ar1d"), verbose = 0L) {
    # see example(modular)
    lmod <- flexFormula(formula, data, family, specials = specials, verbose = verbose)
    devfun <- do.call(mkGlmerDevfun, lmod)
    opt <- optimizeGlmer(devfun, verbose = verbose)
    # FIXME: this should not be necessary...
    environment(devfun)$pp$theta <- opt$par
    list(model = mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr), 
        opt = opt, devfun = devfun)

}

flexFormula <- function(formula, data, family = NULL, specials = c("cs","d","ar1d"), verbose = 0L) {
    # split off reGenerator terms:
    frmlterms <- terms(formula, specials = specials)
    termnames <- attr(frmlterms, "variables")
    where.specials <- unlist(attr(frmlterms, "specials")) - 1  #offset by one since response is counted
    lmerformula <- formula(paste(formula[[2]], "~", paste(attr(frmlterms, "term.labels")[-where.specials],
        collapse = "+")))
    reGenerators <- as.formula(paste("~", paste(attr(frmlterms, "term.labels")[where.specials], 
        collapse = "+")))
    # see example(modular)
    if(is.null(family)) {
        return( lFormula(lmerformula, data, reGenerators = reGenerators))
    } else {
        return(glFormula(lmerformula, data, reGenerators = reGenerators, family = family))
    }
}

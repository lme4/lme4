##' Flexible lmer
##'
##' @param formula Formula
##' @param data Data
##' @param specials Special functions in \code{formula}
##' @param verbose Verbose
##' @return \code{flexMerMod} object
flexLmer <- function(formula, data, specials = c("cs", "d", "ar1d"),
                     verbose = 0L, control = lmerControl(),
                     weights = NULL) {
    # see example(modular)
    lmod <- flexFormula(formula, data, family = NULL, specials = specials,
                        verbose = verbose, control = control, weights = weights)
    devfun <- do.call(mkLmerDevfun,
		      c(lmod,
			list(verbose=verbose,control=control)))
    opt <- optimizeLmer(devfun, verbose = verbose)
    # FIXME: this should not be necessary...
    environment(devfun)$pp$theta <- opt$par
    out <- list(model = mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr), 
                opt = opt, devfun = devfun, reGenerators = lmod$reGenerators)
    class(out) <- "flexMerMod"
    return(out)
}

flexGlmer <- function(formula, data, family = gaussian(), specials = c("cs", "d", "ar1d"), verbose = 0L) {
    # see example(modular)
    lmod <- flexFormula(formula, data, family, specials = specials, verbose = verbose)
    devfun <- do.call(mkGlmerDevfun, lmod)
    opt <- optimizeGlmer(devfun, verbose = verbose)
    # FIXME: this should not be necessary...
    environment(devfun)$pp$theta <- opt$par
    out <- list(model = mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr), 
                opt = opt, devfun = devfun, reGenerators = lmod$reGenerators)
    class(out) <- "flexMerMod"
    return(out)
}

splitregen <- function(formula, specials = c("d", "cs", "ar1d")) {

    if(TRUE) { ## new procedure for spliting specials
        
        ## Recursive function: (f)ind (b)ars (a)nd (s)pecials
        ## cf. fb function in findbars
        fbas <- function(term) {
            if (is.name(term) || !is.language(term)) return(NULL)
            for (sp in specials) if (term[[1]] == as.name(sp)) return(term)
            if (term[[1]] == as.name("(")) return(term)
            stopifnot(is.call(term))
            if (term[[1]] == as.name('|')) return(term)
            if (length(term) == 2) return(fbas(term[[2]]))
            c(fbas(term[[2]]), fbas(term[[3]]))
        }
        splits <- fbas(formula)
        splitID <- sapply(lapply(splits, "[[", 1), as.character)
        splitStan <- splits[splitID == "("]
        splitSpec <- splits[!(splitID == "(")]

        if(length(splitSpec) == 0) stop(
                     "no special covariance structures",
                     "please use lmer, not flexLmer")
        reGenerators <- as.formula(paste("~ ", paste(splitSpec, collapse = " + ")))   

        if(length(splitStan) == 0) {
            lmerformula <- nobars(formula)
        } else {
            lmerformula <- formula(paste(formula[[2]],                        ## response
                                         "~",
                                         as.character(nobars(formula))[[3]],  ## fixed effects
                                         " + ",
                                         paste(splitStan, collapse = " + "))) ## standard random effects
        }
        
    }
    else {
        frmlterms <- terms(formula, specials = specials)
        termnames <- attr(frmlterms, "variables")
        where.specials <- unlist(attr(frmlterms, "specials")) - 1  # offset by one since response is counted
        lmerformula <- formula(paste(formula[[2]], "~",
                                     paste(attr(frmlterms, "term.labels")[-where.specials],
                                           collapse = "+")))
        reGenerators <- as.formula(paste("~",
                                         paste(attr(frmlterms, "term.labels")[where.specials], 
                                               collapse = "+")))
    }
    return(list(lmerformula = lmerformula,
                reGenerators = reGenerators))
}


flexFormula <- function(formula, data, family = NULL, specials = c("cs","d","ar1d"),
                        verbose = 0L, control, weights = NULL) {
                                        # split off reGenerator terms:
    splt <- splitregen(formula, specials)
                                        # see example(modular)
    if(is.null(family)) {
        if(missing(control)) control = lmerControl()
        return( lFormula(splt$lmerformula, data, reGenerators = splt$reGenerators,
                         control = control, weights = weights))
    } else {
        if(missing(control)) control = glmerControl()
        return(glFormula(splt$lmerformula, data, reGenerators = splt$reGenerators,
                         control = control, weights = weights,
                         family = family))
    }
}


print.flexMerMod <- function(x, ...) {
    print(x$model, ...)
    lapply(x$reGenerators, function(reG) {
        print(environment(reG)$RETypeName)
    })
    invisible(NULL)
}

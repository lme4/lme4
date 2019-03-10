## convert vcov from dpoMatrix to regular matrix (protect against car methods)
vv <- function(x) Matrix::as.matrix(vcov(x))

## minimal "influence" function (to make broom::augment_columns work)
influence.merMod <- function(model, groups, data, maxfun=1000, do.coef = TRUE,
                             ...) {

    .groups <- NULL  ## avoid false-positive code checks
    
    if (length(list(...))>0) warning("disregarded additional arguments")
    if (!do.coef) {
        ## simple/quick/trivial results
        result <- list(hat=hatvalues(model))
        class(result) <- "influence.merMod"
        return(result)
    }
    if (missing(data)){
        data <- getCall(model)$data
        data <- if (!is.null(data)) eval(data, parent.frame())
                else stop("model did not use the data argument")
    }
    if (missing(groups)) {
        groups <- ".case"
        data$.case <- rownames(data)
    }
    else if (length(groups) > 1){
        del.var <- paste0(groups, collapse=".")
        data[, del.var] <- apply(data, 1, paste0, collapse=".")
        groups <- del.var
    }
    unique.del <- unique(data[, groups])
    data[[".groups"]] <- data[, groups]
    par <- list(theta=getME(model, "theta"))
    if (inherits(model, "glmerMod")) par$fixef <- fixef(model)
    fixed <- fixef(model)
    fixed.1 <- matrix(0, length(unique.del), length(fixed))
    rownames(fixed.1) <- unique.del
    colnames(fixed.1) <- names(fixed)
    Vs <- VarCorr(model)
    nms <- names(Vs)
    sep <- ":"
    if (length(nms) == 1) {
        nms <- ""
        sep <- ""
    }
    vc <- getME(model, "sigma")^2
    names(vc) <- "sigma^2"
    for (i in 1:length(Vs)){
        V <- Vs[[i]]
        c.names <- colnames(V)
        e.names <- outer(c.names, c.names, function(a, b) paste0("C[", a, ",", b, "]"))
        diag(e.names) <- paste0("v[", c.names, "]")
        v <- V[lower.tri(V, diag=TRUE)]
        names(v) <- paste0(nms[i], sep, e.names[lower.tri(e.names, diag=TRUE)])
        vc <- c(vc, v)
    }
    vc.1 <- matrix(0, length(unique.del), length(vc))
    rownames(vc.1) <- unique.del
    colnames(vc.1) <- names(vc)
    feval <- numeric(length(unique.del))
    converged <- logical(length(unique.del))
    vcov.1 <- vector(length(unique.del), mode="list")
    names(vcov.1) <- names(feval) <- names(converged) <- unique.del
    # control <- if (one.step){
    #     if (inherits(model, "lmerMod")) lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B", maxit=1))
    #     else if (inherits(model, "glmerMod")) glmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B", maxit=1))
    # } else {
    control <- if (inherits(model, "lmerMod")) lmerControl(optCtrl=list(maxfun=maxfun))
        else if (inherits(model, "glmerMod")) glmerControl(optCtrl=list(maxfun=maxfun))
                                        # }
    ## FIXME: parallelize?
    for (del in unique.del){
        data$del <- del
        mod.1 <- suppressWarnings(update(model, data=data,
                                         subset=(.groups != del),
                                         start=par,
                                         control=control))

        opt <- mod.1@optinfo
        feval[del] <- opt$feval
        converged[del] <- opt$conv$opt == 0 && length(opt$warnings) == 0
        fixed.1[del, ] <- fixef(mod.1)
        Vs.1 <- VarCorr(mod.1)
        vc.0 <- getME(mod.1, "sigma")^2
        for (V in Vs.1){
            vc.0 <- c(vc.0, V[lower.tri(V, diag=TRUE)])
        }
        vc.1[del, ] <- vc.0

        vcov.1[[del]] <- vv(mod.1)
    }
    left <- "[-"
    right <- "]"
    if (groups == ".case") {
        groups <- "case"
    }
    nms <- c("fixed.effects", paste0("fixed.effects", left, groups, right),
             "var.cov.comps", paste0("var.cov.comps", left, groups, right),
             "vcov", paste0("vcov", left, groups, right),
             "groups", "deleted", "converged", "function.evals")
    result <- list(fixed, fixed.1, vc, vc.1,
                   vv(model), vcov.1, groups, unique.del, converged, feval)
    names(result) <- nms
    class(result) <- "influence.merMod"
    result
}

## lookup table for influence.merMod elements
ipos <- c("fixed"=1, "fixed.sub"=2,
          "var.cov"=3, "var.cov.sub"=4,
          "vcov"=5, "vcov.sub"=6)

dfbeta.influence.merMod <- function(model, which=c("fixed", "var.cov"), ...){
    which <- match.arg(which)
    b <- model[[ipos[sprintf("%s.sub",which)]]]
    b0 <- model[[ipos[which]]]
    b - matrix(b0, nrow=nrow(b), ncol=ncol(b), byrow=TRUE)
}

dfbetas.influence.merMod <- function(model, ...){
    vList <- model[[ipos["vcov.sub"]]]
    n <- nrow(vList[[1]])
    vmat <- t(vapply(vList, function(x) sqrt(diag(x)),
                     numeric(n)))
    dfbeta(model)/vmat
}

cooks.distance.merMod <- function(model, ...) {
    p <- Matrix::rankMatrix(getME(model,"X"))
    hat <- hatvalues(model)
    ## FIXME: check dispersion
    dispersion <- sigma(model)^2
    res <- residuals(model,type="pearson")
    res <-  (res/(1 - hat))^2 * hat/(dispersion * p)
    res[is.infinite(res)] <- NaN
    res
}

cooks.distance.influence.merMod <- function(model, ...) { 
    db <- dfbeta(model)
    n <- nrow(db)
    p <- ncol(db)
    d <- numeric(n)
    vcovs <- model[["vcov.sub"]]
    sig.sq <- model[["var.cov.sub"]][, 1]
    for (i in 1:n){
        d[i] <- (db[i, ] %*% solve(vcovs[[i]]) %*% db[i, ])/(p*sig.sq[i])
    }
    d
}

## from ?lm.influence:
##  wt.res: a vector of _weighted_ (or for class 'glm' rather _deviance_)
##          residuals.
##
## residuals.lm gives r*sqrt(object$weights) (if non-NULL weights)
##   for type %in% c("deviance","pearson")
##
## residuals.glm gives  (y - mu) * sqrt(wts)/sqrt(object$family$variance(mu))
##
rstudent.merMod <- function (model, ...) {
    r <- residuals(model, type="deviance")
    hat <- hatvalues(model)
    pr <- residuals(model, type="pearson")
    r <- sign(r) * sqrt(r^2 + (hat * pr^2)/(1 - hat))
    r[is.infinite(r)] <- NaN
    r/sigma(model)
}

##' @S3method sigma merMod
sigma.merMod <- function(object, ...) {
    dc <- object@devcomp
    dd <- dc$dims
    if(dd[["useSc"]])
        dc$cmp[[if(dd[["REML"]]) "sigmaREML" else "sigmaML"]] else 1.
}

##' @importFrom stats terms
##' @S3method terms merMod
terms.merMod <- function(x, fixed.only=TRUE, random.only=FALSE, ...) {
    if (missing(fixed.only) && random.only) fixed.only <- FALSE
    if (fixed.only && random.only) stop("can't specify 'only fixed' and 'only random' terms")
    tt <- attr(x@frame,"terms")
    if (fixed.only) {
        tt <- terms.formula(formula(x,fixed.only=TRUE))
        attr(tt,"predvars") <- attr(terms(x@frame),"predvars.fixed")
    }
    if (random.only) {
        tt <- terms.formula(subbars(formula(x,random.only=TRUE)))
        ## FIXME: predvars should be random-only
        attr(tt,"predvars") <- attr(terms(x@frame),"predvars.random")
    }
    tt
}

##' @importFrom stats update
##' @S3method update merMod
update.merMod <- function(object, formula., ..., evaluate = TRUE) {
    if (is.null(call <- getCall(object)))
        stop("object should contain a 'call' component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
        call$formula <- update.formula(formula(object), formula.)
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate) {
        ff <- environment(formula(object))
        pf <- parent.frame()  ## save parent frame in case we need it
        sf <- sys.frames()[[1]]
        tryCatch(eval(call, envir=ff),
                 error=function(e) {
                     tryCatch(eval(call, envir=sf),
                              error=function(e) {
                                  eval(call, pf)
                              })
                 })
    } else call
}

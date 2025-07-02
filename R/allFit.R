meth.tab.0 <- cbind(optimizer=
                      rep(c("bobyqa",
                            "Nelder_Mead",
                            "nlminbwrap",
                            "nmkbw",
                            "optimx",
                            "nloptwrap" ),
                          c(rep(1,5),2)),
                    method= c(rep("",4), "L-BFGS-B",
                            "NLOPT_LN_NELDERMEAD", "NLOPT_LN_BOBYQA"))

## ugh: hardcoded list (incomplete?) of allowable *control* options by optimizer
## could make more of an effort to match max-iterations/evals,
## (x|f)*(abs|rel) tolerance, ...
opt.ctrls <- list(bobyqa=c("npt","rhobeg","rhoend","iprint","maxfun"),
                  Nelder_Mead=c("iprint","maxfun","FtolAbs",
                                "FtolRel","XtolRel","MinfMax",
                                "xst","xt","verbose","warnOnly"),
                  nlminbwrap=c("eval.max","iter.max","trace","abs.tol",
                               "rel.tol","x.tol","xf.tol","step.min",
                               "step.max","sing.tol","scale.init",
                               "diff.g"),
                  optimx=c("trace","fnscale","parscale","ndeps",
                          "maxit","abstol","reltol","method"),
                  nloptwrap=c("minf_max","ftol_rel","ftol_abs",
                              "xtol_rel", "xtol_abs", "maxeval", "maxtime",
                              "algorithm"),
                  nmkbw=c("tol","maxfeval","regsimp","maximize",
                          "restarts.max","trace","maxfun"))

## name of 'max function evaluations' for each optimizer
maxfun_arg <- c(bobyqa = "maxfun",
                Nelder_Mead = "maxfun",
                nlminbwrap = "eval.max",
                optimx = "maxit",
                nloptwrap = "maxeval",
                nmkbw = "maxfeval")
                           
                
nmkbw <- function(fn, par, lower, upper, control) {
    if (length(par)==1) {
        res <- optim(fn=fn,par=par,lower=lower,upper=100*par,
                     method="Brent")
    } else {
        res <- dfoptim::nmkb(fn=fn,par=par,
                             lower=lower,upper=upper,control=control)
    }
    res$fval <- res$value
    res
}


##' Attempt to re-fit a [g]lmer model with a range of optimizers.
##' The default is to use all known optimizers for R that satisfy the
##' requirements (do not require explicit gradients, allow
##' box constraints), in three categories; (i) built-in
##' (minqa::bobyqa, lme4::Nelder_Mead, nlminbwrap), (ii) wrapped via optimx
##' (most of optimx's optimizers that allow box constraints require
##' an explicit gradient function to be specified; the two provided
##' here are really base R functions that can be accessed via optimx,
##' (iii) wrapped via nloptr; (iv)
##'
##' @param object a fitted model
##' @param meth.tab a matrix (or data.frame) with columns
##' - method  the name of a specific optimization method to pass to the optimizer
##'           (leave blank for built-in optimizers)
##' - optimizer  the \code{optimizer} function to use
##' @param verbose print progress messages?
##' @param catch.err catch errors?
##' @return a list of fitted \code{merMod} objects
##' @seealso slice, slice2D in the bbmle package
##' @examples
##' library(lme4)
##' gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'                  data = cbpp, family = binomial)
##' gm_all <- allFit(gm1, parallel=TRUE)
##' ss <- summary(gm_all)
##' ss$fixef               ## extract fixed effects
##' ss$llik                ## log-likelihoods
##' ss$sdcor               ## SDs and correlations
##' ss$theta               ## Cholesky factors
##' ss$which.OK            ## which fits worked

allFit <- function(object, meth.tab = NULL,
                   data=NULL,
                   verbose=TRUE,
                   show.meth.tab = FALSE,
                   maxfun = 1e5,
                   parallel = c("no", "multicore", "snow"),
                   ncpus = getOption("allFit.ncpus", 1L),
                   cl = NULL,
                   catch.errs = TRUE) {

    if (is.null(meth.tab)) {
        meth.tab <- meth.tab.0
    }
    optvec <- meth.tab.0[,"optimizer"]
    if (!requireNamespace("dfoptim", quietly = TRUE)) {
        optvec <- setdiff(optvec, "nmkbw")
        meth.tab <- meth.tab[meth.tab[,"optimizer"] %in% optvec, ]
    }
    if (!requireNamespace("optimx", quietly = TRUE)) {
        optvec <- setdiff(optvec, "optimx")
        meth.tab <- meth.tab[meth.tab[,"optimizer"] %in% optvec,]
    }
    if (show.meth.tab) {
        return(meth.tab)
    }
    stopifnot(length(dm <- dim(meth.tab)) == 2, dm[1] >= 1, dm[2] >= 2,
	      is.character(optimizer <- meth.tab[,"optimizer"]),
	      is.character(method    <- meth.tab[,"method"]))

    parallel <- match.arg(parallel)
    do_parallel <- have_mc <- have_snow <- NULL # "-Wall"
    eval(initialize.parallel) # (parallel, ncpus)   --> ./utilities.R
    ## |--> (do_parallel, have_mc, have_snow)

    fit.names <- gsub("\\.$", "", paste(optimizer, method, sep="."))
    ffun <- local({
        ## required local vars
        object
        verbose
        fit.names
        optimizer
        method
        maxfun
        function(..i) {
            if (verbose) cat(fit.names[..i],": ")
            ctrl <- getCall(object)$control
            ## NB:  'ctrl' must become a correct *argument* list for g?lmerControl()
            if(is.null(ctrl)) {
                ctrl <- list(optimizer=optimizer[..i])
            } else {
                if(is.call(ctrl)) # typically true
                    ctrl <- lapply(as.list(ctrl)[-1], eval)
                ctrl$optimizer <- optimizer[..i]
            }
            ## add method/algorithm to optCtrl as necessary
            mkOptCtrl <- function(...) {
                x <- list(...)
                cc <- ctrl$optCtrl
                if (is.null(cc)) cc <- list()
                cc[[maxfun_arg[[optimizer[[..i]]]]]] <- maxfun
                for (n in names(x)) { ## replace existing values
                    cc[[n]] <- x[[n]]
                }
                cc
            }
            sanitize <- function(x,okvals) {
                if (is.null(x)) return(NULL)
                if (is.null(okvals)) return(x)
                x <- x[intersect(names(x),okvals)]
                ## ?? having names(control) be character(0)
                ##  screws up nmkbw ... ??
                if (length(names(x))==0) names(x) <- NULL
                x
            }
            ctrl$optCtrl <- switch(optimizer[..i],
                                   optimx    = mkOptCtrl(method   = method[..i]),
                                   nloptwrap = mkOptCtrl(algorithm= method[..i]),
                                   mkOptCtrl())
                                       
            ctrl$optCtrl <- sanitize(ctrl$optCtrl,
                                     opt.ctrls[[optimizer[..i]]])
            ctrl <- do.call(if(isGLMM(object)) glmerControl else lmerControl, ctrl)
            ## need to stick ctrl in formula env so it can be found ...
            assign("ctrl", ctrl, environment(formula(object)))
            if (catch.errs) {
              tt <- system.time(rr <- tryCatch(update(object, control = ctrl),
                                               error = function(e) e))
            } else {
              tt <- system.time(rr <- update(object, control = ctrl))
            }
            attr(rr, "optCtrl") <- ctrl$optCtrl # contains crucial info here
            attr(rr, "time") <- tt  # store timing info
            if (verbose) {
              if (inherits(rr,"error")) {
                cat("[failed]\n")
              } else {
                cat("[OK]\n")
              }
            }
            return(rr)
        }
    })

    seq_fit <- seq_along(fit.names)
    res <- if (do_parallel) {
               if (have_mc) {
                   parallel::mclapply(seq_fit,
                                      ffun, mc.cores = ncpus)
               } else if(have_snow) {
                   if(is.null(cl)) {
                       cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                       ## consider exporting data/package ?
                       ## parallel::clusterEvalQ(cl,library("lme4"))
                       ## try to get data and export it?
                       ## parallel::clusterExport(cl,??)
                       res <- parallel::parLapply(cl, seq_fit, ffun)
                       parallel::stopCluster(cl)
                       res
                   } else parallel::parLapply(cl, seq_fit, ffun)
               } else {
                   warning("'do_parallel' is true, but 'have_mc' and 'have_snow' are not.  Should not happen!")
                   ## or stop()  or  we could silently use lapply(..)
                   setNames(as.list(fit.names), fit.names)
               }
           } else
               lapply(seq_fit, ffun)

    names(res) <- fit.names
    structure(res, class = "allFit", fit = object, sessionInfo =  sessionInfo(),
              data = data # is dropped if NULL
              )
}

print.allFit <- function(x, width=80, ...) {
    cat("original model:\n")
    f <- attr(x,"fit")
    ss <- function(x) {
        if (nchar(x)>width) {
            strwrap(paste0(substr(x,1,width-3),"..."))
        } else x
    }
    ff <- ss(deparse1(formula(f)))
    cat(ff,"\n")
    cat("data: ",deparse(getCall(f)$data),"\n")
    cat("optimizers (",length(x),"): ",
        ss(paste(names(x),collapse=", ")),"\n",
        sep="")
    which.bad <- vapply(x,FUN=is,"error",FUN.VALUE=logical(1))
    if ((nbad <- sum(which.bad))>0) {
        cat(nbad,"optimizer(s) failed\n")
    }
    cat("differences in negative log-likelihoods:\n")
    nllvec <- -vapply(x[!which.bad],logLik,numeric(1))
    cat("max=",signif(max(nllvec-min(nllvec)),3),
        "; std dev=",signif(sd(nllvec),3), "\n")
    ## FIXME: print magnitudes of parameter diffs
    ## cat("differences in parameters:\n")
    ## ss <- summary(x)
    ## allpars <- cbind(ss$fixef, ss$sdcor)
    ## par_max <-
    invisible(x)
}

summary.allFit <- function(object, ...) {
    afun <- function(x, FUN, ...) {
        f1 <- FUN(x[[1]], ...)
        nm <- names(f1)
        n <- length(f1)
        res <- vapply(x, FUN, numeric(n), ...)
        if (!is.null(dim(res))) {
            res <- t(res)
        } else {
            res <- as.matrix(res)
            colnames(res) <- nm
        }
        res
    }
    which.OK <- !vapply(object, is, "error", FUN.VALUE=logical(1))
    objOK <- object[which.OK]
    msgs <- lapply(objOK, function(x) x@optinfo$conv$lme4$messages)
    fixef <- afun(objOK, fixef)
    llik <- vapply(objOK, logLik, numeric(1))
    times <- afun(objOK, attr, "time")
    feval <- vapply(objOK, function(x) x@optinfo$feval, numeric(1))
    vfun <- function(x) as.data.frame(VarCorr(x))[["sdcor"]]
    sdcor <- afun(objOK, vfun)
    theta <- afun(objOK, getME, name="theta")
    cnm <- tnames(objOK[[1]])
    if (sigma(objOK[[1]])!=1) cnm <- c(cnm,"sigma")
    colnames(sdcor) <- unname(cnm)
    sdcor <- as.data.frame(sdcor)
    res <- namedList(which.OK, msgs, fixef, llik, sdcor, theta, times, feval)
    class(res) <- "summary.allFit"
    res
}

## should add a print method for summary:
##  * fixed effects, random effects: summary of differences?

## not yet ...
plot.allFit <- function(x, abbr=16, ...) {
    values <- opt <- NULL ## R code check/non-standard evaluation
     if (! (requireNamespace("ggplot2"))) {
         stop("ggplot2 package must be installed to plot allFit objects")
     }
     aes <- NULL ## code check
     ss <- summary(x)
     ff <- stack(as.data.frame(ss$fixef))
     ff$opt <- rep(rownames(ss$fixef),length.out=nrow(ff))
     if (!is.null(abbr)) ff$opt <- abbreviate(ff$opt, minlength=abbr)
     (ggplot2::ggplot(ff, aes(values, opt, colour=opt))
         + ggplot2::geom_point()
         + ggplot2::facet_wrap(~ind,scale="free")
         + ggplot2::theme(legend.position="none")
     )
}

# Plot the results from the fixed effects produced by different optimizers. This function
# takes the output from lme4::allFit(), tidies it, selects fixed effects and plots them.
plot.fixef.allFit <- function(allFit_output,
                              # Set the same Y axis limits in every plot
                              shared_y_axis_limits = TRUE,
                              # Multiply Y axis limits by a factor (only
                              # available if shared_y_axis_limits = TRUE)
                              multiply_y_axis_limits = 1,
                              # Number of decimal points
                              decimal_points = NULL,
                              # Select predictors
                              select_predictors = NULL,
                              # Number of rows
                              nrow = NULL,
                              # Y axis title
                              y_title = 'Fixed effect',
                              # Alignment of the Y axis title
                              y_title_hjust = NULL,
                              # Add number to the names of optimizers
                              number_optimizers = TRUE,
                              # Plot height (useful for many predictors)
                              height = NULL) {

    # Check for required packages
    pkgs <- c("ggplot2", "patchwork")
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(missing_pkgs) > 0) {
        stop("The following packages are required: ",
             paste(missing_pkgs, collapse = ", "),
             ". Please install them.", call. = FALSE)
    }

    # Extract fixed effects from the allFit() output
    allFit_fixef_mat <- summary(allFit_output)$fixef
    allFit_fixef <- data.frame(
        Optimizer = rownames(allFit_fixef_mat),
        allFit_fixef_mat,
        check.names = FALSE
    )
    allFit_fixef <- stats::reshape(allFit_fixef,
                                   varying = colnames(allFit_fixef_mat),
                                   v.names = "value",
                                   timevar = "fixed_effect",
                                   times = colnames(allFit_fixef_mat),
                                   direction = "long",
                                   idvar = "Optimizer",
                                   ids = allFit_fixef$Optimizer)
    rownames(allFit_fixef) <- NULL # remove row names

    # If number_optimizers, assign number to each optimizer and place it before its name
    if (number_optimizers) {
        allFit_fixef$Optimizer <- paste0(as.numeric(factor(allFit_fixef$Optimizer)), '. ', allFit_fixef$Optimizer)
    }

    # If select_predictors were specified, select them along with the intercept
    if (!is.null(select_predictors)) {
        allFit_fixef <- allFit_fixef[allFit_fixef$fixed_effect %in% c('(Intercept)', select_predictors), ]
    }

    # Order variables
    allFit_fixef <- allFit_fixef[, c('Optimizer', 'fixed_effect', 'value')]

    # Warn if multiply_y_axis_limits is set without shared_y_axis_limits
    if (multiply_y_axis_limits != 1 && !shared_y_axis_limits) {
        message('The argument `multiply_y_axis_limits` has not been used because it requires `shared_y_axis_limits` set to TRUE.')
    }

    # Warn for extreme y_title_hjust values
    if (!is.null(y_title_hjust) && (y_title_hjust < 0.5 || y_title_hjust > 6)) {
        message('NOTE: For y_title_hjust, a working range of values is between 0.6 and 6.')
    }

    # Convert decimal_points to accuracy for scales::label_number
    accuracy <- if (!is.null(decimal_points)) 10^(-decimal_points) else NULL

    # Define common theme elements
    strip_style <- ggplot2::element_text(size = 10, margin = ggplot2::margin(t = 4, b = 6))
    strip_background_style <- ggplot2::element_rect(fill = 'grey96')
    blank_x_axis <- list(
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
    )

    # First row: intercept_plot
    intercept <- allFit_fixef[allFit_fixef$fixed_effect == '(Intercept)', ]
    intercept_plot <- ggplot2::ggplot(intercept, ggplot2::aes(fixed_effect, value, colour = Optimizer)) +
        ggplot2::geom_point(position = ggplot2::position_dodge(1)) +
        ggplot2::facet_wrap(~fixed_effect, scale = 'free') +
        ggplot2::guides(colour = ggplot2::guide_legend(title.position = 'left')) +
        blank_x_axis +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                       strip.text = strip_style,
                       strip.background = strip_background_style,
                       legend.margin = ggplot2::margin(0.3, 0, 0.8, 1, 'cm'),
                       legend.title = ggplot2::element_text(size = ggplot2::unit(15, 'pt'), angle = 90, hjust = 0.5))

    # Second row: predictors_plot
    predictors <- allFit_fixef[allFit_fixef$fixed_effect != '(Intercept)', ]
    predictors$fixed_effect <- factor(predictors$fixed_effect, levels = unique(predictors$fixed_effect))

    # Set number of rows for predictors
    n_predictors <- length(unique(predictors$fixed_effect))
    predictors_plot_nrow <- if (!is.null(nrow)) nrow - 1 else ceiling(n_predictors / 2)

    predictors_plot <- ggplot2::ggplot(predictors, ggplot2::aes(fixed_effect, value, colour = Optimizer)) +
        ggplot2::geom_point(position = ggplot2::position_dodge(1)) +
        ggplot2::facet_wrap(~fixed_effect, scale = 'free', nrow = predictors_plot_nrow) +
        ggplot2::labs(y = y_title) +
        blank_x_axis +
        ggplot2::theme(axis.title.y = ggplot2::element_text(size = 14, margin = ggplot2::margin(0, 15, 0, 5, 'pt')),
                       strip.text = strip_style,
                       strip.background = strip_background_style,
                       legend.position = 'none')

    # Apply y-axis scaling
    y_scale <- NULL
    if (shared_y_axis_limits) {
        y_range <- range(allFit_fixef$value)
        y_abs_max <- max(abs(allFit_fixef$value))
        y_limits <- c(y_range[1] - y_abs_max / 7 * multiply_y_axis_limits,
                      y_range[2] + y_abs_max / 7 * multiply_y_axis_limits)
        y_scale <- ggplot2::scale_y_continuous(limits = y_limits, labels = if (!is.null(accuracy)) scales::label_number(accuracy = accuracy) else ggplot2::waiver())
    } else if (!is.null(accuracy)) {
        y_scale <- ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = accuracy))
    }

    if (!is.null(y_scale)) {
        intercept_plot <- intercept_plot + y_scale
        predictors_plot <- predictors_plot + y_scale
    }

    # Adjust y-axis title justification
    if (is.null(y_title_hjust)) {
        # Automatic adjustment based on number of rows
        y_title_hjust <- 5 / (predictors_plot_nrow + 1.5)
    }
    predictors_plot <- predictors_plot + ggplot2::theme(axis.title.y = ggplot2::element_text(hjust = y_title_hjust))

    if (predictors_plot_nrow > 5) {
        message('Many rows! Consider distributing predictors into several plots using argument `select_predictors`')
    }

    # Define layout based on number of rows
    intercept_height <- 5.3
    row_height <- if (!is.null(height)) height / (predictors_plot_nrow + 1.2) else 5
    predictors_height <- intercept_height + 0.5 + predictors_plot_nrow * row_height
    layout <- c(
      patchwork::area(t = 1.5, r = 8.9, b = 1.5 + intercept_height, l = 0),  # intercept row
      patchwork::area(t = 1.5 + intercept_height + 0.5, r = 9, b = predictors_height, l = 0) # predictors row(s)
    )

    # Return matrix of plots
    patchwork::wrap_plots(intercept_plot, predictors_plot, design = layout, nrow = 2)
}

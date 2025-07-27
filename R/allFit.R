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
##' @param start_from_mle (logical) initialize the refitting process with starting values from the fitted model?
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
                   catch.errs = TRUE,
                   start_from_mle = TRUE) {

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
            # Using the MLE as a starting point
            if (start_from_mle) {
              if (isGLMM(object)) {
                pars <- getME(object, c("theta", "fixef"))
              } else {
                pars <- getME(object, "theta")
                if(isNLMM(object)){
                  warning("results are not guaranteed when using nlmer")
                }
              }
              assign("pars", pars, environment(formula(object)))
            }
            
            fit_call <- if (start_from_mle) {
              quote(update(object, start = pars, control = ctrl))
            } else {
              quote(update(object, control = ctrl))
            }
            
            tt <- system.time(
              rr <- if (catch.errs) {
                      tryCatch(eval(fit_call), error = function(e) e)
                    } else {
                      eval(fit_call)
                    }
            )

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
                              # Select predictors
                              select_predictors = NULL,
                              # Number of rows
                              nrow = NULL,
                              # X axis title
                              x_title = "Estimate",
                              # Y axis title
                              y_title = "Optimizer",
                              # Show x-axis title
                              show_x_title = TRUE,
                              # Show y-axis title
                              show_y_title = TRUE,
                              # Show x-axis title on inner plots
                              inner_x_title = FALSE,
                              # Show y-axis title on inner plots
                              inner_y_title = FALSE,
                              # Show intercept
                              show_intercept = TRUE,
                              # Add point ranges (confidence intervals)
                              point_ranges = TRUE,
                              # Confidence level for point ranges
                              conf.level = 0.95,
                              # Decimal points for rounding
                              decimal_points = NULL,
                              # Shared x-axis limits across subplots
                              shared_x_axis_limits = FALSE,
                              # Additional arguments passed to tinyplot
                              ...) {
  
    # Check for required packages
    pkgs <- c("tinyplot")
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(missing_pkgs) > 0) {
        stop("The following packages are required: ",
             paste(missing_pkgs, collapse = ", "),
             ". Please install them.", call. = FALSE)
    }

    # Extract fixed effects from successful fits
    ss <- summary(allFit_output)
    successful_fits <- allFit_output[ss$which.OK]
    
    if (length(successful_fits) == 0) {
        stop("No successful fits found in allFit output", call. = FALSE)
    }
    
    # Extract fixed effects and standard errors using approach similar to developer's draft
    get_fe <- function(x, opt) {
        cc <- coef(summary(x))
        dd <- data.frame(
            optimizer = opt, 
            term = rownames(cc), 
            estimate = cc[,"Estimate"], 
            sd = cc[,"Std. Error"],
            stringsAsFactors = FALSE
        )
        
        # Add confidence intervals if requested
        if (point_ranges && !is.null(conf.level) && !is.na(conf.level)) {
            qq <- qnorm((1 + conf.level) / 2)
            dd <- transform(dd,
                           lwr = estimate - qq * sd,
                           upr = estimate + qq * sd)
            
            # Filter out extreme confidence intervals that might be due to convergence issues
            # If CI width is more than 100 times the estimate, treat as suspect
            ci_width <- dd$upr - dd$lwr
            estimate_abs <- abs(dd$estimate)
            extreme_ci <- ci_width > 100 * pmax(estimate_abs, 0.01)  # Use minimum threshold
            
            if (any(extreme_ci)) {
                dd$lwr[extreme_ci] <- NA
                dd$upr[extreme_ci] <- NA
            }
        }
        rownames(dd) <- NULL
        dd
    }
    
    # Get data for all successful fits
    fe_list <- Map(get_fe, successful_fits, names(successful_fits))
    allFit_fixef <- do.call(rbind, fe_list)
    
    # Apply rounding if specified
    if (!is.null(decimal_points)) {
        allFit_fixef$estimate <- round(allFit_fixef$estimate, decimal_points)
        if ("lwr" %in% colnames(allFit_fixef)) {
            allFit_fixef$lwr <- round(allFit_fixef$lwr, decimal_points)
            allFit_fixef$upr <- round(allFit_fixef$upr, decimal_points)
        }
    }
    
    # Filter predictors if specified
    if (!is.null(select_predictors)) {
        allFit_fixef <- allFit_fixef[allFit_fixef$term %in% select_predictors, ]
    }

    # Remove intercept if requested
    if (!show_intercept) {
        allFit_fixef <- allFit_fixef[allFit_fixef$term != "(Intercept)", ]
    }
    
    # Calculate layout
    n_predictors <- length(unique(allFit_fixef$term))
    if (n_predictors == 0) {
        stop("No predictors left to plot. Check select_predictors or show_intercept.", call. = FALSE)
    }
    plot_nrow <- if (!is.null(nrow)) nrow else ceiling(sqrt(n_predictors))
    plot_ncol <- ceiling(n_predictors / plot_nrow)
    
    # Calculate shared x-axis limits if requested
    x_limits <- NULL
    if (shared_x_axis_limits) {
        if (point_ranges && "lwr" %in% colnames(allFit_fixef)) {
            # Use finite values only for axis limits
            finite_lwr <- allFit_fixef$lwr[is.finite(allFit_fixef$lwr)]
            finite_upr <- allFit_fixef$upr[is.finite(allFit_fixef$upr)]
            if (length(finite_lwr) > 0 && length(finite_upr) > 0) {
                x_limits <- range(c(finite_lwr, finite_upr), na.rm = TRUE)
            } else {
                x_limits <- range(allFit_fixef$estimate, na.rm = TRUE)
            }
        } else {
            x_limits <- range(allFit_fixef$estimate, na.rm = TRUE)
        }
        # Add some padding
        x_range <- diff(x_limits)
        x_limits <- x_limits + c(-0.05 * x_range, 0.05 * x_range)
    }
    
    # Set up plotting layout with adjusted margins and spacing
    old_par <- par(mfrow = c(plot_nrow, plot_ncol), 
                 mar = c(3, 12, 2, 2),    # Increased right margin for better column spacing
                 oma = c(0, 0, 0, 0),
                 mgp = c(2, 0.7, 0))      # Adjust axis label positioning
    on.exit(par(old_par))
  
    # Get additional arguments from ...
    extra_args <- list(...)
    
    # Plot each fixed effect in a separate panel
    unique_effects <- unique(allFit_fixef$term)
    for (i in seq_along(unique_effects)) {
        effect <- unique_effects[i]
        effect_data <- allFit_fixef[allFit_fixef$term == effect, ]
        
        # Calculate position in grid
        current_row <- ceiling(i / plot_ncol)
        current_col <- ((i - 1) %% plot_ncol) + 1
        
        # Determine if this is the bottom-most plot in its column
        plots_in_this_col <- seq(current_col, n_predictors, by = plot_ncol)
        is_bottom_in_col <- i == max(plots_in_this_col)
        
        # Determine if this plot should show axis titles
        show_x_for_this_plot <- show_x_title && (inner_x_title || is_bottom_in_col)
        show_y_for_this_plot <- show_y_title && (inner_y_title || current_col == 1)
        
        # Create y-axis positions for optimizers
        y_pos <- seq_len(nrow(effect_data))
        
        # Prepare default plotting arguments
        default_args <- list(
          x = effect_data$estimate,
          y = y_pos,
          pch = 16,
          xlab = if (show_x_for_this_plot) x_title else "",
          ylab = "",  # Remove ylab from tinyplot, handle separately with mtext
          main = effect,
          yaxt = "n"
        )
        
        # Add x-axis limits if shared
        if (!is.null(x_limits)) {
            default_args$xlim <- x_limits
        } else if (point_ranges && "lwr" %in% colnames(effect_data)) {
            # Calculate individual plot limits based on confidence intervals
            # Use a more robust approach that excludes extreme outliers
            finite_lwr <- effect_data$lwr[is.finite(effect_data$lwr) & !is.na(effect_data$lwr)]
            finite_upr <- effect_data$upr[is.finite(effect_data$upr) & !is.na(effect_data$upr)]
            
            if (length(finite_lwr) > 0 && length(finite_upr) > 0) {
                # Use quantiles to exclude extreme outliers
                all_values <- c(finite_lwr, finite_upr, effect_data$estimate)
                all_values <- all_values[is.finite(all_values)]
                
                if (length(all_values) > 0) {
                    q01 <- quantile(all_values, 0.01, na.rm = TRUE)
                    q99 <- quantile(all_values, 0.99, na.rm = TRUE)
                    plot_xlim <- c(q01, q99)
                    plot_range <- diff(plot_xlim)
                    plot_xlim <- plot_xlim + c(-0.1 * plot_range, 0.1 * plot_range)
                    default_args$xlim <- plot_xlim
                }
            }
        }
        
        # Merge with user arguments, giving precedence to user arguments
        plot_args <- utils::modifyList(default_args, extra_args)
        
        # Create the plot
        tinyplot_fun <- get("tinyplot", envir = asNamespace("tinyplot"))
        do.call(tinyplot_fun, plot_args)
        
        # Add y-axis labels to all plots by default
        axis(2, at = y_pos, labels = effect_data$optimizer, las = 1, cex.axis = 0.7, tick = FALSE)
        
        # Add y-axis title with proper positioning (only for appropriate plots)
        if (show_y_for_this_plot && y_title != "") {
            mtext(y_title, side = 2, line = 11, cex = 0.8)
        }
        
        # Add point ranges if available and requested
        if (point_ranges && "lwr" %in% colnames(effect_data)) {
            for (j in seq_len(nrow(effect_data))) {
                if (!is.na(effect_data$lwr[j]) && !is.na(effect_data$upr[j]) && 
                    is.finite(effect_data$lwr[j]) && is.finite(effect_data$upr[j])) {
                    # Add horizontal line for confidence interval
                    segments(x0 = effect_data$lwr[j],
                            x1 = effect_data$upr[j],
                            y0 = y_pos[j], 
                            y1 = y_pos[j],
                            lwd = 1)
                    # Add small vertical caps at the ends
                    segments(x0 = effect_data$lwr[j], x1 = effect_data$lwr[j],
                            y0 = y_pos[j] - 0.1, y1 = y_pos[j] + 0.1)
                    segments(x0 = effect_data$upr[j], x1 = effect_data$upr[j],
                            y0 = y_pos[j] - 0.1, y1 = y_pos[j] + 0.1)
                }
            }
        }
    }
    
    if (plot_nrow > 3) {
        message("Many rows! Consider distributing predictors into several plots using argument `select_predictors`")
    }
    
    invisible(allFit_fixef)
}

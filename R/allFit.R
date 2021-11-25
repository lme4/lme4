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

nmkbw <- function(fn,par,lower,upper,control) {
    if (length(par)==1) {
        res <- optim(fn=fn,par=par,lower=lower,upper=100*par,
                     method="Brent")
    } else {
        if (!is.null(control$maxfun)) {
            control$maxfeval <- control$maxfun
            control$maxfun <- NULL
        }
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
                   maxfun=1e5,
                   parallel = c("no", "multicore", "snow"),
                   ncpus = getOption("allFit.ncpus", 1L),
                   cl = NULL,
                   catch.errs = TRUE) {

    if (is.null(meth.tab)) {
        meth.tab <- meth.tab.0
    }
    if (!requireNamespace("dfoptim")) {
        meth.tab <- meth.tab[meth.tab.0[,"optimizer"] != "nmkbw",]
    }
    if (!requireNamespace("optimx")) {
        meth.tab <- meth.tab[meth.tab.0[,"optimizer"] != "optimx",]
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
                if (is.null(cc)) return(x)  ## easy! no controls specified
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
                                   mkOptCtrl("maxfun"=maxfun))
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

plot.fixef.allFit = function(allFit_output,
                             ## Set the same Y axis limits in every plot
                             shared_y_axis_limits = TRUE,
                             ## Multiply Y axis limits by a factor (only
                             ## available if shared_y_axis_limits = TRUE)
                             multiply_y_axis_limits = 1,
                             ## Number of decimal points
                             decimal_points = NULL,
                             ## Select predictors
                             select_predictors = NULL,
                             ## Number of rows
                             nrow = NULL,
                             ## Y axis title
                             y_title = 'Fixed effect',
                             ## Alignment of the Y axis title
                             y_title_hjust = NULL,
                             ## Add number to the names of optimizers
                             number_optimizers = TRUE,
                             ## Replace colon in interactions with x
                             interaction_symbol_x = TRUE) {

  extra_pkgs <- c('dplyr', 'reshape2', 'stringr', 'scales', 'ggplot2')

  # Data wrangling
  if (!requireNamespace('dplyr')) install.packages('dplyr')
  if (!requireNamespace('reshape2')) install.packages('reshape2')
  # Text processing
  if (!requireNamespace('stringr')) install.packages('stringr')
  # Set number of decimal points
  if (!requireNamespace('scales')) install.packages('scales')
  # Plotting
  if (!requireNamespace('ggplot2')) install.packages('ggplot2')
  # Matrix of plots
  if (!requireNamespace('patchwork')) install.packages('patchwork')

  require(dplyr)
  require(reshape2)
  require(stringr)
  require(scales)
  require(ggplot2)
  require(patchwork)

  # Tidy allFit output

  # Extract fixed effects from the allFit() output
  allFit_fixef = summary(allFit_output)$fixef %>%  # Select fixed effects in the allFit results
    reshape2::melt() %>%  # Structure the output as a data frame
    rename('Optimizer' = 'Var1', 'fixed_effect' = 'Var2')  # set informative names

  # If number_optimizers = TRUE, assign number to each optimizer and place it before its name
  if(number_optimizers == TRUE) {
    allFit_fixef$Optimizer = paste0(as.numeric(allFit_fixef$Optimizer), '. ', allFit_fixef$Optimizer)
  }

  # If select_predictors were specified, select them along with the intercept (the latter required)
  if(!is.null(select_predictors)) {
    allFit_fixef = allFit_fixef %>% filter(fixed_effect %in% c('(Intercept)', select_predictors))
  }

  # Order variables
  allFit_fixef = allFit_fixef[, c('Optimizer', 'fixed_effect', 'value')]

  # PLOT. The overall plot is formed of a first row containing the intercept and the legend
  # (intercept_plot), and a second row containing the predictors (predictors_plot),
  # which may in turn occupy several rows.

  # If multiply_y_axis_limits was specified but shared_y_axis_limits = FALSE,
  # warn that shared_y_axis_limits is required.
  if(!multiply_y_axis_limits == 1 & shared_y_axis_limits == FALSE) {
    message('The argument `multiply_y_axis_limits` has not been used because it requires `shared_y_axis_limits` set to TRUE.')
  }

  # If extreme values were entered in y_title_hjust, show warning
  if(!is.null(y_title_hjust)) {
    if(y_title_hjust < 0.5 | y_title_hjust > 6) {
      message('NOTE: For y_title_hjust, a working range of values is between 0.6 and 6.')
    }
  }

  # If decimal_points were specified, convert number to the format used in 'scales' package
  if(!is.null(decimal_points)) {
    decimal_points =
      ifelse(decimal_points == 1, 0.1,
             ifelse(decimal_points == 2, 0.01,
                    ifelse(decimal_points == 3, 0.001,
                           ifelse(decimal_points == 4, 0.0001,
                                  ifelse(decimal_points == 5, 0.00001,
                                         ifelse(decimal_points == 6, 0.000001,
                                                ifelse(decimal_points == 7, 0.0000001,
                                                       ifelse(decimal_points == 8, 0.00000001,
                                                              ifelse(decimal_points == 9, 0.000000001,
                                                                     ifelse(decimal_points == 10, 0.0000000001,
                                                                            ifelse(decimal_points == 11, 0.00000000001,
                                                                                   ifelse(decimal_points == 12, 0.000000000001,
                                                                                          ifelse(decimal_points == 13, 0.0000000000001,
                                                                                                 ifelse(decimal_points == 14, 0.00000000000001,
                                                                                                        ifelse(decimal_points > 15, 0.000000000000001,
                                                                                                               0.001
                                                                                                        )))))))))))))))
  }

  # First row: intercept_plot

  # Select intercept data only
  intercept = allFit_fixef %>% filter(fixed_effect == '(Intercept)')

  intercept_plot = intercept %>%
    ggplot(., aes(fixed_effect, value, colour = Optimizer)) +
    geom_point(position = position_dodge(1)) +
    facet_wrap(~fixed_effect, scale = 'free') +
    guides(colour = guide_legend(title.position = 'left')) +
    theme_bw() + theme(axis.title = element_blank(), axis.ticks.x = element_blank(),
                       axis.text.x = element_blank(),
                       strip.text = element_text(size = 10, margin = margin(t = 4, b = 6)),
                       strip.background = element_rect(fill = 'grey96'),
                       legend.margin = margin(0.3, 0, 0.8, 1, 'cm'),
                       legend.title = element_text(size = unit(15, 'pt'), angle = 90, hjust = 0.5))

  # Second row: predictors_plot

  # Select all predictors except intercept
  predictors = allFit_fixef %>% filter(!fixed_effect == '(Intercept)')

  # If interaction_symbol_x = TRUE (default), replace colon with times symbol x between spaces
  if(interaction_symbol_x == TRUE) {
    # Replace colon in interactions with \u00D7, i.e., x; then set factor class
    predictors$fixed_effect = predictors$fixed_effect %>% str_replace_all(':', ' \u00D7 ') %>% factor()
  }

  # Order predictors as in the original output from lme4::allFit()
  predictors$fixed_effect = factor(predictors$fixed_effect, levels = unique(predictors$fixed_effect))

  # Set number of rows for the predictors excluding the intercept.
  # First, if nrow argument specified, use it
  if(!is.null(nrow)) {
    predictors_plot_nrow = nrow - 1  # Subtract 1 as intercept row not considered

    # Else, if nrow argument not specified, calculate sensible number of rows: i.e., divide number of
    # predictors (exc. intercept) by 2 and round up the result. For instance, 7 predictors --> 3 rows
  } else predictors_plot_nrow = (length(unique(predictors$fixed_effect)) / 2) %>% ceiling()

  predictors_plot = ggplot(predictors, aes(fixed_effect, value, colour = Optimizer)) +
    geom_point(position = position_dodge(1)) +
    facet_wrap(~fixed_effect, scale = 'free',
               # Note that predictors_plot_nrow was defined a few lines above
               nrow = predictors_plot_nrow) +
    labs(y = y_title) +
    theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.y = element_text(size = 14, margin = margin(0, 15, 0, 5, 'pt')),
                       strip.text = element_text(size = 10, margin = margin(t = 4, b = 6)),
                       strip.background = element_rect(fill = 'grey96'), legend.position = 'none')

  # Below, the function scale_y_continuous is applied conditionally to avoid overriding settings. First,
  # if shared_y_axis_limits = TRUE and decimal_points were specified, set the same Y axis limits in
  # every plot and set decimal_points. By default, also expand limits by a seventh of its original limit,
  # and allow further multiplication of limits through multiply_y_axis_limits.
  if(shared_y_axis_limits == TRUE & !is.null(decimal_points)) {

    intercept_plot = intercept_plot +
      scale_y_continuous(limits = c(min(allFit_fixef$value) - allFit_fixef$value %>% abs %>% max / 7 * multiply_y_axis_limits,
                                    max(allFit_fixef$value) + allFit_fixef$value %>% abs %>% max / 7 * multiply_y_axis_limits),
                         # Set number of decimal points
                         labels = scales::label_number(accuracy = decimal_points))

    predictors_plot = predictors_plot +
      scale_y_continuous(limits = c(min(allFit_fixef$value) - allFit_fixef$value %>% abs %>% max / 7 * multiply_y_axis_limits,
                                    max(allFit_fixef$value) + allFit_fixef$value %>% abs %>% max / 7 * multiply_y_axis_limits),
                         # Set number of decimal points
                         labels = scales::label_number(accuracy = decimal_points))

    # Else, if shared_y_axis_limits = TRUE but decimal_points were not specified, do as above but without
    # setting decimal_points.
  } else if(shared_y_axis_limits == TRUE & is.null(decimal_points)) {

    intercept_plot = intercept_plot +
      scale_y_continuous(limits = c(min(allFit_fixef$value) - allFit_fixef$value %>% abs %>% max / 7 * multiply_y_axis_limits,
                                    max(allFit_fixef$value) + allFit_fixef$value %>% abs %>% max / 7 * multiply_y_axis_limits))

    predictors_plot = predictors_plot +
      scale_y_continuous(limits = c(min(allFit_fixef$value) - allFit_fixef$value %>% abs %>% max / 7 * multiply_y_axis_limits,
                                    max(allFit_fixef$value) + allFit_fixef$value %>% abs %>% max / 7 * multiply_y_axis_limits))

    # Else, if shared_y_axis_limits = FALSE and decimal_points were specified, set decimal_points.
  } else if(shared_y_axis_limits == FALSE & !is.null(decimal_points)) {

    intercept_plot = intercept_plot +
      scale_y_continuous(labels = scales::label_number(accuracy = decimal_points))

    predictors_plot = predictors_plot +
      scale_y_continuous(labels = scales::label_number(accuracy = decimal_points))
  }

  # Plot matrix: based on number of predictors_plot_nrow, adjust height of Y axis title
  # (unless specified by user), and assign space to intercept_plot and predictors_plot
  if(predictors_plot_nrow == 1) {

    # If y_title_hjust specified by user, use it
    if(!is.null(y_title_hjust)) {
      predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = y_title_hjust))
      # Otherwise, set a sensible height
    } else predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = 3.6))

    layout = c(
      area(t = 1.5, r = 8.9, b = 6.8, l = 0),  # intercept row
      area(t = 7.3, r = 9, b = 11, l = 0)      # predictors row(s)
    )

  } else if(predictors_plot_nrow == 2) {

    # If y_title_hjust specified by user, use it
    if(!is.null(y_title_hjust)) {
      predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = y_title_hjust))
      # Otherwise, set a sensible height
    } else predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = 1.4))

    layout = c(
      area(t = 1.5, r = 8.9, b = 6.8, l = 0),  # intercept row
      area(t = 7.3, r = 9, b = 16, l = 0)      # predictors row(s)
    )

  } else if(predictors_plot_nrow == 3) {

    # If y_title_hjust specified by user, use it
    if(!is.null(y_title_hjust)) {
      predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = y_title_hjust))
      # Otherwise, set a sensible height
    } else predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = 0.92))

    layout = c(
      area(t = 1.5, r = 8.9, b = 6.8, l = 0),  # intercept row
      area(t = 7.3, r = 9, b = 21, l = 0)      # predictors row(s)
    )

  } else if(predictors_plot_nrow == 4) {

    # If y_title_hjust specified by user, use it
    if(!is.null(y_title_hjust)) {
      predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = y_title_hjust))
      # Otherwise, set a sensible height
    } else predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = 0.8))

    layout = c(
      area(t = 1.5, r = 8.9, b = 6.8, l = 0),  # intercept row
      area(t = 7.3, r = 9, b = 26, l = 0)      # predictors row(s)
    )

  } else if(predictors_plot_nrow == 5) {

    # If y_title_hjust specified by user, use it
    if(!is.null(y_title_hjust)) {
      predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = y_title_hjust))
      # Otherwise, set a sensible height
    } else predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = 0.73))

    layout = c(
      area(t = 1.5, r = 8.9, b = 6.8, l = 0),  # intercept row
      area(t = 7.3, r = 9, b = 31, l = 0)      # predictors row(s)
    )

  } else if(predictors_plot_nrow > 5) {

    # If y_title_hjust specified by user, use it
    if(!is.null(y_title_hjust)) {
      predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = y_title_hjust))
      # Otherwise, set a sensible height
    } else predictors_plot = predictors_plot + theme(axis.title.y = element_text(hjust = 0.65))

    layout = c(
      area(t = 1.5, r = 8.9, b = 6.8, l = 0),  # intercept row
      area(t = 7.3, r = 9, b = 36, l = 0)      # predictors row(s)
    )

    # Also, advise user to consider distributing predictors into several plots
    message('Many rows! Consider distributing predictors into several plots using argument `select_predictors`')
  }

  # Return matrix of plots
  wrap_plots(intercept_plot, predictors_plot, design = layout,
             # The 2 below corresponds to intercept_plot and predictors_plot
             nrow = 2)

}

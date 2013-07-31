\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{lme4 News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.0 (2013-08-01)}{
  \subsection{MAJOR USER-VISIBLE CHANGES}{
    \itemize{
      \item Because the internal computational machinery has changed, 
      results from the newest version of `lme4` will not be numerically 
      identical to those from previous versions.  For reasonably well-
      defined fits, they will be extremely close (within numerical 
      tolerances of 1e-4 or so), but for unstable or poorly-defined fits 
      the results may change, and very unstable fits may fail when they 
      (apparently) succeeded with previous versions. Similarly, some fits 
      may be slower with the new version, although on average the new 
      version should be faster and more stable. There are more numerical 
      tuning options available than before (see below); non-default settings 
      may restore the speed and/or ability to fit a particular model without 
      an error. If you notice significant or disturbing changes when fitting 
      a model with the new version of `lme4`, **please notify the maintainers**.
      \item In the past, `lmer` automatically called `glmer` when `family` was
      specified. It still does so, but now warns the user that they should
      preferably use `glmer` directly.
      \item `VarCorr` returns its results in the same format as before (as a 
      list of variance-covariance matrices with `correlation` and `stddev` 
      attributes, plus a `sc` attribute giving the residual standard 
      deviation/scale parameter when appropriate), but prints them in a 
      different (nicer) way.
      \item By default `residuals` gives deviance (rather than Pearson) 
      residuals when applied to `glmer` fits (a side effect of matching `glm`
      behaviour more closely).
    }
  }
  \subsection{MINOR USER-VISIBLE CHANGES}{
    \itemize{
      \item More use is made of S3 rather than S4 classes and methods: one 
      side effect is that methods such as `fixef` no longer conflict with 
      the `nlme` package.
      \item The internal optimizer has changed. `[gn]lmer` now has an 
      `optimizer` argument; `"Nelder_Mead"` is the default for `[n]lmer`, 
      while a combination of `"bobyqa"` (an alternative derivative-free 
      method) and `"Nelder_Mead"` is the default for `glmer`; to use the 
      `nlminb` optimizer as in the old version of `lme4`, you can use 
      `optimizer="optimx"` with `control=list(method="nlminb")` (you will 
      need the `optimx` package to be installed and loaded). See the help 
      pages for details.
      \item Families in GLMMs are no longer restricted to built-in/hard-
      coded families; any family described in `?family`, or following that 
      design, is usable (although there are some hard-coded families, which 
      will be faster)
      \item `[gn]lmer` now produces objects of class `merMod` rather than 
      class `mer` as before
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item A general-purpose `getME()` accessor method has been used to 
      allow extraction of a wide variety of components of a mixed-model 
      fit. `getME()` also allows a vector of objects to be returned as
      a list of mixed-model componenets. This has been backported to
      be compatible with older versions of `lme4` that still produce `mer` 
      objects rather than `merMod` objects. However, backporting is incomplete 
      and some objects are still not extractable using `getME()` in older versions
      of `lme4`.
      \item Optimization information (convergence info, warnings, etc.) now
      stored in an `@optinfo` slot.
      \item `bootMer`, a framework for obtaining parameter confidence 
      intervals by parametric bootstrapping.
      \item `plot` methods similar to those from the `nlme` package 
      (although missing `augPred`).
      \item A `predict` method, allowing a choice of which random effects are 
      included in the prediction.
      \item Likelihood profiling (and profile confidence intervals) for 
      `lmer` and `glmer` results.
      \item `nAGQ=0`, an option to do fast (but inaccurate) fitting of GLMMs.
      \item Using `devFunOnly=TRUE` allows the user to extract a deviance 
      function for the model, allowing further diagnostics/customization of 
      model results.
      \item The internal structure of [gn]lmer is now more modular, allowing 
      finer control of the different steps of argument checking; construction 
      of design matrices and data structures; parameter estimation; and construction 
      of the final `merMod` object (see `?modular`).
    }
  }
  \subsection{STILL NON-EXISTENT FEATURES}{
    \itemize{
      \item Automatic MCMC sampling based on the fit turns out to be very difficult
      to implement in a way that is really broadly reliable and robust; `mcmcsamp`
      will not be implemented in the near future. We recommend parametric boostrapping 
      via `bootMer`, or the Kenward-Roger approximation implemented in the `pbkrtest`
      package and leveraged by the `lmerTest` package and the `Anova` function in the 
      `car` package.
      \item "R-side" structures (within-block correlation and heteroscedasticity) are 
      not on the current timetable
      \item Adaptive Gaussian quadrature (AGQ) is not available for multiple and/or
      non-scalar random effects.
      \item 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item In a development version, prior weights were not being used properly in 
      the calculation of the residual standard deviation, but this has been fixed.
      Thank you to Simon Wood for pointing this out.
      \item In a development version, the step-halving component of the penalized 
      iteratively reweighted least squares algorithm was not working, but
      this is now fixed.
      \item In a development version, square `RZX` matrices would lead to a 
      `pwrssUpdate did not converge in 30 iterations` error. This has been fixed
      by adding an extra column of zeros to `RZX`.
    }
  }
}
\section{CHANGES IN VERSION 0.999375-16 (2008-06-23)}{
  \subsection{MAJOR USER-VISIBLE CHANGES}{
    \itemize{
       \item The underlying algorithms and representations for all the
       mixed-effects models fit by this package have changed - for
       the better, we hope. The class "mer" is a common
       mixed-effects model representation for linear, generalized
       linear, nonlinear and generalized nonlinear mixed-effects
       models.
       \item ECME iterations are no longer used at all, nor are analytic
       gradients. Components named 'niterEM', 'EMverbose', or
       'gradient' can be included in the 'control' argument to
       lmer(), glmer() or nlmer() but have no effect.
       \item PQL iterations are no longer used in glmer() and nlmer().
       Only the Laplace approximation is currently available. AGQ,
       for certain classes of GLMMs or NLMMs, is being added.
       \item The 'method' argument to lmer(), glmer() or nlmer() is
       deprecated. Use the 'REML = FALSE' in lmer() to obtain ML
       estimates. Selection of AGQ in glmer() and nlmer() will be
       controlled by the argument 'nAGQ', when completed.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item The representation of mixed-effects models has been
      dramatically changed to allow for smooth evaluation of the
      objective as the variance-covariance matrices for the random
      effects approach singularity. Beta testers found this
      representation to be more robust and usually faster than
      previous versions of lme4.
      \item The mcmcsamp function uses a new sampling method for the
      variance-covariance parameters that allows recovery from
      singularity. The update is not based on a sample from the
      Wishart distribution. It uses a redundant parameter
      representation and a linear least squares update.
      \item CAUTION: Currently the results from mcmcsamp look peculiar and
      are probably incorrect. I hope it is just a matter of my
      omitting a scaling factor but I have seen patterns such as
      the parameter estimate for some variance-covariance parameters
      being the maximum value in the chain, which is highly
      unlikely.
      \item The 'verbose' argument to lmer(), glmer() and nlmer() can be
      used instead of 'control = list(msVerbose = TRUE)'.
    }
  }
}

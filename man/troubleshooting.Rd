\name{troubleshooting}
\alias{troubleshooting}
\title{Troubleshooting}
\description{
  This page attempts to summarize some of the common problems with
  fitting \code{[gn]lmer} models and how to troubleshoot them.
  \itemize{
    \item \code{failure to converge in (xxxx) evaluations} The
    optimizer hit its maximum limit of function evaluations.  To
    increase this, use the \code{optControl} argument of
    \code{[g]lmerControl} -- for \code{Nelder_Mead} and \code{bobyqa}
    the relevant parameter is \code{maxfun}; for \code{optim} and
    \code{optimx}-wrapped optimizers, including
    \code{nlminbwrap}, it's \code{maxit}; for
    \code{nloptwrap}, it's \code{maxeval}.
    \item \code{Model failed to converge with max|grad| ...} The scaled
    gradient at the fitted (RE)ML estimates
    is worryingly large. Try
    \itemize{
      \item refitting the parameters starting at the current
      estimates: getting consistent results (with no warning)
      suggests a false positive
      \item switching optimizers: getting consistent
      results suggests there is not really a problem; getting a similar
      log-likelihood with different parameter estimates suggests that
      the parameters are poorly determined (possibly the result of
      a misspecified or overfitted model)
      \item compute values of the deviance in the neighbourhood of the estimated
      parameters to double-check that \code{lme4} has really found a
      local optimum.
    }
    \item \code{Hessian is numerically singular: parameters are not
      uniquely determined}
    The Hessian (inverse curvature matrix) at the maximum likelihood or REML
    estimates has a very large eigenvalue, indicating that (within numerical
    tolerances) the surface is completely flat in some direction.
    The model may be misspecified, or extremely badly scaled (see
    \code{"Model is nearly unidentifiable"}).
    \item \code{Model is nearly unidentifiable ... Rescale variables?}
    The Hessian (inverse curvature matrix) at the maximum likelihood or REML
    estimates has a large eigenvalue, indicating that the surface is
    nearly flat in some direction. Consider centering and/or scaling
    continuous predictor variables.
    \item \code{Contrasts can be applied only to factors with 2 or more levels}
    One or more of the categorical predictors in the model has fewer than two
    levels. This may be due to user error when converting these predictors to 
    factors prior to modeling, or it may result from some factor levels being
    eliminated due to \code{NA}s in other predictors. Double-check the number
    of data points in each factor level to see which one is the culprit:
    \code{lapply(na.omit(df[,vars]), table)} (where \code{df} is the 
    \code{data.frame} and \code{vars} are the column names of your predictor
    variables).
  }
}


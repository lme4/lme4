\name{troubleshooting}
\alias{troubleshooting}
\title{Troubleshooting}
\description{
  This page attempts to summarize some of the common problems with
  fitting \code{[gn]lmer} models and how to troubleshoot them.
  \itemize{
    \item \code{Model failed to converge with max|grad| ...} The scaled
    gradient at the fitted (RE)ML estimates
    is worryingly large. Try switching optimizers to see if you
    get consistent results, or
    compute values of the deviance in the neighbourhood of the estimated
    parameters to double-check that \code{lme4} has really found a
    local optimum.
    \item \code{Model is nearly unidentifiable ... Rescale variables?}
    The Hessian (curvature matrix) at the maximum likelihood or REML
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


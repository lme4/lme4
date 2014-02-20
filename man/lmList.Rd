\name{lmList}
\alias{lmList}
\title{List of lm Objects with a Common Model}
\usage{
  lmList(formula, data, family, subset, weights, na.action,
    offset, pool, ...)
}
\arguments{
  \item{formula}{a linear formula object of the form
  \code{y ~ x1+...+xn | g}. In the formula object, \code{y}
  represents the response, \code{x1,...,xn} the covariates,
  and \code{g} the grouping factor specifying the
  partitioning of the data according to which different
  \code{lm} fits should be performed.}
  \item{family}{an optional family specification for a
  generalized linear model.}
  \item{pool}{logical scalar, should the variance estimate
  pool the residual sums of squares}
  \item{...}{additional, optional arguments to be passed to
  the model function or family evaluation.}
  \item{data}{an optional data frame containing the
  variables named in \code{formula}. By default the
  variables are taken from the environment from which
  \code{lmer} is called. See Details.}
  \item{subset}{an optional expression indicating the
  subset of the rows of \code{data} that should be used in
  the fit. This can be a logical vector, or a numeric
  vector indicating which observation numbers are to be
  included, or a character vector of the row names to be
  included.  All observations are included by default.}
  \item{weights}{an optional vector of \sQuote{prior
  weights} to be used in the fitting process.  Should be
  \code{NULL} or a numeric vector.}
  \item{na.action}{a function that indicates what should
  happen when the data contain \code{NA}s.  The default
  action (\code{na.omit}, inherited from the 'factory
  fresh' value of \code{getOption("na.action")}) strips any
  observations with any missing values in any variables.}
  \item{offset}{this can be used to specify an \emph{a
  priori} known component to be included in the linear
  predictor during fitting. This should be \code{NULL} or a
  numeric vector of length equal to the number of cases.
  One or more \code{\link{offset}} terms can be included in
  the formula instead or as well, and if more than one is
  specified their sum is used.  See
  \code{\link{model.offset}}.}
}
\details{
  \itemize{
    \item \code{data} should be a data frame (not, e.g.
    a \code{groupedData} object from the \pkg{nlme} package);
    use \code{\link{as.data.frame}} first to convert the data.
    \item While \code{data} is optional, the
    package authors \emph{strongly} recommend its use,
    especially when later applying methods such as
    \code{update} and \code{drop1} to the fitted model
    (\emph{such methods are not guaranteed to work properly
      if \code{data} is omitted}). If \code{data} is omitted,
    variables will be taken from the environment of
    \code{formula} (if specified as a formula) or from the
    parent frame (if specified as a character vector).
  }
}
\description{
  List of lm Objects with a Common Model
}


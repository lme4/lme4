\name{mkSimulateTemplate}
\title{Make templates suitable for guiding mixed model simulations}
\alias{mkParsTemplate}
\alias{mkDataTemplate}
\description{
  Make data and parameter templates suitable for guiding mixed model
  simulations, by specifying a model formula and other information
  (EXPERIMENTAL).  Most useful for simulating balanced designs and for 
  getting started on unbalanced simulations.
}
\usage{
mkParsTemplate(formula, data)
mkDataTemplate(formula, data, nGrps = 2, nPerGrp = 1, rfunc = NULL, ...)
}
\arguments{
  \item{formula}{A mixed model formula (see \code{\link{lmer}}).}
  \item{data}{A data frame containing the names in \code{formula}.}
  \item{nGrps}{Number of levels of a grouping factor.}
  \item{nPerGrp}{Number of observations per level.}
  \item{rfunc}{Function for generating covariate data (e.g. \code{\link{rnorm}}.}
  \item{...}{Additional parameters for \code{rfunc}.}
}
\seealso{
  These functions are designed to be used with \code{\link{simulate.merMod}}.
}

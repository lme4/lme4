\name{merMCMC-class}
\docType{class}
% Classes
\alias{merMCMC-class}
% Methods
\alias{coerce,merMCMC,matrix-method}
%% FIXME: Document these, including their unexpected arguments
\alias{qqmath,merMCMC-method}
\alias{densityplot,merMCMC-method}
\alias{xyplot,merMCMC-method}
\alias{as.matrix,merMCMC-method}
\alias{as.data.frame,merMCMC-method}
\alias{coerce,merMCMC,data.frame-method} 
\title{Mixed-model Markov chain Monte Carlo results}
\description{

  Objects of class \code{"merMCMC"} are Markov chain Monte Carlo samples
  from the distribution of the parameters of a fitted mixed-effects
  model.

} 
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("merMCMC", ...)}
  or, more commonly, via the \code{mer} method for the generic
  \code{\link{mcmcsamp}} function.
}
\section{Slots}{
  \describe{
    \item{\code{Gp}:}{a copy of the \code{Gp} slot of the original
      \code{\linkS4class{mer}} object}
    \item{\code{ST}:}{matrix of samples from the parameters determining
      the \code{ST} slot of the \code{\linkS4class{mer}} object}
    \item{\code{call}:}{Matched call from the original
      \code{\linkS4class{mer}} object}
    \item{\code{deviance}:}{vector of samples of the (ML) deviance}
    \item{\code{dims}:}{a copy of the \code{dims} slot of the original
      \code{\linkS4class{mer}} object}
    \item{\code{fixef}:}{matrix of samples of the fixed-effects parameters}
    \item{\code{nc}:}{Integer vector of length \code{dims["nf"]}.  The
      number of columns of random effects in each term.}
    \item{\code{ranef}:}{matrix of samples of the random effects.  This
      matrix has zero columns unless \code{saveb = TRUE} is specified in
      the call to \code{\link{mcmcsamp}}.  Consider the size of this
      matrix, which could be very large, before setting \code{saveb = TRUE}.}
    \item{\code{sigma}:}{vector of samples of the common scale
      parameter or \code{numeric(0)} if \code{dims["useSc"]} is \code{FALSE}.}
  }
}
\section{Methods}{
  \describe{
    \item{HPDinterval}{\code{signature(object = "merMCMC")}: use the
      chain to calculate Highest Posterior Density (HPD) intervals of a
      given empirical probability content for the model parameters.  See
      \code{\link{HPDinterval}}.}
    \item{VarCorr}{\code{signature(x = "merMCMC")}: transform the
      \code{ST} and \code{sigma} slots to some combination of variances, 
      covariances, standard deviations and correlations.  See
      \code{\link{VarCorr}} for details.}
    \item{as.data.frame}{\code{signature(x = "merMCMC")}: returns the
      fixef-effects and variance-covariance parameters from the chain in
      the form of a data frame.  The \code{type} argument for the
      \code{\link{VarCorr}} method can be passed to this method to
      select the type of variance-covariance parameters returned.}
    \item{as.matrix}{\code{signature(x = "merMCMC")}: Same as the
      \code{as.data.frame} method described above but returning a matrix.}
    \item{coerce}{\code{signature(from = "merMCMC", to = "data.frame")}:
      Same as the \code{as.data.frame} method.}
    \item{densityplot}{\code{signature(object = "merMCMC")}: plot
      empirical densities for the parameters from the chain.  See also
      \code{\link{densityplot}}.}
    \item{qqmath}{\code{signature(object = "merMCMC")}: plot
      quantile-quantile plots for the parameters from the sample in the
      chain. See also \code{\link{qqmath}}.}
    \item{xyplot}{\code{signature(object = "merMCMC")}: plot
      traces of the parameter samples in the chain.}
  }
}
\seealso{
  \code{\link{mcmcsamp}} produces these objects,
  \code{\link{lmer}}, \code{\link{glmer}} and \code{\link{nlmer}}
  produce the \code{\linkS4class{mer}} objects.
}
\examples{
showClass("merMCMC")
}
\keyword{classes}

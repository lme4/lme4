\name{VarCorr}
\docType{genericFunction}
% Generic
\alias{VarCorr}
% Methods
\alias{VarCorr,mer-method}
\alias{VarCorr,merMCMC-method}
\title{Extract variance and correlation components}
\usage{
\S4method{VarCorr}{mer}(x, \dots)
\S4method{VarCorr}{merMCMC}(x, type = c("raw", "varcov", "sdcorr", "logs"), \dots)
}
\arguments{
  \item{x}{a fitted model object, usually an object inheriting from
    class \code{\linkS4class{mer}}.
  }
  \item{type}{character string indicating the type of result to be
    returned, either \code{"raw"}, the raw representation as in the
    \code{ST} slot of the \code{\linkS4class{mer}} class, or
    \code{"varcov"}, variances and covariances, or \code{"sdcorr"},
    standard deviations and correlations, or \code{"logs"}, logs of the
    standard deviations and Fisher's z transformation of the correlations.}
  \item{\dots}{Additional, optional arguments for some methods.  At
    present none are used.}
}
\description{
  Extract the estimated variances, standard deviations, and correlations
  of the random-effects terms in a mixed-effects model, of class
  \code{mer}.

  When appropriate, the within-group error variance and standard
  deviation are also calculated.
}
\value{
  an object of class \code{VarCorr}.
}
\seealso{
  the \code{\link{lmer}} function and \code{\linkS4class{mer}} class;
  the result class \code{\linkS4class{VarCorr}}.
}
\examples{
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
             data = sleepstudy))
(VC <- VarCorr(fm2))
}
\keyword{models}

\name{mcmcsamp}
\docType{genericFunction}
\alias{mcmcsamp}
\alias{mcmcsamp,mer-method}
\title{Generate an MCMC sample}
\description{
  This generic function generates a sample from the posterior
  distribution of the parameters of a fitted model using Markov Chain
  Monte Carlo methods.
}
\usage{
mcmcsamp(object, n, verbose, \dots)
\S4method{mcmcsamp}{mer}(object, n, verbose, saveb, \dots)
}
\arguments{
  \item{object}{An object of a suitable class - usually an
    \code{"\linkS4class{mer}"} object.}
  \item{n}{integer - number of samples to generate. Defaults to 1.}
  \item{verbose}{logical - if \code{TRUE} verbose output is printed.
    Defaults to \code{FALSE}.}
  \item{saveb}{logical - if \code{TRUE} the values of the random effects
    are saved as part of the chain.  Default is \code{FALSE}.  Note that
    saving the random effects can require a considerable amount of
    memory.  Use with caution.}
  \item{\dots}{Some methods may take additional, optional arguments.}
}
\value{
  An object of class \code{"\linkS4class{merMCMC}"} for which many
  methods are defined.
}
\examples{
(fm1 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
set.seed(101); samp0 <- mcmcsamp(fm1, n = 1000)
str(samp0)
}
\keyword{methods}
\keyword{datagen}

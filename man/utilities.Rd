\name{prt-utilities}
\title{Print and Summary Method Utilities for Mixed Effects}
\alias{.prt.methTit}
\alias{.prt.VC}
\alias{.prt.aictab}
\alias{.prt.call}
\alias{.prt.family}
\alias{.prt.grps}
\alias{.prt.methTit}
\alias{.prt.resids}
\alias{.prt.warn}
\alias{llikAIC}
\alias{methTit}
\description{
  The \code{\link{print}}, \code{\link{summary}} methods (including the
  \code{print} for the \code{summary()} result) in \pkg{lme4} are
  modular, using about ten small utility functions.  Other packages,
  building on \pkg{lme4} can use the same utilities for ease of
  programming and consistency of output.

  \code{llikAIC()} extracts the  log likelihood, AIC, and related
  statics from a Fitted LMM.

  \code{.prt.*()} all use \code{\link{cat}} and \code{\link{print}} to
  produce output.
}
\usage{
llikAIC(object, devianceFUN = devCrit, chkREML = TRUE,
        devcomp = object@devcomp)

methTitle(dims)

.prt.methTit(mtit, class)
}
\arguments{
  %% llikAIC() :
  \item{object}{a LMM model fit}
  \item{devianceFUN}{the function to be used for computing the deviance;
    should not be changed for \pkg{lme4} created objects.}
  \item{chkREML}{optional logical indicating \code{object} maybe a REML
    fit.}
  \item{devcomp}{}

  %% methTitle():
  \item{mtit}{the result of methTitle(obj)}
  \item{class}{typically class(obj)}
}
\value{
  \code{llikAIC()} returns a \code{\link{list}} with components
  \item{logLik} and \code{AICtab} where the first is
  \code{\link{logLik}(object)} and \code{AICtab}
  is a \dQuote{table} of \code{\link{AIC}}, \code{\link{BIC}},
  \code{\link{logLik}}, deviance and \code{\link{df.residual}()}
  values.
}
\keyword{utilities}

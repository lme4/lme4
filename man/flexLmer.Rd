\name{flexLmer}
\alias{flexLmer}
\alias{flexGlmer}
\title{Flexible [g]lmer}
\usage{
flexLmer(formula, data, specials = c("cs", "d", "ar1d"), verbose = 0L)
flexGlmer(formula, data, family, specials = c("cs", "d", "ar1d"), verbose = 0L)
flexFormula <- function(formula, data, family, specials = c("cs","d","ar1d"), verbose = 0L)
}
\arguments{
  \item{formula}{Formula}
  \item{data}{Data}
  \item{family}{Family}
  \item{specials}{Special functions in \code{formula}}
  \item{verbose}{Verbose}
}
\value{
  For \code{flex[G]lmer}: \code{merMod} object.
  For \code{flexFormula}: similar to \code{[g]lFormula}.
}

\name{refitML}
\alias{refitML}
\alias{refitML.merMod}
\title{Refit a model by maximum likelihood criterion}
\usage{
  refitML(x, ...)

  \method{refitML}{merMod} (x, optimizer = "bobyqa", ...)
}
\arguments{
  \item{x}{a fitted model, usually of class
  \code{"\linkS4class{lmerMod}"}, to be refit according to
  the maximum likelihood criterion}

  \item{...}{optional additional parameters.  None are used
  at present.}

  \item{optimizer}{a string indicating the optimizer to be
  used.}
}
\value{
  an object like \code{x} but fit by maximum likelihood
}
\description{
  Refit a model using the maximum likelihood criterion
}
\details{
  This function is primarily used to get a maximum
  likelihood fit of a linear mixed-effects model for an
  \code{\link{anova}} comparison.
}


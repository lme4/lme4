\name{refitML}
\title{Refit a Model by Maximum Likelihood Criterion}
\alias{refitML}
\alias{refitML.merMod}
\description{
  Refit a (\code{merMod}) model using the maximum likelihood criterion.
}
\usage{
refitML(x, \dots)
\method{refitML}{merMod} (x, optimizer = "bobyqa", \dots)
}
\arguments{
  \item{x}{a fitted model, usually of class
    \code{"\linkS4class{lmerMod}"}, to be refit according to the maximum
    likelihood criterion.}
  \item{\dots}{optional additional parameters.  None are used
    at present.}
  \item{optimizer}{a string indicating the optimizer to be used.}
}
\value{
  an object like \code{x} but fit by maximum likelihood
}
\details{
  This function is primarily used to get a maximum likelihood fit of a
  linear mixed-effects model for an \code{\link{anova}} comparison.
}
\seealso{ \code{\link{refit}} and \code{\link{update.merMod}} for more
  extensive refitting.
}



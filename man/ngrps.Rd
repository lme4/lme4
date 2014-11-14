\name{ngrps}
\alias{ngrps}
\title{Number of Levels of a Factor or a "merMod" Model}
\description{
  Returns the number of levels of a \code{\link{factor}} or a set of
  factors, currently e.g., for each of the grouping factors of
  \code{\link{lmer}()}, \code{\link{glmer}()}, etc.
}
\usage{
ngrps(object, \dots)
}
\arguments{
  \item{object}{an \code{R} object, see Details.}
  \item{\dots}{currently ignored.}
}
\value{
  The number of levels (of a factor) or vector of number of levels for
  each \dQuote{grouping factor} of a
}
\details{
  Currently there are methods for objects of class
  \code{\linkS4class{merMod}}, i.e., the result of \code{\link{lmer}()} etc,
  and \code{\link{factor}} objects.
}
\examples{
ngrps(factor(seq(1,10,2)))
ngrps(lmer(Reaction ~ 1|Subject, sleepstudy))

## A named vector if there's more than one grouping factor :
ngrps(lmer(strength ~ (1|batch/cask), Pastes))
## cask:batch      batch
##         30         10

methods(ngrps) # -> "factor" and "merMod"
}

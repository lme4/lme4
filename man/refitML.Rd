\name{refitML}
\title{Refit a Mixed-Effect Model by Maximum Likelihood (ML)}
\alias{refitML}
\alias{refitML.merMod}
\description{
  Refits a fitted mixed-effecs model by maximum likelihood,
  typically a \emph{linear} ME model which was fit by REML.
}
\usage{
refitML(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{merMod}}.}
  \item{\dots}{optional additional arguments (none used at present)}
}
\value{
  an \R object, typically of class  \code{\linkS4class{lmerMod}},
  similar to \code{x}, but refit by ML (if needed).
}
\seealso{
  \code{\link{lmer}()}; class \code{\linkS4class{lmerMod}}.
}
\examples{
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm1M <- refitML(fm1))
## The random effects, differ, but the fixed effects are the same:
stopifnot(all.equal(fixef(fm1), fixef(fm1M)))
}
\keyword{models}

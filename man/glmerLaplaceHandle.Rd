\name{glmerLaplaceHandle}
\alias{glmerLaplaceHandle}
\title{Handle for \code{glmerLaplace}}
\usage{
glmerLaplaceHandle(pp, resp, nAGQ, tol, maxit, verbose)
}
\arguments{
  \item{pp}{\code{\link{merPredD}} object}
  \item{resp}{\code{\link{lmResp}} object}
  \item{nAGQ}{see \code{\link{glmer}}}
  \item{tol}{tolerance}
  \item{maxit}{maximum number of pwrss iterations}
  \item{verbose}{display optimizer progress}
}
\description{
  Handle for calling the \code{glmerLaplace} \code{C++} function.  Not
  intended for routine use.
}
\value{
  Value of the objective function
}
\keyword{utilities}

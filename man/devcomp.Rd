\name{devcomp}
\alias{devcomp}
\title{Extract the deviance component list}
\usage{
  devcomp(x)
}
\arguments{
  \item{x}{a fitted model of class
  \code{\linkS4class{merMod}}}
}
\value{
  a list with components \item{dims}{a named integer vector
  of various dimensions} \item{cmp}{a named numeric vector
  of components of the deviance}
}
\description{
  Return the deviance component list
}
\details{
  A fitted model of class \code{\linkS4class{merMod}} has a
  \code{devcomp} slot as described in the value section.
}
\note{
  This function is deprecated, use \code{getME(.,
  "devcomp")}
}


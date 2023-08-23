\name{namedList}
\alias{namedList}
\title{Self-naming list function}
\description{this function takes a list of arguments and combines them
  into a list; any \emph{unnamed} arguments are automatically named
  to match their symbols.  The \code{tibble::lst()} function offers
  similarly functionality.
}
\usage{
namedList(...)
}
\arguments{
  \item{...}{comma-separated arguments}
}
\examples{
a <- 1
b <- 2
c <- 3
str(namedList(a, b, d = c))
}

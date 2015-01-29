\name{factorize}
\alias{factorize}
\title{Attempt to convert grouping variables to factors}
\usage{
factorize(x,frloc)
}
\arguments{
  \item{x}{a formula}
  \item{frloc}{a data frame}
}
\description{
  If variables within a data frame are not factors, try to
  convert them.  Not intended for end-user use; this is
  a utility function that needs to be exported, for technical reasons.
}
\value{
  a copy of the data frame with factors converted
}
\keyword{models}
\keyword{utilities}


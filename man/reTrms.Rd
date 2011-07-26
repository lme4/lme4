\name{reTrms}
\alias{reTrms}
\docType{class}
\title{
  An C++ class for information related to random effects based on model
  formula terms.
}
\description{
  This C++Class encapsulates the information associated with the random
  effects in a mixed-effects model from random effects terms in the
  model formula.
}
\section{Constructors}{
  \describe{
    \item{\code{reModule(Rcpp::S4)}:}{from an S4
      \code{\linkS4class{reModule}} object}
  }
}
\section{Fields}{
  \describe{
    \item{\code{assign}}{integer - assignment of terms to grouping
      factors (read only)}
    \item{\code{cnms}}{list of column names per term (read only)}
    \item{\code{flist}}{list of grouping factors (read only)}
    \item{\code{ncols}}{integer - number of columns per term (read
      only)}
    \item{\code{nctot}}{integer - total number of columns per grouping
      factor (read only)} 
    \item{\code{nlevs}}{integer - number of levels per grouping factor
      (after removing unused levels) (read only)}
    \item{\code{offsets}}{integer - offsets per term into random effects
      (read only)}
  }
}
\section{Methods}{
  \describe{
    \item{\code{condVar(scale)}}{returns a list of 3D arrays, argument is scalar scale factor.}
    \item{\code{terms(i)}}{returns the terms associated with a grouping
      factor (argument is a 0-based index)} 
  }
}
%\details{}
%\source{}
%\references{}
%\examples{}
\keyword{datasets}

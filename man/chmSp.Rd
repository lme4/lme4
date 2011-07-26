\name{chmSp}
\alias{chmSp}
\docType{class}
\title{
  An C++ class for compressed column-oriented sparse matrices.
}
\description{
  A C++Class for \code{\linkS4class{dgCMatrix}} objects.
}
\section{Constructors}{
  \describe{
    \item{\code{chmSp(Rcpp::S4)}:}{from an 
      \code{\linkS4class{dgCMatrix}} object}
  }
}
\section{Fields}{
  \describe{
    \item{\code{ncol}}{integer - number of columns (read only)}
    \item{\code{nnz}}{integer - number of non-zeros (read only)}
    \item{\code{nrow}}{integer - number of rows (read only)}
    \item{\code{nzmax}}{integer - maximum number of non-zeros (read
      only)}
  }
}
%\details{}
%\source{}
%\references{}
%\examples{}
\keyword{datasets}

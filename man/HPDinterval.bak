\name{HPDinterval}
\alias{HPDinterval}
\docType{genericFunction}
\alias{HPDinterval,merMCMC-method}
\alias{HPDinterval,matrix-method}
\title{Highest Posterior Density intervals}
\description{
  Create Highest Posterior Density (HPD) intervals for the parameters in
  an MCMC sample.
}
\usage{
HPDinterval(object, prob = 0.95, \dots)
}
\arguments{
  \item{object}{The object containing the MCMC sample - usually of class
    \code{\linkS4class{merMCMC}} or a numeric matrix.}
  \item{prob}{A numeric scalar in the interval (0,1) giving the target
    probability content of the intervals.  The nominal probability
    content of the intervals is the multiple of \code{1/nrow(obj)}
    nearest to \code{prob}.}
  \item{\dots}{Optional additional arguments for methods.  None are used
    at present.}
}
\details{
  For each parameter the interval is constructed from the empirical cdf
  of the sample as the shortest interval for which the difference in
  the ecdf values of the endpoints is the nominal probability.  Assuming
  that the distribution is not severely multimodal, this is the HPD interval.
}

\value{
  For an \code{\linkS4class{merMCMC}} object, a list of matrices with
  columns \code{"lower"} and \code{"upper"} and rows corresponding to
  the parameters.  The attribute \code{"Probability"} is the nominal
  probability content of the intervals.  
}
%\author{Douglas Bates}
%\examples{}
\keyword{univar}
\keyword{htest}

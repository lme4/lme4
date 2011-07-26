\name{glmFamily}
\alias{glmFamily}
\docType{class}
\title{
 An C++ class for glm family objects
}
\description{
  An Rcpp module for a reference class representing a glm family
  object.  The constructor takes a \code{\link{family}} object.  The
  instance provides methods for functions in the family, making these
  available in \R{} and in \code{C++}.
}
\section{Constructors}{
  \describe{
    \item{\code{glmFamily(Rcpp::List)}:}{from an (S3)
      \code{\link{family}} object}
  }
}
\section{Fields}{
  \describe{
    \item{\code{family}}{character - name of the family}
    \item{\code{link}}{character - name of the link}
  }
}
\section{Methods}{
  \describe{
    \item{\code{devResid(mu, weights, y)}}{numeric - evalute the vector
      of deviance residuals}
    \item{\code{linkFun(mu)}}{numeric - evaluate \code{eta}}
    \item{\code{linkInv(eta)}}{numeric - evaluate \code{mu}}
    \item{\code{muEta(eta)}}{numeric - evaluate the derivative of
      \code{mu} with respect to \code{eta}}
    \item{\code{variance(mu)}}{numeric - evaluate the variance}
  }
}
\details{
  A glm \code{\link{family}} object is a list of functions including the
  link function (\code{linkfun}), the inverse link function
  (\code{linkinv}), the \code{variance} function, the \code{mu.eta}
  function that calculates the derivative of the mean with respect to
  the linear predictor, and a function that evaluates the (squared)
  deviance residuals.

  A \code{Rcpp_glmFamily} object allows these functions to be accessed
  within \code{C++} code or within \R{} code.

  For several common families the transformation functions:
  \code{linkFun}, \code{linkInv}, \code{muEta} and \code{variance}, are
  evaluated in compiled code.  }

\note{
  This class is primarily used for testing.  It is not used directly.
}
\seealso{\code{\link{glmerResp}}, \code{\link{glmer}}}
%\source{}
%\references{}
\examples{
bFam <- new(glmFamily, binomial())
bFam$family
bFam$link
set.seed(1)
mu <- sort(runif(10))	        # mean vector with elements in (0,1)
eta <- bFam$linkFun(mu)         # linear predictor
cbind(mu, eta, muEta=bFam$muEta(eta), variance=bFam$variance(mu))
}
\keyword{methods}

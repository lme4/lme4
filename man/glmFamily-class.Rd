\name{glmFamily-class}
\Rdversion{1.1}
\docType{class}
\alias{glmFamily-class}
\title{Class \code{"glmFamily"} - a reference class for \code{\link{family}}}
\description{
  This class is a wrapper class for \code{\link{family}} objects
  specifying a distibution family and link function for a generalized
  linear model (\code{\link{glm}}).  The reference class contains an
  external pointer to a C++ object representing the class.  For common
  families and link functions the functions in the family are
  implemented in compiled code so they can be accessed from other
  compiled code and for a speed boost.
}
\section{Extends}{
  All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
%\references{}
%\author{}
\note{
  Objects from this reference class correspond to objects in a C++ class.
  Methods are invoked on the C++ class using the external pointer in
  the \code{Ptr} field.  When saving such an object the external pointer
  is converted to a null pointer, which is why there is a redundant
  field \code{ptr} that is an active-binding function returning the
  external pointer.  If the \code{Ptr} field is a null pointer, the
  external pointer is regenerated for the stored \code{family} field.
}
\seealso{
  \code{\link{family}}, \code{\link{glmFamily}}
}
\examples{
str(glmFamily$new(family=poisson()))
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{Ptr}:}{Object of class \code{externalptr} ~~ }
    \item{\code{family}:}{Object of class \code{family}}
    \item{\code{ptr}:}{Object of class \code{activeBindingFunction} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{\code{devResid(mu, weights, y)}:}{returns the deviance
      residuals (numeric vector).  All arguments are numeric vectors.}
    \item{\code{variance(mu)}:}{returns the component variances from the
      mean vector, \code{mu}, and the prior weights}
    \item{\code{link(mu)}:}{applies the link function to the mean,
      \code{mu}.  The result is a numeric vector, \code{eta}.}
    \item{\code{muEta(eta)}:}{returns the diagonal of the Jacobian
      matrix of the inverse link function (numeric vector).}
    \item{\code{linkInv(eta)}:}{applies the inverse link function.}
  }
}

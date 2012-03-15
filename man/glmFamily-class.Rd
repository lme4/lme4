\docType{class}
\name{glmFamily-class}
\alias{glmFamily-class}
\title{Class \code{"glmFamily"} - a reference class for \code{\link{family}}}
\description{
  This class is a wrapper class for \code{\link{family}}
  objects specifying a distibution family and link function
  for a generalized linear model (\code{\link{glm}}).  The
  reference class contains an external pointer to a C++
  object representing the class.  For common families and
  link functions the functions in the family are
  implemented in compiled code so they can be accessed from
  other compiled code and for a speed boost.
}
\note{
  Objects from this reference class correspond to objects
  in a C++ class.  Methods are invoked on the C++ class
  using the external pointer in the \code{Ptr} field.  When
  saving such an object the external pointer is converted
  to a null pointer, which is why there is a redundant
  field \code{ptr} that is an active-binding function
  returning the external pointer.  If the \code{Ptr} field
  is a null pointer, the external pointer is regenerated
  for the stored \code{family} field.
}
\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\examples{
str(glmFamily$new(family=poisson()))
}
\seealso{
  \code{\link{family}}, \code{\link{glmFamily}}
}
\keyword{classes}


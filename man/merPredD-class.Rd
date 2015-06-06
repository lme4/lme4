\name{merPredD-class}
\alias{merPredD-class}
\title{Class \code{"merPredD"} - a Dense Predictor Reference Class}
\description{
  A reference class (see mother class definition
  \code{"\linkS4class{envRefClass}"} for a mixed-effects model predictor
  module with a dense model matrix for the fixed-effects parameters.
  The reference class is associated with a C++ class of the same name.
  As is customary, the generator object, \code{\link{merPredD}}, for the
  class has the same name as the class.
}
\note{
  Objects from this reference class correspond to objects in a C++
  class.  Methods are invoked on the C++ class object using the external
  pointer in the \code{Ptr} field.  When saving such an object the
  external pointer is converted to a null pointer, which is why there
  are redundant fields containing enough information as \R objects to be
  able to regenerate the C++ object.  The convention is that a field
  whose name begins with an upper-case letter is an \R object and the
  corresponding field, whose name begins with the lower-case letter is a
  method.  References to the external pointer should be through the
  method, not directly through the \code{Ptr} field.
}
\examples{
showClass("merPredD")
pp <- slot(lmer(Yield ~ 1|Batch, Dyestuff), "pp")
stopifnot(is(pp, "merPredD"))
str(pp) # an overview of all fields and methods' names.
}
\seealso{
  \code{\link{lmer}}, \code{\link{glmer}},
  \code{\link{nlmer}}, \code{\link{merPredD}},
  \code{\linkS4class{merMod}}.
}
\keyword{classes}


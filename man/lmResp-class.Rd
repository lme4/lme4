\name{lmResp-class}
\title{Reference Classes for Response Modules, \code{"(lm|glm|nls|lmer)Resp"}}
\alias{glmResp-class}
\alias{lmerResp-class}
\alias{lmResp-class}
\alias{nlsResp-class}
\description{
  Reference classes for response modules, including linear
  models, \code{"lmResp"}, generalized linear models,
  \code{"glmResp"}, nonlinear models, \code{"nlsResp"} and
  linear mixed-effects models, \code{"lmerResp"}. Each
  reference class is associated with a C++ class of the
  same name.  As is customary, the generator object for
  each class has the same name as the class.
}
\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.  Furthermore,
  \code{"glmResp"}, \code{"nlsResp"} and \code{"lmerResp"}
  all extend the \code{"lmResp"} class.
}
\note{
  Objects from these reference classes correspond to
  objects in C++ classes.  Methods are invoked on the C++
  classes using the external pointer in the \code{ptr}
  field.  When saving such an object the external pointer
  is converted to a null pointer, which is why there are
  redundant fields containing enough information as R
  objects to be able to regenerate the C++ object.  The
  convention is that a field whose name begins with an
  upper-case letter is an R object and the corresponding
  field whose name begins with the lower-case letter is a
  method.  Access to the external pointer should be through
  the method, not through the field.
}
\examples{
showClass("lmResp")
str(lmResp$new(y=1:4))
showClass("glmResp")
str(glmResp$new(family=poisson(), y=1:4))
showClass("nlsResp")
showClass("lmerResp")
str(lmerResp$new(y=1:4))
}
\seealso{
  \code{\link{lmer}}, \code{\link{glmer}},
  \code{\link{nlmer}}, \code{\linkS4class{merMod}}.
}
\keyword{classes}


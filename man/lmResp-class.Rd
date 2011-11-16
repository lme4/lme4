\name{lmResp-class}
\alias{lmResp-class}
\alias{glmResp-class}
\alias{lmerResp-class}
\alias{nlsResp-class}
\title{Classes \code{"lmResp"}, \code{"glmResp"}, \code{"nlsResp"} and \code{"lmerResp"}}
\description{
  Reference classes for response modules, including linear models,
  \code{"lmResp"}, generalized linear models, \code{"glmResp"},
  nonlinear models, \code{"nlsResp"} and linear mixed-effects models,
  \code{"lmerResp"}.  Each reference class is associated with a C++
  class of the same name.  As is customary, the generator object for
  each class has the same name as the class.
}
\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.  Furthermore, \code{"glmResp"},
  \code{"nlsResp"} and \code{"lmerResp"} all extend the \code{"lmResp"}
  class. 
}
%\references{}
%\author{}
\note{
  Objects from these reference classes correspond to objects in C++ classes.
  Methods are invoked on the C++ classes using the external pointer in
  the \code{ptr} field.  When saving such an object the external pointer
  is converted to a null pointer, which is why there are redundant
  fields containing enough information as R objects to be able to
  regenerate the C++ object.  The convention is that a field whose name begins
  with an upper-case letter is an R object and the corresponding field
  whose name begins with the lower-case letter is an active-binding
  function.  The values should be accessed through the active-binding
  function, not directly as the R object field.
}
\seealso{
  \code{\link{lmer}}, \code{\link{glmer}}, \code{\link{nlmer}}, \code{\linkS4class{merMod}}.
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
\keyword{classes}
\section{Fields}{
  For class \code{"lmResp"},
  \describe{
    \item{\code{y}:}{\code{numeric} response vector }
    \item{\code{Ptr}:}{\code{externalptr} to the instance of the C++ class}
    \item{\code{ptr}:}{\code{activeBindingFunction} that returns the
      external pointer, \code{Ptr}, regenerating it if necessary.  All
      access to the external pointer should be through this field.}
    \item{\code{Offset}:}{\code{numeric} vector representing the offset
      in the model.  A numeric vector of length 0 indicates no offset.}
    \item{\code{offset}:}{\code{activeBindingFunction} returning the
      offset.  This always returns a vector of the same length as the
      linear predictor.}
    \item{\code{Weights}:}{\code{numeric} vector of (prior) case
      weights. A numeric vector of length 0 indicates equal prior weights.}
    \item{\code{weights}:}{\code{activeBindingFunction} returning the
      prior case weights.  This always returns a vector of the same
      length as the response.}
  }
  Additional fields in class \code{"glmResp"},
  \describe{
    \item{\code{family}:}{A glm \code{family}, see \code{\link{family}}}
    \item{\code{N}:}{\code{integer} vector of number of individual
      observations per collapsed observation.  Only meaningful for the
      \code{\link{binomial}} family.}
    \item{\code{n}:}{\code{activeBindingFunction} to extract the
      \code{n} vector.  Always returns an integer vector of the same
      length as the response.}
  }
}
\section{Methods}{
  For class \code{"lmResp"},
  \describe{
    \item{\code{sqrtrwt()}:}{returns the square root of the residual weights}
    \item{\code{updateMu(gamma)}:}{update the mean response given the
      linear predictor \code{gamma}.  Returns the weighted residual sum of
      squares.}
    \item{\code{wrss()}:}{returns the weighted residual sum of
      squares.}
    \item{\code{fitted()}:}{returns the current value of the mean response.}
    \item{\code{wtres()}:}{returns the vector of the weighted residuals.}
    \item{\code{sqrtXwt()}:}{returns a vector of weights for the
      predictor model matrices.}
  }
  Additional methods for class \code{"glmResp"}
  \describe{
    \item{\code{sqrtWrkWt()}:}{\code{numeric} matrix of one column
      giving the square root of the working weights for the model matrices.}
    \item{\code{resDev()}:}{the residual deviance (i.e. the sum of the
      deviance residuals). }
    \item{\code{wrkResids()}:}{\code{numeric} vector of working residuals.}
    \item{\code{link()}:}{\code{character} string - the name of link.} 
    \item{\code{allInfo()}:}{a data frame whose columns are all the
      quantities of interest that are represented as numeric vectors of
      the same length as the response.}
    \item{\code{updateWts()}:}{cause the residual weights and the model
      matrix weights to be updated to the values corresponding to the
      current mean response.}
    \item{\code{muEta()}:}{the diagonal of the Jacobian matrix
      \eqn{d\mu/d\eta}{d mu/d eta}.}
    \item{\code{variance()}:}{\code{numeric} vector of case variances.}
    \item{\code{devResid()}:}{\code{numeric} vector of deviance residuals.}
    \item{\code{fam()}:}{\code{character} stirng - the name of the family.}
    \item{\code{wrkResp()}:}{\code{numeric} vector of the working responses.}
    \item{\code{eta()}:}{\code{numeric} vector of the linear predictor values.}
    \item{\code{Laplace(ldL2, ldRX2, sqrL)}:}{returns the Laplace
      approximation to the deviance.  The arguments are the logarithm of
      the square of the determinant of the sparse Cholesky factor, \eqn{\mathbf{L}}{L},
      the logarithm of the square of the determinant of the sparse Cholesky
      factor, \eqn{\mathbf{R}_X}{RX}, (not used) and the squared length of
      the orthogonal random effects vector, \eqn{\mathbf{u}}{u}.}
  }
}

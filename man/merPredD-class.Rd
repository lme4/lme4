\name{merPredD-class}
\alias{merPredD-class}
\title{Class \code{"merPredD"} -  a dense predictor reference class}
\description{
  A reference class for a mixed-effects model predictor module with a
  dense model matrix for the fixed-effects parameters.  The reference
  class is associated with a C++ class of the same name.  As is
  customary, the generator object, \code{\link{merPredD}}, for the class
  has the same name as the class.
}
\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
%\references{}
%\author{}
\note{
  Objects from this reference class correspond to objects in a C++ class.
  Methods are invoked on the C++ class object using the external pointer in
  the \code{ptr} field.
  % When saving such an object the external pointer
  % is converted to a null pointer, which is why there are redundant
  % fields containing enough information as R objects to be able to
  % regenerate the C++ object.  The convention is that a field whose name begins
  % with an upper-case letter is an R object and the corresponding field
  % whose name begins with the lower-case letter is an active-binding
  % function.  The values should be accessed through the active-binding
  % function, not directly as the R object field.
}
\seealso{
  \code{\link{lmer}}, \code{\link{glmer}}, \code{\link{nlmer}},
  \code{\link{merPredD}}, \code{\linkS4class{merMod}}.
}
\examples{
showClass("merPredD")
str(lmer(Yield ~ 1|Batch, Dyestuff)@pp)
}
\keyword{classes}
\section{Fields}{
  For class \code{"lmResp"},
  \describe{
    \item{\code{X}:}{dense model matrix of class \code{\linkS4class{ddenseModelMatrix}}}
    \item{\code{Zt}:}{transposed sparse model matrix of class \code{\linkS4class{dgCMatrix}}}
    \item{\code{Lambdat}:}{transposed sparse relative covariance factor
      of class \code{\linkS4class{dgCMatrix}}}
    \item{\code{Lind}:}{integer vector of 1-based indices into the
      \code{theta} field.}
    \item{\code{theta}:}{numeric vector of variance component
      parameters.  Assignment of this field causes an update of values
      in the C++ class.}
    \item{\code{ptr}:}{\code{externalptr} to the instance of the C++ class}
    % \item{\code{ptr}:}{\code{activeBindingFunction} that returns the
    %   external pointer, \code{Ptr}, regenerating it if necessary.  All
    %   access to the external pointer should be through this field.}
    \item{\code{beta0}:}{numeric vector of offsets for the fixed-effects
      coefficients}
    \item{\code{u0}:}{numeric vector of offsets for the orthogonal
      random effects}
  }
}
\section{Methods}{
  \describe{
    \item{\code{b()}:}{returns the conditional mode of the random
      effects (numeric vector).}
    \item{\code{beta()}:}{returns the fixed-effects parameters (numeric
      vector).}
    \item{\code{CcNumer()}:}{returns the numerator of the orthogonality
      convergence criterion (numeric scalar).}
    \item{\code{delb()}:}{returns the full increment for the
      fixed-effects parameters (numeric vector).}
    \item{\code{delu()}:}{returns the full increment for the
      orthogonal random effects (numeric vector).}
    \item{\code{getLambdat()}:}{returns the relative covariance factor
      at the current value of \code{theta} (a
      \code{\linkS4class{dgCMatrix}} object).  The \code{Lambdat} field
      is not updated with each change in the \code{theta} field but the
      value returned by this method is.}
    \item{\code{initialize()}:}{the standard \code{initialize} method is
      overridden in this class.}
    \item{\code{installPars(fac)}:}{updates the offset for the
      fixed-effects parameters and the orthogonal random effects using
      step factor \code{fac}.}
    \item{\code{L()}:}{returns the sparse Cholesky factor for the random
      effects (a \code{\linkS4class{dCHMsimpl}} object).}
    \item{\code{LamtUt()}:}{returns the product of \code{Lambdat} and
      \code{Ut}, a \code{linkS4class{dgCMatrix}} object.}
    \item{\code{ldL2()}:}{returns the log of the square of the
      determinant of the sparse Cholesky factor, \code{L} (numeric
      scalar).}
    \item{\code{ldRX2()}:}{returns the log of the square of the
      determinant of the dense Cholesky factor, \code{RX} (numeric
      scalar).}
    \item{\code{linPred(fac)}:}{returns the linear predictor for step
      factor \code{fac} (numeric vector).}
    \item{\code{P()}:}{returns the 0-based permutation vector of the
      sparse Cholesky factor (integer vector).}
    \item{\code{RX()}:}{returns the upper-triangular dense Cholesky
      factor for the random effects (numeric matrix).}
    \item{\code{RXdiag()}:}{returns the diagonal of \code{RX} (numeric
      vector)}
    \item{\code{RZX()}:}{returns the off-diagonal block of full Cholesky
      factor (numeric matrix).}
    \item{\code{solve()}:}{updates both increments, \code{delu} and
      \code{delb}, as the least squares solutions given the updated
      residuals and X-weights.}
    \item{\code{solveU()}:}{updates only the increment \code{delu}
      the least squares solutions given the updated residuals and
      X-weights.}
    \item{\code{u()}:}{returns the conditional mode of the orthogonal
      random effects (numeric vector).}
    \item{\code{unsc()}:}{returns the unscaled variance-covariance
      matrix of the fixed-effects parameters (numeric matrix).}
    \item{\code{updateDecomp()}:}{causes an update of \code{L},
      \code{RZX} and \code{RX}, the parts of the full Cholesky
      decomposition in the C++ object.  Also updates \code{ldL2} and
      \code{ldRX2} (returns NULL).}
    \item{\code{updateRes(wtres)}:}{updates the C++ object for the
      numeric vector of weighted residuals, \code{wtres} (returns
      NULL).}
    \item{\code{updateXwts(wts)}:}{updates the C++ object for the
      numeric matrix of square roots of the X weights \code{wts}
      (returns NULL).}
    \item{\code{Ut()}:}{the transposed, weighted random effects model
      matrix, a \code{\linkS4class{dgCMatrix}} object.}
    \item{\code{Utr()}:}{the product of \code{Ut} and the weighted
      residuals (numeric vector).}
    \item{\code{V()}:}{the weighted fixed effects model matrix (numeric
      matrix).}
    \item{\code{Vtr()}:}{the cross-product of \code{V} and the weighted
      residuals (numeric vector).}
    \item{\code{VtV()}:}{the cross-product of \code{V} (numeric
      matrix).}
  }
}

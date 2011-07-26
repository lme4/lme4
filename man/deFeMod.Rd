\name{deFeMod}
\alias{deFeMod}
\alias{coerce,Rcpp_deFeMod,deFeMod-method}
\alias{show,Rcpp_deFeMod-method}
\docType{class}
\title{
  An C++ class for information related to the fixed-effects parameters.
}
\description{
  This C++Class encapsulates the information associated with the fixed
  effects in a mixed-effects model.
}
\section{Constructors}{
  \describe{
    \item{\code{deFeMod(Rcpp::S4, int)}:}{from an S4
      \code{\linkS4class{deFeMod}} object and the number of
      observations, \code{n}}
    \item{\code{deFeMod(X, n, q)}}{from the dense model matrix, \code{X}
      and dimensions}
  }
}
\section{Read-Only Fields}{
  \describe{
    \item{\code{CcNumer}}{numeric scalar - contribution to numerator of
      the convergence criterion}
    \item{\code{RX}}{\code{\linkS4class{Cholesky}} - dense Cholesky
      factor}
    \item{\code{RZX}}{\code{\linkS4class{dgeMatrix}} - Off-diagonal
      block in large Cholesky factor}
    \item{\code{UtV}}{\code{\linkS4class{dgeMatrix}} - Weighted
      crossproduct}
    \item{\code{V}}{\code{\linkS4class{dgeMatrix}} - Scaled model matrix}
    \item{\code{VtV}}{\code{\linkS4class{dpoMatrix}} - Weighted
      crossproduct}
    \item{\code{Vtr}}{numeric - Weighted crossproduct of model matrix
      and residuals}
    \item{\code{X}}{\code{\linkS4class{dgeMatrix}} - original model
      matrix}
    \item{\code{coef}}{numeric - coefficient vector after increment}
    \item{\code{ldRX2}}{numeric - log of the square of the determinant
      of \code{RX}}
  }
}
\section{Read-Write Fields}{
  \describe{
    \item{\code{coef0}}{numeric - base coefficient vector}
    \item{\code{incr}}{numeric - increment for coefficient vector}
  }
}
\section{Methods}{
  \describe{
    \item{\code{linPred1(fac)}}{numeric - update \code{coef} to
      \code{coef0 + fac*incr} and return linear predictor contribution}
    \item{\code{reweight(sqrtXwt, wtres)}}{update V and Vtr for new
      weights.}
    \item{\code{solveIncr()}}{update \code{fac} from \code{VtV} and
      solve for \code{incr} from \code{Vtr}}
    \item{\code{updateIncr(cu)}}{update \code{incr} given 
      \code{cu}}
    \item{\code{updateRzxpRxp(Lambdap, Lp)}}{update the
      triangular factor sections given external pointers to
      \code{Lambda} and \code{L}}
    \item{\code{updateUtV(Utp)}}{update \code{UtV} given an external
      pointer to \code{Ut}}
  }
}
%\details{}
%\source{}
%\references{}
%\examples{}
\keyword{datasets}

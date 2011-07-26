\name{reModule}
\alias{reModule}
\alias{coerce,Rcpp_reModule,reModule-method}
\alias{show,Rcpp_reModule-method}
\docType{class}
\title{
  An C++ class for information related to random effects
}
\description{
  This C++Class encapsulates the information associated with the random
  effects in a mixed-effects model.
}
\section{Constructors}{
  \describe{
    \item{\code{reModule(Rcpp::S4)}:}{from an S4
      \code{\linkS4class{reModule}} object}
    \item{\code{reModule(Zt, Lambda, L, Lind, lower)}}{from components}
  }
}
\section{Read-Only Fields}{
  \describe{
    \item{\code{CcNumer}}{numeric scalar - contribution to numerator of
      the convergence criterion}
    \item{\code{L}}{\code{\linkS4class{CHMfactor}} - sparse Cholesky
      factor} 
    \item{\code{Lambda}}{\code{\linkS4class{dgCMatrix}} - relative
      covariance factor} 
    \item{\code{Lambdap}}{external pointer to \code{Lambda}}
    \item{\code{Lind}}{integer - 1-based index vector into theta for
      \code{Lambda@x}}
    \item{\code{Lp}}{external pointer to \code{L}}
    \item{\code{Ut}}{external pointer to weighted, orthogonal design
      matrix}
    \item{\code{Zt}}{\code{\linkS4class{dgCMatrix}} - transpose of the
      model matrix for the random effects}
    \item{\code{b}}{numeric - random effects on original scale}
    \item{\code{cu}}{numeric - intermediate solution for \code{u}}
    \item{\code{ldL2}}{numeric - logarithm of the square of the
      determinant of \code{L}}
    \item{\code{linPred}}{numeric - linear predictor contribution based
      on current \code{u}}
    \item{\code{lower}}{numeric - lower bounds on components of
      \code{theta}}
    \item{\code{sqrLenU}}{numeric - squared length of u, the orthogonal
      random effects}
    \item{\code{u}}{numeric - orthogonal random effects vector}
  }
}
\section{Read-Write Fields}{
  \describe{
    \item{\code{incr}}{numeric - increment for \code{u}}
    \item{\code{theta}}{numeric - variance component parameters}
    \item{\code{u0}}{numeric - base orthogonal random effects vector}
  }
}
\section{Methods}{
  \describe{
    \item{\code{installU0()}}{install the current value of \code{u} as \code{u0}}
    \item{\code{linPred1(fac)}}{update \code{u} as \code{u0 + fac*incr}
      and return linear predictor contribution}
    \item{\code{reweight(sqrtXwt, wtres, useU0)}}{update L, Ut and cu
      for new weights.  The value of \code{u} is taken from \code{u0} or
      from \code{u0} according to whether or not \code{useU0} is TRUE}
    \item{\code{solveIncr()}}{solve for \code{incr} only. Returns squared length of \code{c1}}
    \item{\code{updateIncr(cu)}}{solve for \code{incr} given the updated
      \code{cu} from the feModule's \code{updateBeta} method}
  }
}
%\details{}
%\source{}
%\references{}
%\examples{}
\keyword{datasets}

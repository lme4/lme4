\name{lmerResp}
\alias{lmerResp}
\alias{coerce,Rcpp_lmerResp,lmerResp-method}
\alias{show,Rcpp_lmerResp-method}
\docType{class}
\title{
  An C++ class for \code{lmerResp} objects
}
\description{
  The \code{\linkS4class{lmerResp}} S4 class represents the information
  associated with the response in an \code{\link{lm}} or
  \code{\link{lmer}} model. This C++Class encapsulates the information
  and provides methods accessible in \R{} or in \code{C++}.
}
\section{Constructors}{
  \describe{
    \item{\code{lmerResp(Rcpp::S4)}:}{from a
      \code{\linkS4class{lmerResp}} object}
    \item{\code{lmerResp("integer", "numeric")}}{from the REML indicator
      and the response \code{y}}
    \item{\code{lmerResp("integer", "numeric", "numeric")}}{from the REML
      indicator, the response \code{y} and the \code{weights}}
    \item{\code{lmerResp("integer", "numeric", "numeric",
	"numeric")}}{from the REML indicator, the response \code{y}, the
      \code{weights} and the \code{offset}}
  }
}
\section{Read-Only Fields}{
  \describe{
    \item{\code{REML}}{scalar integer - 0L for ML estimation and p for
      REML}
    \item{\code{mu}}{numeric - mean vector (read only)}
    \item{\code{offset}}{numeric - the offset vector, always present
      even if zero (read/write)}    
    \item{\code{sqrtXwt}}{numeric - square root of the weights applied
      to the model matrix (read only)}
    \item{\code{sqrtrwt}}{numeric - square root of the weights applied
      to the residuals (read only)}
    \item{\code{weights}}{numeric - prior weights (read/write)}
    \item{\code{wrss}}{numeric scalar - weighted sum of squared 
      residuals (read only)}
    \item{\code{wtres}}{numeric - weighted residuals (read only)}
    \item{\code{y}}{numeric - response vector (read only)}
  }
}
\section{Methods}{
  \describe{
    \item{\code{Laplace(ldL2, ldRX2, sqrLenU)}}{Return the profiled REML
      criterion or the profiled deviance, according to the setting of
      the REML field}
    \item{\code{updateMu(eta)}}{Update \code{mu} and derived quantities
      given a value of the linear predictor, \code{eta}.}
    \item{\code{updateWts()}}{Update \code{sqrtrwt} and
      \code{sqrtXwt}. (Not generally used for this class.)}
  }
}
\details{
  The methods in this class parallel those in the
  \code{"\link{glmerResp}"} class but typically \code{updateWts} is
  not called.  The \code{updateMu} method is called after every change
  in the variance component parameter, providing an updated value of the
  residual sum of squares.

  The first argument to the constructors (other than from an S4 object)
  is the REML indicator, which is an integer value of 0, for
  ML estimates, or \code{p}, the rank of the fixed-effects model matrix,
  for REML estimates.
}
%\source{}
%\references{}
\examples{
ll <- with(sleepstudy, new(lmerResp, 0L, Reaction))
comp <- c("y", "mu", "weights", "offset", "sqrtrwt", "wtres")
names(comp) <- comp
head(sapply(comp, function(nm) ll[[nm]]), n = 10)
ll$wrss
}
\keyword{methods}

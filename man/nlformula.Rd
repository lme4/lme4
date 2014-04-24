\name{nlformula}
\alias{nlformula}
\title{Manipulate a Nonlinear Model Formula}
\description{
  Check and manipulate the formula for a nonlinear model, such as
  specified in \code{\link{nlmer}}.
}
\usage{
  nlformula(mc)
}
\arguments{
  \item{mc}{matched call from the calling function, typically
    \code{\link{nlmer}()}.  Should have arguments named
    \describe{
      \item{formula:}{a formula of the form \code{resp ~ nlmod ~ meform}
	where \code{resp} is an expression for the response,
	\code{nlmod} is the nonlinear model expression and
	\code{meform} is the mixed-effects model formula.
	\code{resp} can be omitted when, e.g., optimizing a design.}
      \item{data:}{a data frame in which to evaluate the model function}
      \item{start:}{either a numeric vector containing initial estimates
	for the nonlinear model parameters or a list with components
	\describe{
	  \item{\code{nlpars}:}{the initial estimates of the nonlinear
	    model parameters}
	  \item{\code{theta}:}{the initial estimates of the variance
	    component parameters} 
	}
      }
    }
  }
}
\value{
  a list with components
  \item{"respMod"}{a response module of class \code{"\linkS4class{nlsResp}"}} 
  \item{"frame"}{the model frame, including a terms attribute}
  \item{"X"}{the fixed-effects model matrix}
  \item{"reTrms"}{the random-effects terms object}
}
\details{
  The model formula for a nonlinear mixed-effects model is
  of the form \code{resp ~ nlmod ~ mixed} where
  \code{resp} is an expression (usually just a name) for
  the response, \code{nlmod} is the call to the nonlinear
  model function, and \code{mixed} is the mixed-effects
  formula defining the linear predictor for the parameter
  matrix.  If the formula is to be used for optimizing
  designs, the \code{resp} part can be omitted.
}
\seealso{
  Other utilities: \code{\link{findbars}},
  \code{\link{mkRespMod}}, \code{\link{mkReTrms}},
  \code{\link{nobars}}, \code{\link{subbars}}
}


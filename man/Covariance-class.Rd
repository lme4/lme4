\docType{class}
\name{Covariance-class}
\alias{Covariance-class}
\alias{Covariance.us-class}
\alias{Covariance.diag-class}
\alias{Covariance.cs-class}
\alias{Covariance.ar1-class}

\title{Virtual Class \sQuote{Covariance} of Covariance Matrices}

\description{
  An \code{S4} class that describes the within-group variance-covariance 
  structures available for the \CRANpkg{lme4} package. The syntax for use is 
  similar to the \CRANpkg{glmmTMB} package. See the \emph{Examples} section 
  on how to use.
  
  There are four main covariance classes: \code{Covariance.us} (unstructured), 
  \code{Covariance.diag} (diagonal), \code{Covariance.cs} (compound symmetry),
  and \code{Covariane.ar1} (autoregressive order 1).
}

\usage{
\S4method{getPar}{Covariance}(object)
\S4method{getParLength}{Covariance}(object)
\S4method{getParNames}{Covariance}(object, cnm, gnm)
\S4method{setPar}{Covariance}(object, value)
\S4method{getTheta}{Covariance}(object)
\S4method{getThetaLength}{Covariance}(object)
\S4method{getThetaNames}{Covariance}(object, cnm, gnm)
\S4method{getThetaIndex}{Covariance}(object)
\S4method{getThetaIndexLength}{Covariance}(object)
\S4method{setTheta}{Covariance}(object, value)
\S4method{getLower}{Covariance}(object)
\S4method{getUpper}{Covariance}(object)
\S4method{getLambda}{Covariance}(object)
\S4method{getLambdat.dp}{Covariance}(object)
\S4method{getLambdat.i}{Covariance}(object)
\S4method{getVC}{Covariance}(object)
\S4method{getVCNames}{Covariance}(object, cnm, gnm)
\S4method{setVC}{Covariance}(object, vcomp, ccomp)
\S4method{getProfPar}{Covariance}(object)
\S4method{setProfPar}{Covariance}(object, value)
\S4method{getProfLower}{Covariance}(object)
\S4method{getProfUpper}{Covariance}(object)
}

\details{
  In matrix notation, a linear mixed model can be represented as:
  \deqn{
  \mathbf{y} = X \beta + Z \mathbf{b} + \boldsymbol{\epsilon}
  }
  Where \eqn{\mathbf{b}} represents an unknown vector of random effects. The
  \code{Covariance} class helps define the structure of the 
  variance-covariance matrix as specified as \eqn{Var(\mathbf{b})}. Typically, 
  we denote the variance-covariance matrix as \eqn{\Sigma}.
  
  Proper details can be found in \code{vignette("lmer")}.
  Our main focus is the restriction on \eqn{\Sigma}. Let \eqn{q} represent 
  the number of columns of the random-effects model matrix \eqn{Z}. 
  \eqn{\Sigma} can be expressed in terms of a relative cofactor 
  \eqn{\Lambda_{\theta}} which is a \eqn{q \times q} that depends on the 
  variance-component parameter \eqn{\theta}:
  \deqn{
  \mathbf{\Sigma}_{\mathbf{\theta}} = \sigma^{2} \Lambda_{\theta} 
  \Lambda_{\theta}^{T}.
  }
  Where \eqn{\sigma} is a scale parameter of the variance of a linear mixed model.
  
  If the within-group variance-covariance structure is unspecified, the default 
  is unstructured (\code{Covariance.us}), meaning we only ensure that the 
  variance-covariance matrix \eqn{\Sigma} to be \emph{positive definite}. 
  
}
\section{Objects from the Class}{
Available standard classes:
  \describe{
    \item{\code{Covariance.us}}{Unstructured (general positive definite). This is the 
    default version.}
    \item{\code{Covariance.diag}}{Diagonal; only the diagonal entries are 
    nonzero, indicating no covariances between variables. If \code{hom = TRUE}, 
    all diagonal entries (variances) are equal across groups.}
    \item{\code{Covariance.cs}}{Compound symmetry. If \code{hom = TRUE}, 
    all diagonal entries (variances) are equal across groups.}
    \item{\code{Covariance.ar1}}{Autoregessive process of order 1. If 
    \code{hom = TRUE}, then... TODO}
  }
Besides \code{Covariance.us}, the remaining classes contain a logical slot 
called \code{hom}, described in the \emph{Slots} section.
}
\section{Slots}{
  \describe{
    \item{\code{nc}}{An integer value giving the number of columns (or components).}
    \item{\code{par}}{A double vector that stores values of \code{theta}, which 
    specifies the covariance parameters for the model.}
    \item{\code{hom}}{A logical which represents whether the variance-covariance 
    matrix has a \emph{homogenous} structure.}
  }
}

\section{Methods}{
\code{vcomp} contains parameters of the variance component of the 
covariance matrix (typically standard deviations). \code{ccomp} contains 
parameters of the correlation component of the covariance matrix. That is,
\itemize{
  \item \code{Covariance.us}: the entries of the lower triangle of the 
  correlation matrix in column-major order.
  \item \code{diag}: should be empty; we only care about the standard deviations.
  \item \code{cs}: \eqn{\rho}
  \item \code{ar1}: \eqn{\rho}
  }
}



\section{S4 Methods}{
\describe{
    \item{\code{getPar}:}{The vector of parameters that the optimizer cares 
    about; contains \code{c(vcomp, ccomp)}.}
    \item{\code{getParLength}:}{Obtains the length of par.}
    \item{\code{getParNames}:}{TODO.}
    \item{\code{setPar}:}{Sets the value of par.}
    \item{\code{getTheta}:}{Obtains the variance-component parameter. These are 
    also the column-wise elements of \code{Lambda}, or the row-wise elements of
    \code{t(Lambda)}.}
    \item{\code{getThetaLength}:}{Obtains the length of the variance-component 
    parameter.}
    \item{\code{getThetaNames}:}{TODO.}
    \item{\code{getThetaIndex}:}{TODO.}
    \item{\code{getThetaIndexLength}:}{TODO.}
    \item{\code{setTheta}:}{Sets the value of the variance-component parameter.}
    \item{\code{getLower}:}{Retrieves the lower bound for \code{vcomp} and 
    \code{ccomp}, in that order.}
    \item{\code{getUpper}:}{Retrieves the upper bound for \code{vcomp} and 
    \code{ccomp}, in that order.}
    \item{\code{getLambda}:}{The relative covariance factor.}
    \item{\code{getLambdat.dp}:}{TODO.}
    \item{\code{getLambda.i}:}{TODO.}
    \item{\code{getVC}:}{Retrieves both \code{vcomp} and \code{ccomp}.}
    \item{\code{getVCNames}:}{TODO.}
    \item{\code{setVC}:}{Sets the values for VC, which is composed of \code{vcomp}
    and \code{ccomp}.}
    \item{\code{getProfPar}:}{TODO.}
    \item{\code{setProfPar}:}{TODO.}
    \item{\code{getProfLower}:}{TODO.}
    \item{\code{getProfUpper}:}{TODO.}
  }
}


\examples{
## Unstructured
fm1 <- lmer(Reaction ~ Days + us(Days | Subject), sleepstudy)
## Diagional
fm1 <- lmer(Reaction ~ Days + diag(Days | Subject), sleepstudy)
## Compound symmetry
fm1 <- lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy)
## TODO: EXAMPLE FOR AR1
}
\keyword{classes}
\keyword{internal}

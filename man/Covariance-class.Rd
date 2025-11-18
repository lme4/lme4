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
}

\details{
  In matrix notation, a linear mixed model can be represented as:
  \deqn{
  \mathbf{y} = X \beta + Z \mathbf{b} + \boldsymbol{\epsilon}
  }
  Where \eqn{\mathbf{b}} represents an unknown vector of random effects. The
  \code{Covariance} class helps define the structure of the 
  variance-covariance matrix as specified as \eqn{Var(\mathbf{b})}.
  
  If the within-group variance-covariance structure is unspecified, the default 
  is unstructured (\code{Covariance.us}), meaning we only ensure that the 
  variance-covariance matrix to be \emph{positive definite}. 
}
\section{Objects from the Class}{
Available standard classes:
  \describe{
    \item{\code{Covariance.us}}{Unstructured (general positive definite). This is the 
    default version.}
    \item{\code{Covariance.diag}}{Diagonal; only the diagonal entries are 
    nonzero, indicating no covariances between variables. If \code{hom = TRUE}, 
    all diagonal entries (variances) are equal across groups.}
    \item{\code{Covariance.cs}}{Compound symmetry. (TODO)}
    \item{\code{Covariance.ar1}}{Autoregessive process of order 1. (TODO)}
  }
Besides \code{Covariance.us}, the remaining classes contain a logical slot 
called \code{hom} which additionally represents whether the 
variance-covariance matrix has a \emph{homogenous} structure.
}
\section{Slots}{
  \describe{
    \item{\code{nc}}{An integer value giving the number of columns (or components).}
    \item{\code{par}}{A double vector that stores values of \code{theta}, which 
    specifies the covariance parameters for the model.}
    \item{\code{hom}}{TODO.}
  }
}

\section{Methods}{
  \describe{
    \item{\code{getPar}}{TODO.}
    \item{\code{getParLength}}{TODO.}
    \item{\code{getParNames}}{TODO.}
    \item{\code{setPar}}{TODO.}
    \item{\code{getTheta}}{TODO.}
    \item{\code{getThetaLength}}{TODO.}
    \item{\code{getThetaNames}}{TODO.}
    \item{\code{getThetaIndex}}{TODO.}
    \item{\code{getThetaIndexLength}}{TODO.}
    \item{\code{setTheta}}{TODO.}
    \item{\code{getLower}}{TODO.}
    \item{\code{getUpper}}{TODO.}
    \item{\code{getLambda}}{TODO.}
    \item{\code{getLambda.dp}}{TODO.}
    \item{\code{getLambda.i}}{TODO.}
    \item{\code{getVC}}{TODO.}
    \item{\code{getVCNames}}{TODO.}
    \item{\code{setVC}}{TODO.}
    \item{\code{getProfPar}}{TODO.}
    \item{\code{setProfPar}}{TODO.}
    \item{\code{getProfLower}}{TODO.}
    \item{\code{getProfUpper}}{TODO.}
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

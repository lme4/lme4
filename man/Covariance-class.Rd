\docType{class}
\name{Covariance-class}
\alias{Covariance-class}
\alias{Covariance.us-class}
\alias{Covariance.diag-class}
\alias{Covariance.cs-class}
\alias{Covariance.ar1-class}
\alias{initialize,Covariance.us-method}
\alias{initialize,Covariance.diag-method}
\alias{initialize,Covariance.cs-method}
\alias{initialize,Covariance.ar1-method}
\alias{getPar}
\alias{getPar,Covariance-method}
\alias{getPar,merMod-method}
\alias{getParLength}
\alias{getParLength,Covariance-method}
\alias{getParLength,merMod-method}
\alias{getParNames}
\alias{getParNames,Covariance.us,character,character-method}
\alias{getParNames,Covariance.diag,character,character-method}
\alias{getParNames,Covariance.cs,character,character-method}
\alias{getParNames,Covariance.ar1,character,character-method}
\alias{getParNames,merMod,missing,missing-method}
\alias{setPar}
\alias{setPar,Covariance,numeric-method}
\alias{getTheta}
\alias{getTheta,Covariance.us-method}
\alias{getTheta,Covariance.diag-method}
\alias{getTheta,Covariance.cs-method}
\alias{getTheta,Covariance.ar1-method}
\alias{getTheta,merMod-method}
\alias{getThetaLength}
\alias{getThetaLength,Covariance.us-method}
\alias{getThetaLength,Covariance.diag-method}
\alias{getThetaLength,Covariance.cs-method}
\alias{getThetaLength,Covariance.ar1-method}
\alias{getThetaLength,merMod-method}
\alias{getThetaNames}
\alias{getThetaNames,Covariance.us,character,character-method}
\alias{getThetaNames,Covariance.diag,character,character-method}
\alias{getThetaNames,Covariance.cs,character,character-method}
\alias{getThetaNames,Covariance.ar1,character,character-method}
\alias{getThetaNames,merMod,missing,missing-method}
\alias{getThetaIndex}
\alias{getThetaIndex,Covariance.us-method}
\alias{getThetaIndex,Covariance.diag-method}
\alias{getThetaIndex,Covariance.cs-method}
\alias{getThetaIndex,Covariance.ar1-method}
\alias{getThetaIndexLength}
\alias{getThetaIndexLength,Covariance.us-method}
\alias{getThetaIndexLength,Covariance.diag-method}
\alias{getThetaIndexLength,Covariance.cs-method}
\alias{getThetaIndexLength,Covariance.ar1-method}
\alias{setTheta}
\alias{setTheta,Covariance.us,numeric-method}
\alias{setTheta,Covariance.diag,numeric-method}
\alias{setTheta,Covariance.cs,numeric-method}
\alias{setTheta,Covariance.ar1,numeric-method}
\alias{getLower}
\alias{getLower,Covariance.us-method}
\alias{getLower,Covariance.diag-method}
\alias{getLower,Covariance.cs-method}
\alias{getLower,Covariance.ar1-method}
\alias{getLower,merMod-method}
\alias{getUpper}
\alias{getUpper,Covariance.us-method}
\alias{getUpper,Covariance.diag-method}
\alias{getUpper,Covariance.cs-method}
\alias{getUpper,Covariance.ar1-method}
\alias{getUpper,merMod-method}
\alias{getLambda}
\alias{getLambda,Covariance.us-method}
\alias{getLambda,Covariance.diag-method}
\alias{getLambda,Covariance.cs-method}
\alias{getLambda,Covariance.ar1-method}
\alias{getLambda,merMod-method}
\alias{getLambdat.dp}
\alias{getLambdat.dp,Covariance.us-method}
\alias{getLambdat.dp,Covariance.diag-method}
\alias{getLambdat.dp,Covariance.cs-method}
\alias{getLambdat.dp,Covariance.ar1-method}
\alias{getLambdat.i}
\alias{getLambdat.i,Covariance.us-method}
\alias{getLambdat.i,Covariance.diag-method}
\alias{getLambdat.i,Covariance.cs-method}
\alias{getLambdat.i,Covariance.ar1-method}
\alias{getVC}
\alias{getVC,Covariance.us-method}
\alias{getVC,Covariance.diag-method}
\alias{getVC,Covariance.cs-method}
\alias{getVC,Covariance.ar1-method}
\alias{getVCNames}
\alias{getVCNames,Covariance.us,character,character-method}
\alias{getVCNames,Covariance.diag,character,character-method}
\alias{getVCNames,Covariance.cs,character,character-method}
\alias{getVCNames,Covariance.ar1,character,character-method}
\alias{getVCNames,merMod,missing,missing-method}
\alias{setVC}
\alias{setVC,Covariance.us,numeric,numeric-method}
\alias{setVC,Covariance.diag,numeric,numeric-method}
\alias{setVC,Covariance.cs,numeric,numeric-method}
\alias{setVC,Covariance.ar1,numeric,numeric-method}
\alias{getProfPar}
\alias{getProfPar,Covariance.us-method}
\alias{getProfPar,Covariance.diag-method}
\alias{getProfPar,Covariance.cs-method}
\alias{getProfPar,Covariance.ar1-method}
\alias{setProfPar}
\alias{setProfPar,Covariance.us,numeric-method}
\alias{setProfPar,Covariance.diag,numeric-method}
\alias{setProfPar,Covariance.cs,numeric-method}
\alias{setProfPar,Covariance.ar1,numeric-method}
\alias{getProfPar,merMod-method}
\alias{getProfLower}
\alias{getProfLower,Covariance.us-method}
\alias{getProfLower,Covariance.diag-method}
\alias{getProfLower,Covariance.cs-method}
\alias{getProfLower,Covariance.ar1-method}
\alias{getProfLower,merMod-method}
\alias{getProfUpper}
\alias{getProfUpper,Covariance.us-method}
\alias{getProfUpper,Covariance.diag-method}
\alias{getProfUpper,Covariance.cs-method}
\alias{getProfUpper,Covariance.ar1-method}
\alias{getProfUpper,merMod-method}

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
getPar(object)
getParLength(object)
getParNames(object, cnm, gnm, prf = NULL, old = TRUE)
setPar(object, value, pos = 0L)
getTheta(object)
getThetaLength(object)
getThetaNames(object, cnm, gnm, prf = NULL, old = TRUE)
getThetaIndex(object)
getThetaIndexLength(object)
setTheta(object, value, pos = 0L)
getLower(object)
getUpper(object)
getLambda(object)
getLambdat.dp(object)
getLambdat.i(object)
getVC(object)
getVCNames(object, cnm, gnm, prf = NULL, old = TRUE)
setVC(object, vcomp, ccomp)
getProfPar(object, profscale, sc = NULL)
setProfPar(object, profscale, sc = NULL, value, pos = 0L)
getProfLower(object, profscale, sc = NULL)
getProfUpper(object, profscale, sc = NULL)
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
    nonzero, indicating no covariances between variables.}
    \item{\code{Covariance.cs}}{Compound symmetry.}
    \item{\code{Covariance.ar1}}{Autoregessive process of order 1.}
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
    matrix has a \emph{homogenous} structure. That is, whether the variances 
    are equal across groups.}
  }
}

\section{S4 Methods}{
\code{vcomp} contains parameters of the variance component of the 
covariance matrix (typically standard deviations). \code{ccomp} contains 
parameters of the correlation component of the covariance matrix. That is,
\itemize{
  \item \code{Covariance.us}: the entries of the lower triangle of the 
  correlation matrix in column-major order.
  \item \code{diag}: should be empty; we only care about the standard deviations.
  \item \code{cs}, \code{ar1}: \eqn{\rho}
  }
Methods that obtain names (i.e., \code{getParNames}, \code{getThetaNames}, 
and \code{getVCNames}), \code{cnm} represents the column (of the model matrix) 
name, and \code{gnm} represents the group name.

For methods that involve \code{getProf} (i.e., \code{getProfPar}, 
\code{setProfPar}, \code{getProfLower}, \code{getProfUpper}), the user may 
set the argument \code{profscale = "varcor"}, which transforms the variance 
components andcorrelation components into the corresponding variance-covariance 
matrix. Furthermore, the user may scale the standard deviations by a factor
of the \code{sc} argument.

\describe{
    \item{\code{getPar}:}{The vector of parameters that the optimizer cares 
    about; contains \code{c(vcomp, ccomp)}.}
    \item{\code{getParLength}:}{Obtains the length of \code{par}.}
    \item{\code{getParNames}:}{Obtains the names of \code{par} in terms of 
    \code{cnm} and \code{gnm}.}
    \item{\code{setPar}:}{Sets the value of par.}
    \item{\code{getTheta}:}{Obtains the variance-component parameter. These are 
    also the column-wise elements of \code{Lambda}, or the row-wise elements of
    \code{t(Lambda)}.}
    \item{\code{getThetaLength}:}{Obtains the length of the variance-component 
    parameter \code{theta}.}
    \item{\code{getThetaNames}:}{Obtains the names of \code{theta} in terms of 
    \code{cnm} and \code{gnm}.}
    \item{\code{getThetaIndex}:}{A mapping of the elements of \code{theta} into 
    \code{Lambda} in column-major order.}
    \item{\code{getThetaIndexLength}:}{Obtains the length of \code{getThetaIndex};
    which corresponds to the number of non-zero entries of \code{Lambda}.}
    \item{\code{setTheta}:}{Sets the value of the variance-component parameter.}
    \item{\code{getLower}:}{Retrieves the lower bound for \code{vcomp} and 
    \code{ccomp}, in that order.}
    \item{\code{getUpper}:}{Retrieves the upper bound for \code{vcomp} and 
    \code{ccomp}, in that order.}
    \item{\code{getLambda}:}{The relative covariance factor.}
    \item{\code{getLambdat.dp}:}{
    Returns a vector of length \code{nc} where each element contains the count 
    of nonzero entries in the corresponding column. Together with 
    \code{Lambdat.i}, defines the compressed sparse column representation 
    of \code{t(Lambda)}.
    }
    \item{\code{getLambdat.i}:}{
    Returns a vector of \code{0}-based row indices (\code{length <= nc*(n+1)/2}) 
    for nonzero entries in \code{t(Lambda)}, stored in column-major order. 
    Together with \code{Lambdat.dp}, defines the compressed sparse column 
    representation of \code{t(Lambda)}.
    }
    \item{\code{getVC}:}{Retrieves both \code{vcomp} and \code{ccomp}.}
    \item{\code{getVCNames}:}{Obtains the names of the values associated with
    \code{vcomp} and \code{ccomp} in terms of \code{cnm} and \code{gnm}.}
    \item{\code{setVC}:}{Sets the values for VC, which is composed of \code{vcomp}
    and \code{ccomp}.}
    \item{\code{getProfPar}:}{Converting the values from \code{vcomp}
    and \code{ccomp} (or \code{getVC}) to extract parameters from the object 
    and translates them to a different parameterization 
    suitable for profiling. If \code{profscale = "varcor"} then the 
    variance-covariance matrix is returned.}
    \item{\code{setProfPar}:}{Sets the parameters to use for profiling.}
    \item{\code{getProfLower}:}{Obtains the lower bounds of the values obtained
    by \code{getProfPar} for profiling.}
    \item{\code{getProfUpper}:}{Obtains the upper bounds of the values obtained
    by \code{getProfPar} for profiling.}
  }
}


\examples{
\dontshow{
.nms <- c("getPar", "getTheta", "getLambda", "getVC", "getReCovs")
.fns <- mget(.nms, envir = getNamespace("lme4"), mode = "function")
list2env(.fns, envir = environment())
}
## Unstructured
fm1.us <- lmer(Reaction ~ Days + us(Days | Subject), sleepstudy)
## Diagional
fm1.diag <- lmer(Reaction ~ Days + diag(Days | Subject), sleepstudy)
fm1.diag.hom <- lmer(Reaction ~ Days + diag(Days | Subject, hom = TRUE), 
  sleepstudy)
## Compound symmetry
fm1.cs <- lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy)
fm1.cs.hom <- lmer(Reaction ~ Days + cs(Days | Subject, hom = TRUE), 
  sleepstudy)
## Auto-regressive order 1
sleepstudy$Daysf <- factor(sleepstudy$Days, ordered = TRUE)
fm1.ar1 <- lmer(Reaction ~ Daysf + ar1(0 + Daysf | Subject, hom = TRUE), 
                sleepstudy, REML = FALSE)

## Applying some of the S4 methods for merMod
getPar(fm1.us)
getTheta(fm1.us)


## The same could be applied to covariance
fm1.cov <- getReCovs(fm1.us)[[1]]
getPar(fm1.cov)
getTheta(fm1.cov)
getLambda(fm1.cov)
getVC(fm1.cov)
\dontshow{
rm(list = c(.nms, ".nms", ".fns"))
}
}
\keyword{classes}
\keyword{internal}

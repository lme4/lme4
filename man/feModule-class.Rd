\name{feModule-class}
\Rdversion{1.1}
\docType{class}
\alias{feModule-class}
\alias{deFeMod-class}
\alias{spFeMod-class}
\title{Class "deFeMod"}
\description{
  The \code{"deFeMod"} and \code{"spFeMod"} classes are response module
  classes for the fixed-effects terms in a mixed model.  They inherit
  from the \code{"\linkS4class{predModule}"} classes
  \code{"\linkS4class{dPredModule}"} and
  \code{"\linkS4class{sPredModule}"}, respectively and contain
  additional slots related to the random-effects module.  These
  additional slots are inherited from the
  \code{"\linkS4class{feModule}"} class.
}
\section{Objects from the Class}{
  Objects from these classes are usually created as part of
  \code{\linkS4class{merMod}} object returned by functions
  \code{\link{lmer}}, \code{\link{glmer}} or \code{\link{nlmer}}.
}
\section{Slots}{
  \describe{
    \item{\code{RZX}:}{The off-diagonal part of the Cholesky
      decomposition for the penalized least squares problem - of class 
      \code{"\linkS4class{dgeMatrix}"} for \code{"deFeMod"} objects or
      \code{"\linkS4class{dgCMatrix}"} for \code{"spFeMod"} objects.}
    \item{\code{coef}, \code{X}, \code{fac}:}{Inherited from the
      appropriate \code{"\linkS4class{predModule}"} class.}
  }
}
\section{Methods}{
No methods defined with class "feModule" in the signature.
}
\seealso{\code{\linkS4class{merMod}}}
\examples{
showClass("deFeMod")
showClass("spFeMod")
}
\keyword{classes}

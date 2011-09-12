\name{lmerMod-class}
\Rdversion{1.1}
\docType{class}
\alias{lmerMod-class}

\title{Class \code{"lmerMod"} of Fitted Linear Mixed-Effect Models}

\description{
  A fitted linear mixed-effects model contains specific response module
  plus the parts, i.e., \code{\link{slot}}s of a general mixed-effects
  model, from class \code{\linkS4class{merMod}}, notably with it
  prediction module.
}
\section{Objects from the Class}{
  Objects can in priniciple be created by calls of the form
  \code{new("lmerMod", ...)}, but practically from \code{\link{lmer}()}.
}
\section{Slots}{
  \describe{
    \item{\code{resp}:}{a response module object of class
      \code{"\linkS4class{lmerResp}"}, see there.}

    \item{\code{call}, \code{frame}, ..., \code{pp}:}{all other slots
      are inherited from the super class \code{"\linkS4class{merMod}"}.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{merMod}"}, directly.
}
\section{Methods}{
  No methods defined with class "lmerMod" in the signature.
}
% \references{
% }
\seealso{
  \code{\link{lmer}()}; class \code{\linkS4class{merMod}}.
}
\examples{
showClass("lmerMod")
}
\keyword{classes}

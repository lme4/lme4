\name{lmerResp-class}
\Rdversion{1.1}
\docType{class}
\alias{lmerResp-class}
\title{"lmerResp" class}
\description{
  The \code{"lmerResp"} class is a trivial extension of the
  \code{\linkS4class{respModule}} class% from the \pkg{Matrix} package.
  It adds an integer scalar slot, \code{REML}, which is either
  zero, for maximum likelihood estimation, or \code{p}, the number of
  columns in the fixed-effects model matrix, for REML estimation.
}
\section{Objects from the Class}{
  Objects from these classes are usually created as part of
  \code{\linkS4class{merMod}} object returned by \code{\link{lmer}}.
}
\section{Slots}{
  \describe{
    \item{\code{REML}:}{Integer indicator of whether to estimate
      parameters according to the REML criterion - \code{lmerResp}
      objects only.}
  }
  All other slots are inherited from the \code{\linkS4class{respModule}}
  class of the \pkg{Matrix} package.
}
\section{Methods}{
No methods defined with class "merResp" in the signature.
}
\examples{
showClass("lmerResp")
}
\keyword{classes}

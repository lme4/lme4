\name{VarCorr-class}
\docType{class}
\alias{VarCorr-class}
\alias{show,VarCorr-method}
\title{Class "VarCorr"}
\description{The variance and correlation information for a
  random-effects model.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("VarCorr", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{scale}:}{The estimated standard deviation of the lowest
      level noise term in the model.}
    \item{\code{reSumry}:}{Object of class \code{"list"}, a list
      containing the summary of the positive-definte matrices for the
      random-effects levels.}
  }
}
\section{Methods}{
  \describe{
    \item{\link[methods:show]{show}}{\code{signature(object =
	"VarCorr")}: show the object.}
  }
}
\keyword{classes}

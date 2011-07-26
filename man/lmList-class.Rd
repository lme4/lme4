\name{lmList-class}
\docType{class}
\alias{lmList-class}
% not exported (FIXME?)
\alias{lmList.confint-class}
%
\alias{coef,lmList-method}
\alias{formula,lmList-method}
\alias{confint,lmList-method}
\alias{plot,lmList-method}
\alias{show,lmList-method}
\alias{update,lmList-method}
% (need both alias - depending on R version ?? : )
\alias{plot,lmList.confint-method}
\alias{plot,lmList.confint,ANY-method}
\title{Class "lmList" of  'lm'  Objects on Common Model}
\description{Class \code{"lmList"} is an S4 class with basically a list
  of objects of class \code{\link{lm}} with a common model.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lmList", ...)}
  or, more commonly, by a call to \code{\link{lmList}}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}, a list of
      objects of class \code{lm}}
    \item{\code{call}:}{Object of class \code{"call"}, the function call
      used to create the \code{lmList} object.}
    \item{\code{pool}:}{Object of class \code{"logical"}, if TRUE then
      calculate the pooled standard deviation estimate when displaying
      the object.}
  }
}
\section{Extends}{
Class \code{"list"}, from data part.
Class \code{"vector"}, by class "list".
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "lmList")}: extract the
      coefficients for the linear models.}
    \item{formula}{\code{signature(x = "lmList")}: extract the formula
      used when creating the lmList object.}
    \item{confint}{\code{signature(object = "lmList")}: confidence
      intervals for the fitted linear models.}
%% FIXME: these return class "lmList.confint" (non exported class !!)
%%       and this has a nice plot method
    \item{\link{show}}{\code{signature(object =
	"lmList")}: show the object.}
    \item{\link[stats]{update}}{\code{signature(object = "lmList")}:
      update the fit.}
  }
}
\keyword{classes}

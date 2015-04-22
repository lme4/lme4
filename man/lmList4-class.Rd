\docType{class}
\name{lmList4-class}
\alias{lmList4-class}
\alias{show,lmList4-method}
\title{Class "lmList4" of 'lm' Objects on Common Model}
\description{
  Class \code{"lmList4"} is an S4 class with basically a
  \code{\link{list}} of objects of class \code{\link{lm}} with a common
  model (but different data); see \code{\link{lmList}()} which returns
  these.

  Package \pkg{nlme}'s \code{\link[nlme]{lmList}()} returns objects of
  S3 class \code{"lmList"} and provides methods for them, on which our
  methods partly build.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("lmList4", ...)} or, more commonly, by a call to
  \code{\link{lmList}()}.
}
\section{Methods}{
  A dozen \code{\link{methods}} are provided.  Currently, S4
  methods for \code{\link{show}}, coercion (\code{\link{as}(.,.)}) and
  others inherited via \code{"list"},  and S3 methods for
  \code{\link{coef}}, \code{\link{confint}}, \code{\link{fitted}},
  \code{\link{fixef}}, \code{\link{formula}}, \code{\link{logLik}},
  \code{\link{pairs}}, \code{\link{plot}}, \code{\link{predict}},
  \code{\link{print}}, \code{\link{qqnorm}}, \code{\link{ranef}},
  \code{\link{residuals}}, \code{\link{sigma}}, \code{\link{summary}},
  and \code{\link{update}}.

  \describe{
    \item{sigma(object)}{returns the standard deviation
      \eqn{\hat{\sigma}}{\sigma-hat} (of the errors in the linear models), assuming a
      \emph{common} variance \eqn{\sigma^2} by pooling (even when
      \code{pool = FALSE} was used in the fit).}
  }
}
\seealso{\code{\link{lmList}}}
\examples{
if(getRversion() >= "3.2.0") {
  (mm <- methods(class = "lmList4"))
  ## The S3 ("not S4") ones :
  mm[!attr(mm,"info")[,"isS4"]]
}
## For more examples:  example(lmList)  i.e., ?lmList
}
\keyword{classes}


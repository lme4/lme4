\name{getME}
\alias{getME}
\title{Extract or Get Generalize Components from a Fitted Mixed Effects Model}
\description{
  Extract (or \dQuote{get}) \dQuote{components} -- in a generalized
  sense -- from a fitted mixed    effects model, i.e. (in this version
  of the package) from an \R object of class \code{"\linkS4class{mer}"}.

  The goal is to provide \dQuote{everything a user may want} from a
  fitted \code{"mer"} object \emph{as far} as it is not available by
  methods, such as \code{\link{fixef}}, \code{\link{ranef}},
  \code{\link{vcov}}, etc.
}
\usage{
getME(object, name)
}
\arguments{
  \item{object}{a fitted mixed-effects model of class
    \code{"\linkS4class{mer}"}, i.e. typically the result of
    \code{\link{lmer}()}, \code{\link{glmer}()} or \code{\link{nlmer}()}.
  }
  \item{name}{a character string specifying the name of the
    \dQuote{component}.  Note this may not be the name of
    \code{\link{slot}} of \code{object}. }
}
\value{
  Unspecified, as very much depending on the \code{\link{name}}.
}
\seealso{
%  \code{\link{getCall}()} (in \R >= 2.14.0; otherwise in \pkg{MatrixModels}).
% no \link{} -> no warning ...
  \code{getCall()} (in \R >= 2.14.0; otherwise in \pkg{MatrixModels}).

  More standard methods for \code{*mer()} objects, such as
  \code{\link{ranef}},
  \code{\link{fixef}},
  \code{\link{vcov}}, etc.
}
\examples{
## shows many methds you should consider *before* getME():
showMethods(class = "mer")

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
Z <- getME(fm1, "Z")
stopifnot(is(Z, "CsparseMatrix"),
          c(180,36) == dim(Z),
	  all.equal(fixef(fm1), getME(fm1, "beta"),
		    check.attr=FALSE, tol = 0))

## All that can be accessed [potentially ..]:
(nmME <- eval(formals(getME)$name))

\dontshow{
## internal consistency check ensuring that all work:
## "try(.)" because some are not yet implemented:
str(parts <- sapply(nmME, function(nm) try(getME(fm1, nm)),
                    simplify=FALSE))
}% dont..
}
\keyword{utilities}

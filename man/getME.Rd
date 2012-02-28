\name{getME}
\alias{getL}
\alias{getL,merMod-method}
\alias{getME}
\title{Extract or Get Generalize Components from a Fitted Mixed Effects Model}
\usage{
getME(object,
      name = c("X", "Z", "Zt", "u", "Gp", "L", "Lambda", "Lambdat",
               "Lind", "RX", "RZX", "beta", "theta", "REML", "n_rtrms",
               "is_REML", "devcomp", "offset"))
}
\arguments{
  \item{object}{a fitted mixed-effects model of class
  \code{"\linkS4class{merMod}"}, i.e. typically the result
  of \code{\link{lmer}()}, \code{\link{glmer}()} or
  \code{\link{nlmer}()}.}

  \item{name}{a character string specifying the name of the
  \dQuote{component}.  Possible values are:\cr \describe{
  \item{X}{fixed-effects model matrix}
  \item{Z}{random-effects model matrix} \item{Zt}{transpose
  of random-effects model matrix} \item{u}{conditional mode
  of the \dQuote{spherical} random effects variable}
  \item{Gp}{groups pointer vector.  A pointer to the
  beginning of each group of random effects corresponding
  to the random-effects terms.} \item{L}{sparse Cholesky
  factor of the penalized random-effects model.}
  \item{Lambda}{relative covariance factor of the random
  effects.} \item{Lambdat}{transpose of the relative
  covariance factor of the random effects.}
  \item{Lind}{index vector for inserting elements of
  \eqn{\theta}{theta} into the nonzeros of
  \eqn{\Lambda}{Lambda}} \item{RX}{Cholesky factor for the
  fixed-effects parameters} \item{RZX}{cross-term in the
  full Cholesky factor} \item{beta}{fixed-effects parameter
  estimates (same as the result of \code{\link{fixef}})}
  \item{n_rtrms}{number of random-effects terms}
  \item{is_REML}{same as the result of
  \code{\link{isREML}}} \item{devcomp}{a list consisting of
  a named numeric vector, \dQuote{cmp}, and a named integer
  vector, \dQuote{dims}, describing the fitted model}
  \item{offset}{model offset} }}
}
\value{
  Unspecified, as very much depending on the
  \code{\link{name}}.
}
\description{
  Extract (or \dQuote{get}) \dQuote{components} -- in a
  generalized sense -- from a fitted mixed-effects model,
  i.e. (in this version of the package) from an object of
  class \code{"\linkS4class{merMod}"}.
}
\details{
  The goal is to provide \dQuote{everything a user may
  want} from a fitted \code{"merMod"} object \emph{as far}
  as it is not available by methods, such as
  \code{\link{fixef}}, \code{\link{ranef}},
  \code{\link{vcov}}, etc.
}
\seealso{
  \code{\link{getCall}()}, More standard methods for mer
  objects, such as \code{\link{ranef}},
  \code{\link{fixef}}, \code{\link{vcov}}, etc.
}
\examples{
## shows many methods you should consider *before* getME():
methods(class = "merMod")

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
Z <- getME(fm1, "Z")
stopifnot(is(Z, "CsparseMatrix"),
          c(180,36) == dim(Z),
	  all.equal(fixef(fm1), getME(fm1, "beta"),
		    check.attr=FALSE, tol = 0))

## All that can be accessed [potentially ..]:
(nmME <- eval(formals(getME)$name))
}
\keyword{utilities}


\name{getME}
\alias{getL}
\alias{getL,merMod-method}
\alias{getME}
\title{Extract or Get Generalize Components from a Fitted Mixed Effects Model}
\usage{
  getME(object,
    name = c("X", "Z", "Zt", "Ztlist", "y", "mu", "u", "b", "Gp", "Tp", "L", "Lambda", "Lambdat", "Lind", "A", "RX", "RZX", "sigma", "flist", "beta", "theta", "ST", "REML", "is_REML", "n_rtrms", "n_rfacs", "devcomp", "offset", "lower"))
}
\arguments{
  \item{object}{a fitted mixed-effects model of class
  \code{"\linkS4class{merMod}"}, i.e. typically the result
  of \code{\link{lmer}()}, \code{\link{glmer}()} or
  \code{\link{nlmer}()}.}

  \item{name}{a character vector specifying the name(s) of
  the \dQuote{component}. If \code{length(name)}>1, a named
  list of components will be returned. Possible values
  are:\cr \describe{ \item{X}{fixed-effects model matrix}
  \item{Z}{random-effects model matrix} \item{Zt}{transpose
  of random-effects model matrix.  Note that the structure
  of \code{Zt} has changed since \code{lme4.0}; to get a
  backward-compatible structure, use
  \code{do.call(Matrix::rBind,getME(.,"Ztlist"))}}
  \item{Ztlist}{list of components of the transpose of the
  random-effects model matrix, separated by individual
  variance component} \item{y}{response vector}
  \item{mu}{conditional mean of the response}
  \item{u}{conditional mode of the \dQuote{spherical}
  random effects variable} \item{b}{conditional mode of the
  random effects variable} \item{Gp}{groups pointer vector.
  A pointer to the beginning of each group of random
  effects corresponding to the random-effects terms,
  beginning with 0 and including a final element giving the
  total number of random effects} \item{Tp}{theta pointer
  vector.  A pointer to the beginning of the theta
  sub-vectors corresponding to the random-effects terms,
  beginning with 0 and including a final element giving the
  total number of random effects} \item{L}{sparse Cholesky
  factor of the penalized random-effects model.}
  \item{Lambda}{relative covariance factor of the random
  effects.} \item{Lambdat}{transpose of the relative
  covariance factor of the random effects.}
  \item{Lind}{index vector for inserting elements of
  \eqn{\theta}{theta} into the nonzeros of
  \eqn{\Lambda}{Lambda}} \item{A}{Scaled sparse model
  matrix (class
  \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for
  the unit, orthogonal random effects, \eqn{U}, equal to
  \code{getME(.,"Zt") \%*\% getME(.,"Lambdat")}}
  \item{RX}{Cholesky factor for the fixed-effects
  parameters} \item{RZX}{cross-term in the full Cholesky
  factor} \item{sigma}{residual standard error}
  \item{flist}{a list of the grouping variables (factors)
  involved in the random effect terms}
  \item{beta}{fixed-effects parameter estimates (identical
  to the result of \code{\link{fixef}}, but without names)}
  \item{theta}{random-effects parameter estimates: these
  are parameterized as the relative Cholesky factors of
  each random effect term} \item{ST}{a list of matrices
  giving the relative Cholesky factors for each random
  effect term} \item{n_rtrms}{number of random-effects
  terms} \item{n_rfacs}{number of distinct random-effects
  grouping factors} \item{REML}{restricted maximum
  likelihood} \item{is_REML}{same as the result of
  \code{\link{isREML}}} \item{devcomp}{a list consisting of
  a named numeric vector, \dQuote{cmp}, and a named integer
  vector, \dQuote{dims}, describing the fitted model}
  \item{offset}{model offset} \item{lower}{lower bounds on
  model parameters (random effects parameters only)} }}
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
\examples{
## shows many methods you should consider *before* using getME():
methods(class = "merMod")

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
}\% dont..
}
\seealso{
  \code{\link{getCall}()}, More standard methods for mer
  objects, such as \code{\link{ranef}},
  \code{\link{fixef}}, \code{\link{vcov}}, etc.: see
  \code{methods(class="merMod")}
}
\keyword{utilities}


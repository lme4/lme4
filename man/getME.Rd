\name{getME}
\alias{getL}
\alias{getL,merMod-method}
\alias{getME}
\title{Extract or Get Generalized Components from a Fitted Mixed Effects Model}
\usage{
getME(object,
      name = c("X", "Z", "Zt", "Ztlist", "mmList", "y", "mu", "u", "b",
               "Gp", "Tp", "L", "Lambda", "Lambdat", "Lind",
               "A", "RX", "RZX", "sigma", "flist",
               "fixef", "beta", "theta", "ST", "REML", "is_REML",
               "n_rtrms", "n_rfacs", "N", "n", "p", "q",
               "cnms", "devcomp", "offset", "lower"))
}
\arguments{
  \item{object}{a fitted mixed-effects model of class
    \code{"\linkS4class{merMod}"}, i.e., typically the result of
    \code{\link{lmer}()}, \code{\link{glmer}()} or \code{\link{nlmer}()}.}

  \item{name}{a character vector specifying the name(s) of
    the \dQuote{component}.  If \code{length(name)}>1, a named
    list of components will be returned.  Possible values are:\cr
    \describe{
      \item{X}{fixed-effects model matrix}
      \item{Z}{random-effects model matrix}
      \item{Zt}{transpose
	of random-effects model matrix.  Note that the structure
	of \code{Zt} has changed since \code{lme4.0}; to get a
	backward-compatible structure, use
	\code{do.call(Matrix::rBind,getME(.,"Ztlist"))}}
      \item{Ztlist}{list of components of the transpose of the
	random-effects model matrix, separated by individual
	variance component}
      \item{mmList}{list of raw model matrices associated with random
	effects terms}
      \item{y}{response vector}
      \item{mu}{conditional mean of the response}
      \item{u}{conditional mode of the \dQuote{spherical}
	random effects variable}
      \item{b}{conditional mode of the
	random effects variable}
      \item{Gp}{groups pointer vector.
	A pointer to the beginning of each group of random
	effects corresponding to the random-effects terms,
	beginning with 0 and including a final element giving the
	total number of random effects}
      \item{Tp}{theta pointer vector.  A pointer to the beginning of the theta
	sub-vectors corresponding to the random-effects terms,
	beginning with 0 and including a final element giving the
	number of thetas.}
      \item{L}{sparse Cholesky factor of the penalized random-effects model.}
      \item{Lambda}{relative covariance factor \eqn{\Lambda}{Lambda} of the random effects.}
      \item{Lambdat}{transpose \eqn{\Lambda'}{Lambda'} of \eqn{\Lambda}{Lambda} above.}
      \item{Lind}{index vector for inserting elements of
	\eqn{\theta}{theta} into the nonzeros of \eqn{\Lambda}{Lambda}.}
      \item{A}{Scaled sparse model matrix (class
	\code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for
	the unit, orthogonal random effects, \eqn{U}, equal to
	\code{getME(.,"Zt") \%*\% getME(.,"Lambdat")}}
      \item{RX}{Cholesky factor for the fixed-effects parameters}
      \item{RZX}{cross-term in the full Cholesky factor}
      \item{sigma}{residual standard error; note that \code{sigma(object)} is preferred.}
      \item{flist}{a list of the grouping variables (factors)
	involved in the random effect terms}
      \item{fixef}{fixed-effects parameter estimates}
      \item{beta}{fixed-effects parameter estimates (identical
	to the result of \code{\link{fixef}}, but without names)}
      \item{theta}{random-effects parameter estimates: these
	are parameterized as the relative Cholesky factors of
	each random effect term}
      \item{ST}{A list of S and T factors in the TSST' Cholesky
      factorization of the relative variance matrices of the random
      effects associated with each random-effects term.  The unit lower
      triangular matrix, \eqn{T}, and the diagonal matrix, \eqn{S}, for
      each term are stored as a single matrix with diagonal elements
      from \eqn{S} and off-diagonal elements from \eqn{T}.}
      \item{n_rtrms}{number of random-effects terms}
      \item{n_rfacs}{number of distinct random-effects grouping factors}
      \item{N}{number of rows of \code{X}}
      \item{n}{length of \code{y}}
      \item{p}{number of columns of \code{X}}
      \item{q}{number of columns of \code{Z}}
      \item{cnms}{the \dQuote{component names}, a \code{\link{list}}.}
      \item{REML}{\code{0} indicates the model was fitted by maximum
	likelihood, any other positive integer indicates fitting by 
	restricted maximum likelihood}
      \item{is_REML}{same as the result of \code{\link{isREML}(.)}}
      \item{devcomp}{a list consisting of a named numeric vector,
        \code{cmp}, and a named integer vector, \code{dims}, describing
        the fitted model.  The elements of \code{cmp} are:\cr
	\describe{
          \item{ldL2}{twice the log determinant of \code{L}}
          \item{ldRX2}{twice the log determinant of \code{RX}}
          \item{wrss}{weighted residual sum of squares}
          \item{ussq}{squared length of \code{u}}
          \item{pwrss}{penalized weighted residual sum of squares,
            \dQuote{wrss + ussq}}
	  \item{drsum}{sum of residual deviance (GLMMs only)}
	  \item{REML}{REML criterion at optimum (LMMs fit
            by REML only)}
	  \item{dev}{deviance criterion at optimum
            (models fit by ML only)}
	  \item{sigmaML}{ML estimate of residual standard deviation}
	  \item{sigmaREML}{REML estimate of residual standard deviation}
	  \item{tolPwrss}{tolerance for declaring convergence in the
	  penalized iteratively weighted residual sum-of-squares (GLMMs only)}
	} The elements of \code{dims} are:\cr
	\describe{
	  \item{N}{number of rows of \code{X}}
          \item{n}{length of \code{y}}
	  \item{p}{number of columns of \code{X}}
	  \item{nmp}{\code{n-p}}
	  \item{nth}{length of \code{theta}}
	  \item{q}{number of columns of \code{Z}}
	  \item{nAGQ}{see \code{\link{glmer}}}
          \item{compDev}{see \code{\link{glmerControl}}}
          \item{useSc}{\code{TRUE} if model has a scale parameter}
	  \item{reTrms}{number of random effects terms}
          \item{REML}{\code{0} indicates the model was fitted by maximum
	    likelihood, any other positive integer indicates fitting by 
	    restricted maximum likelihood}
          \item{GLMM}{\code{TRUE} if a GLMM}
	  \item{NLMM}{\code{TRUE} if an NLMM}
	}
      }
      \item{offset}{model offset}
      \item{lower}{lower bounds on model parameters (random effects
	parameters only).}
    }
  }
}
\value{
  Unspecified, as very much depending on the \code{\link{name}}.
}
\description{
  Extract (or \dQuote{get}) \dQuote{components} -- in a
  generalized sense -- from a fitted mixed-effects model,
  i.e., (in this version of the package) from an object of
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
  \code{\link{getCall}()}.  More standard methods for \code{"merMod"}
  objects, such as \code{\link{ranef}}, \code{\link{fixef}},
  \code{\link{vcov}}, etc.: see \code{methods(class="merMod")}
}
\examples{
## shows many methods you should consider *before* using getME():
methods(class = "merMod")

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
Z <- getME(fm1, "Z")
stopifnot(is(Z, "CsparseMatrix"),
          c(180,36) == dim(Z),
	  all.equal(fixef(fm1), getME(fm1, "beta"),
		    check.attributes=FALSE, tolerance = 0))

## A way to get *all* getME()s :
getME.all <- function(obj) {
    sapply(eval(formals(getME)$name), getME, object=obj, simplify=FALSE)
}
## internal consistency check ensuring that all work:
parts <- getME.all(fm1)
str(parts, max=2)
}
\keyword{utilities}


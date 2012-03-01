\name{bootMer}
\alias{bootMer}
\title{Model-based (Semi-)Parametric Bootstrap for Mixed Models}
\usage{
  bootMer(x, FUN, nsim = 1, seed = NULL, use.u = FALSE,
    type = c("parametric", "semiparametric"),
    verbose = FALSE, control = list(), .progress = "none",
    PBargs = list())
}
\arguments{
  \item{x}{fitted \code{*lmer()} model, see
  \code{\link{lmer}}, \code{\link{glmer}}, etc.}

  \item{FUN}{a \code{\link{function}(x)}, computing the
  \emph{statistic} of interest, which must be a numeric
  vector, possibly named.}

  \item{nsim}{number of simulations, positive integer; the
  bootstrap \eqn{B} (or \eqn{R}).}

  \item{seed}{optional argument to \code{\link{set.seed}}.}

  \item{use.u}{logical, indicating, if the spherized random
  effects should be simulated / bootstrapped as well.  If
  \code{FALSE}, they are not changed, and all inference is
  conditional on these.}

  \item{type}{character string specifying the type of
  bootstrap, \code{"parametric"} or
  \code{"semiparametric"}.  Default is \code{"parametric"}.
  Partial matching is allowed.}

  \item{verbose}{logical indicating if progress should
  print output}

  \item{control}{an optional \code{\link{list}}, to be
  passed to the minimizer (of the log-likelihood, or RE
  likelihood).}

  \item{.progress}{character string - type of progress bar
  to display.  Default is \code{"none"}.}

  \item{PBargs}{a list of additional arguments to the
  progress bar function.}
}
\value{
  an object of S3 \code{\link{class}} \code{"boot"},
  compatible with \pkg{boot} package's \code{boot()}
  result.
}
\description{
  Perform model-based (Semi-)parametric bootstrap for mixed
  models.
}
\details{
  The semi-parametric variant is not yet implemented, and
  we only provide a method for \code{\link{lmer}} and
  \code{\link{glmer}} results.

  The working name for bootMer() was
  \dQuote{simulestimate()}, as it is an extension of
  \code{\link{simulate}}, but we rather want to emphasize
  its potential for valid inference.

  In each of the \code{nsim} simulations --- that is what
  the \emph{parametric} bootstrap does, both
  \dQuote{\emph{spherized}} random effects \eqn{u} and the
  i.i.d errors \eqn{\epsilon} are generated, using
  \code{\link{rnorm}()} with parameters corresponding to
  the fitted model \code{x}.
}
\examples{
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
## see ?"profile-methods"
mySumm <- function(.) { s <- sigma(.)
    c(beta =getME(., "beta"), sigma = s, sig01 = s * getME(., "theta")) }
(t0 <- mySumm(fm01ML)) # just three parameters

## 3.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
system.time( boo01 <- bootMer(fm01ML, mySumm, nsim = 100) )

## to "look" at it
if(need.boot <- is.na(match("package:boot", search())))
  require("boot")## is a recommended package, i.e. *must* be there
boo01 # look at estimated bias for sig01 (-9.1, also when nsim = 1000)

## ------ Bootstrap-based confidence intervals ------------

(bCI.1 <- boot.ci(boo01, index=1, type=c("norm", "basic", "perc")))# beta

## Sigma - original scale, and interval calculation on log-scale:
(bCI.2  <- boot.ci(boo01, index=2, type=c("norm", "basic", "perc")))
(bCI.2l <- boot.ci(boo01, index=2, type=c("norm", "basic", "perc"),
                   h = log, hdot = function(.) 1/., hinv = exp))

(bCI.3 <- boot.ci(boo01, index=3, type=c("norm", "basic", "perc")))# sig01

}
\references{
  Davison, A.C. and Hinkley, D.V. (1997) \emph{Bootstrap
  Methods and Their Application}.  Cambridge University
  Press.
}
\seealso{
  For inference, including confidence intervals,
  \code{\link{profile-methods}}.

  \code{\link[boot]{boot}()}, and then
  \code{\link[boot]{boot.ci}} from package \pkg{boot}.
}
\keyword{htest}
\keyword{models}

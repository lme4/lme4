\name{bootMer}
\alias{bootMer}
\title{Model-based (Semi-)Parametric Bootstrap for Mixed Models}
\usage{
bootMer(x, FUN, nsim = 1, seed = NULL, use.u = FALSE,
	type = c("parametric", "semiparametric"),
	verbose = FALSE, .progress = "none", PBargs = list(),
	parallel = c("no", "multicore", "snow"),
	ncpus = getOption("boot.ncpus", 1L), cl = NULL)
}
\arguments{
  \item{x}{a fitted \code{merMod} object: see
    \code{\link{lmer}}, \code{\link{glmer}}, etc.}
  \item{FUN}{a function taking a fitted
    \code{merMod} object as input and returning the
    \emph{statistic} of interest, which must be a (possibly named)
    numeric vector.}
  \item{nsim}{number of simulations, positive integer; the
    bootstrap \eqn{B} (or \eqn{R}).}
  \item{seed}{optional argument to \code{\link{set.seed}}.}
  \item{use.u}{logical, indicating whether the spherical
    random effects should be simulated / bootstrapped as
    well.  If \code{TRUE}, they are not changed, and all
    inference is conditional on these values. If
    \code{FALSE}, new normal deviates are drawn (see
    Details).}
  \item{type}{character string specifying the type of
    bootstrap, \code{"parametric"} or
    \code{"semiparametric"}; partial matching is allowed.}
  \item{verbose}{logical indicating if progress should
    print output}
  \item{.progress}{character string - type of progress bar
    to display.  Default is \code{"none"}; the function will
    look for a relevant \code{*ProgressBar} function, so
    \code{"txt"} will work in general; \code{"tk"} is
    available if the \pkg{tcltk} package is loaded; or
    \code{"win"} on Windows systems. Progress bars are
    disabled (with a message) for parallel operation.}
  \item{PBargs}{a list of additional arguments to the
    progress bar function (the package authors like
    \code{list(style=3)}).}
  \item{parallel}{The type of parallel operation to be used (if any).
    If missing, the
    default is taken from the option \code{"boot.parallel"} (and if that
    is not set, \code{"no"}).}
  \item{ncpus}{integer: number of processes to be used in parallel operation:
    typically one would choose this to be the number of available CPUs.}
  \item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
    \code{parallel = "snow"}.  If not supplied, a cluster on the
    local machine is created for the duration of the \code{boot} call.}
}
\value{
  an object of S3 \code{\link{class}} \code{"boot"},
  compatible with \pkg{boot} package's
  \code{\link[boot]{boot}()} result.
}
\description{
  Perform model-based (Semi-)parametric bootstrap for mixed
  models.
}
\details{
  The semi-parametric variant is only partially implemented, and
  we only provide a method for \code{\link{lmer}} and
  \code{\link{glmer}} results.

  The working name for bootMer() was
  \dQuote{simulestimate()}, as it is an extension of \code{simulate}
  (see \link{simulate.merMod}), but we want to emphasize its potential
  for valid inference.
  \itemize{
    \item If \code{use.u} is \code{FALSE} and \code{type} is
    \code{"parametric"}, each simulation generates new values of both
    the \dQuote{\emph{spherical}} random effects \eqn{u} and the
    i.i.d. errors \eqn{\epsilon}, using \code{\link{rnorm}()}
    with parameters corresponding to the fitted model \code{x}.

    \item If \code{use.u} is \code{TRUE} and \code{type=="parametric"},
    only the i.i.d. errors (or, for GLMMs, response values drawn from
    the appropriate distributions) are resampled, with the values of
    \eqn{u} staying fixed at their estimated values.

    \item If \code{use.u} is \code{TRUE} and \code{type=="semiparametric"},
    the i.i.d. errors are sampled from the distribution of (response)
    residuals.  (For GLMMs, the resulting
    sample will no longer have the same properties as the original
    sample, and the method may not make sense; a warning is generated.)
    The semiparametric bootstrap is currently an experimental feature,
    and therefore may not be stable.

    \item The case where \code{use.u} is \code{FALSE} and
    \code{type=="semiparametric"} is not implemented; Morris (2002)
    suggests that resampling from the estimated values of \eqn{u} is not
    good practice.
  } %% itemize
} %% details
\references{
  Davison, A.C. and Hinkley, D.V. (1997)
  \emph{Bootstrap Methods and Their Application}.
  Cambridge University Press.

  Morris, J. S. (2002).
  The BLUPs Are Not \sQuote{best} When It Comes to Bootstrapping.
  \emph{Statistics & Probability Letters} \bold{56}(4): 425--430.
  doi:10.1016/S0167-7152(02)00041-X.
}
\seealso{
  \itemize{
    \item \code{\link{confint.merMod}},
    for a more specific approach to bootstrap confidence
    intervals on parameters.
    \item \code{\link{refit}()}, or \code{\link[pbkrtest]{PBmodcomp}()}
    from the \pkg{pbkrtest} package, for parametric bootstrap comparison
    of models.
    \item \code{\link[boot]{boot}()}, and then
    \code{\link[boot]{boot.ci}}, from the \pkg{boot} package.
    \item \code{\link{profile-methods}}, for likelihood-based inference,
    including confidence intervals.
    \item \code{\link{pvalues}},
    for more general approaches to inference and p-value computation
    in mixed models.
  }
}
\examples{
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
## see ?"profile-methods"
mySumm <- function(.) { s <- sigma(.)
    c(beta =getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) }
(t0 <- mySumm(fm01ML)) # just three parameters
## alternatively:
mySumm2 <- function(.) {
   c(beta=fixef(.),sigma=sigma(.), sig01=sqrt(unlist(VarCorr(.))))
}

set.seed(101)
## 3.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
system.time( boo01 <- bootMer(fm01ML, mySumm, nsim = 100) )

## to "look" at it
require("boot") ## a recommended package, i.e. *must* be there
boo01
## note large estimated bias for sig01
## (~30\% low, decreases _slightly_ for nsim = 1000)

## extract the bootstrapped values as a data frame ...
head(as.data.frame(boo01))

## ------ Bootstrap-based confidence intervals ------------

## intercept
(bCI.1 <- boot.ci(boo01, index=1, type=c("norm", "basic", "perc")))# beta

## Residual standard deviation - original scale:
(bCI.2  <- boot.ci(boo01, index=2, type=c("norm", "basic", "perc")))
## Residual SD - transform to log scale:
(bCI.2l <- boot.ci(boo01, index=2, type=c("norm", "basic", "perc"),
                   h = log, hdot = function(.) 1/., hinv = exp))

## Among-batch variance:
(bCI.3 <- boot.ci(boo01, index=3, type=c("norm", "basic", "perc")))# sig01

## Graphical examination:
plot(boo01,index=3)

## Check stored values from a longer (1000-replicate) run:
load(system.file("testdata","boo01L.RData",package="lme4"))
plot(boo01L,index=3)
}
\keyword{htest}
\keyword{models}


\name{prt-utilities}
\title{Print and Summary Method Utilities for Mixed Effects}
\alias{.prt.methTit}
\alias{.prt.VC}
\alias{.prt.aictab}
\alias{.prt.call}
\alias{.prt.family}
\alias{.prt.grps}
\alias{.prt.methTit}
\alias{.prt.resids}
\alias{.prt.warn}
\alias{formatVC}
\alias{llikAIC}
\alias{methTitle}
\description{
  The \code{\link{print}}, \code{\link{summary}} methods (including the
  \code{print} for the \code{summary()} result) in \pkg{lme4} are
  modular, using about ten small utility functions.  Other packages,
  building on \pkg{lme4} can use the same utilities for ease of
  programming and consistency of output.

  Notably see the Examples.

  \code{llikAIC()} extracts the  log likelihood, AIC, and related
  statics from a Fitted LMM.

  \code{formatVC()} \dQuote{format()}s the \code{\link{VarCorr}} matrix of the
  random effects -- for \code{\link{print}()}ing and
  \code{\link{show}()}ing; it is also the \dQuote{workhorse} of
  \code{.prt.VC()}, and returns a \code{\link{character}} matrix.

  \code{.prt.*()} all use \code{\link{cat}} and \code{\link{print}} to
  produce output.
}
\usage{
llikAIC(object, devianceFUN = devCrit, chkREML = TRUE,
        devcomp = object@devcomp)

methTitle(dims)

.prt.methTit(mtit, class)
.prt.family (famL)
.prt.resids (resids, digits, title = "Scaled residuals:", \dots)
.prt.call   (call, long = TRUE)
.prt.aictab (aictab, digits = 1)
.prt.grps   (ngrps, nobs)
.prt.warn   (optinfo, summary = FALSE, \dots)

.prt.VC (varcor, digits, comp, formatter = format, \dots)
formatVC(varcor, digits = max(3, getOption("digits") - 2),
         comp = "Std.Dev.", formatter = format,
         useScale = attr(varcor, "useSc"), \dots)
}
\arguments{
  %% llikAIC() :
  \item{object}{a LMM model fit}
  \item{devianceFUN}{the function to be used for computing the deviance;
    should not be changed for \pkg{lme4} created objects.}
  \item{chkREML}{optional logical indicating if \code{object} maybe a REML
    fit.}% use TRUE for \pkg{lme4} fits
  \item{devcomp}{for \pkg{lme4} always the equivalent of
    \code{object@devcomp}; here a \code{\link{list}}}%... FIXME

  %% methTitle():
  \item{dims}{for \pkg{lme4} always the equivalent of
    \code{object@devcomp$dims}, a named vector or list with components
    \code{"GLMM"}, \code{"NLMM"}, \code{"REML"}, and \code{"nAGQ"} of
    which the first two are \code{\link{logical}} scalars, and the latter
    two typically are \code{FALSE} or \code{\link{numeric}}.}

  %% .prt.methTit
  \item{mtit}{the result of \code{methTitle(object)}}
  \item{class}{typically \code{\link{class}(object)}.}

  %% .prt.family (famL)
  \item{famL}{a \code{\link{list}} with components \code{family} and
    \code{link}, each a \code{\link{character}} string; note that standard
    \R \code{\link{family}} objects can be used directly, as well.}

  %% .prt.resids (resids, digits, title = "Scaled residuals:", \dots)
  \item{resids}{numeric vector of model \code{\link{residuals}}.}
  \item{digits}{non-negative integer of (significant) digits to print minimally.}
  \item{title}{\code{\link{character}} string.}
  \item{\dots}{optional arguments passed on, e.g., to \code{\link{residuals}()}.}

  %% .prt.call   (call, long = TRUE)
  \item{call}{the \code{\link{call}} of the model fit; e.g., available
    via (generic) function \code{\link{getCall}()}.}
  \item{long}{logical indicating if the output may be long, e.g.,
    printing the \code{control} part of the call if there is one.}

  %% .prt.aictab (aictab, digits = 1)
  \item{aictab}{typically the \code{AICtab} component of the result of
    \code{llikAIC()}.}

  %% .prt.VC     (varcor, digits, comp, formatter = format, \dots)
  \item{varcor}{typically the result of \code{\link{VarCorr}()}.}
  \item{comp}{optional ...}% ... FIXME
  \item{formatter}{a \code{\link{function}} used for formatting the numbers.}

  %% .prt.grps   (ngrps, nobs)
  \item{ngrps}{integer (vector), typically the result of
    \code{\link{ngrps}(object)}.}
  \item{nobs}{integer; the number of observations, e.g., the result
    of \code{\link{nobs}}.}

  %% .prt.warn   (optinfo, summary = FALSE, \dots)
  \item{optinfo}{typically \code{object @ optinfo}, the optimization
    infos, including warnings if there were.}
  \item{summary}{logical}

  %% formatVC
  \item{useScale}{(logical) whether the parent model estimates a scale
    parameter}
}
\value{
  \code{llikAIC()} returns a \code{\link{list}} with components

  \item{logLik}{which is \code{\link{logLik}(object)}, and}
  \item{AICtab}{ a \dQuote{table} of \code{\link{AIC}}, \code{\link{BIC}},
    \code{\link{logLik}}, deviance and \code{\link{df.residual}()} values.}
}
\examples{
## Create a few "lme4 standard" models ------------------------------
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fmM <- update(fm1, REML=FALSE) # -> Maximum Likelihood

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
gmA <- update(gm1, nAGQ = 5)


(lA1 <- llikAIC(fm1))
(lAM <- llikAIC(fmM))
(lAg <- llikAIC(gmA))

(m1 <- methTitle(fm1 @ devcomp $ dims))
(mM <- methTitle(fmM @ devcomp $ dims))
(mG <- methTitle(gm1 @ devcomp $ dims))
(mA <- methTitle(gmA @ devcomp $ dims))

.prt.methTit(m1, class(fm1))
.prt.methTit(mA, class(gmA))

.prt.family(gaussian())
.prt.family(binomial())
.prt.family( poisson())

.prt.resids(residuals(fm1), digits = 4)
.prt.resids(residuals(fmM), digits = 2)

.prt.call(getCall(fm1))
.prt.call(getCall(gm1))

.prt.aictab ( lA1 $ AICtab ) # REML
.prt.aictab ( lAM $ AICtab ) # ML --> AIC, BIC, ...

V1 <- VarCorr(fm1)
m <- formatVC(V1)
stopifnot(is.matrix(m), is.character(m), ncol(m) == 4)
print(m, quote = FALSE) ## prints all but the first line of .prt.VC() below:
.prt.VC( V1, digits = 4, formatter = format)
## Random effects:
##  Groups   Name        Std.Dev. Corr
##  Subject  (Intercept) 24.740
##           Days         5.922   0.07
##  Residual             25.592

p1 <- capture.output(V1)
p2 <- capture.output( print(m, quote=FALSE) )
pX <- capture.output( .prt.VC(V1, digits = max(3, getOption("digits")-2)) )
stopifnot(identical(p1, p2),
          identical(p1, pX[-1])) # [-1] : dropping 1st line

.prt.grps(ngrps = ngrps(fm1),
          nobs  = nobs (fm1))
## --> Number of obs: 180, groups:  Subject, 18

.prt.warn(fm1 @ optinfo) # nothing .. had no warnings
}
\keyword{utilities}

\name{plots.thpr}
\title{Mixed-Effects Profile Plots (Regular / Density / Pairs)}
\alias{xyplot.thpr}
\alias{densityplot.thpr}
\alias{splom.thpr}
\description{
  Xyplot, Densityplot, and Pairs plot methods for a mixed-effects model
  profile.

  \code{xyplot()} draws \dQuote{zeta diagrams}, also visualizing
  confidence intervals and their asymmetry.

  \code{densityplot()} draws the profile densities.

  \code{splom()} draws profile pairs plots.  Contours are for the
  marginal two-dimensional regions (i.e. using df = 2).
}
\usage{
\method{xyplot}{thpr}(x, data = NULL,
       levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), 1)),
       conf = c(50, 80, 90, 95, 99)/100,
       absVal = FALSE, scales=NULL,
       which=1:nptot, ...)

\method{densityplot}{thpr}(x, data, ...)

\method{splom}{thpr}(x, data,
      levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), 2)),
      conf = c(50, 80, 90, 95, 99)/100, which=1:nptot,
      draw.lower=TRUE,
     draw.upper=TRUE,
     ...)
}
\arguments{
  \item{x}{a mixed-effects profile, i.e., of class \code{"thpr"},
    typically resulting from \code{\link{profile}(fm)} where \code{fm}
    is a fitted model from \code{\link{lmer}} (or its generalizations).}
  \item{data}{unused - only for compatibility with generic.}
  \item{levels}{the contour levels to be shown; usually
    derived from \code{conf}.}
  \item{conf}{numeric vector of confidence levels to be
    shown as contours.}
  \item{absVal}{logical indicating if \code{\link{abs}(.)}olute values
      should be plotted, often preferred for confidence interval
      visualization.}
    \item{scales}{plotting options to be passed to \code{\link{xyplot}}}
    \item{which}{integer or character vector indicating which parameters
    to profile: default is all parameters (see
    \code{\link{profile-methods}}
    for details).}
  \item{draw.lower}{(logical) draw lower-triangle (zeta scale) panels?}
  \item{draw.upper}{(logical) draw upper-triangle (standard dev/cor scale) panels?}
  \item{...}{further arguments passed to \code{\link{xyplot}},
    \code{\link{densityplot}}, or \code{\link{splom}} from package
    \pkg{lattice}, respectively.}
}
\seealso{ \code{\link[=profile-methods]{profile}}, notably for an
  example.
}
\value{
  \describe{
    \item{xyplot:}{a density plot, a \code{"trellis"} object (\pkg{lattice} package)
      which when \code{\link{print}()}ed produces plots on the current
      graphic device.}
    \item{densityplot:}{a density plot, a \code{"trellis"} object, see above.}
    \item{splom:}{a pairs plot, aka \bold{s}catter\bold{plo}t
      \bold{m}atrix, a \code{"trellis"} object, see above.}
  }
}
\examples{
## see   example("profile.merMod")
}

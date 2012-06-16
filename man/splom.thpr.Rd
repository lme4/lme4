\name{splom.thpr}
\alias{splom.thpr}
\title{Profile pairs plot}
\usage{
  \method{splom}{thpr} (x, data,
    levels = sqrt(qchisq(pmax.int(0, pmin.int(1, conf)), 2)),
    conf = c(50, 80, 90, 95, 99)/100, ...)
}
\arguments{
  \item{x}{the result of \code{\link{profile}} (or very
  similar structure)}

  \item{data}{unused - only for compatibility with generic}

  \item{levels}{the contour levels to be shown; usually
  derived from \code{conf}}

  \item{conf}{numeric vector of confidence levels to be
  shown as contours}

  \item{...}{further arguments passed to
  \code{\link{splom}}}
}
\description{
  Draws profile pairs plots.  Contours are for the marginal
  two-dimensional regions (i.e. using df = 2).
}


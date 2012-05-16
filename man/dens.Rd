\name{dens}
\alias{dens}
\title{Approximate densities from profiles}
\usage{
  dens(pr, npts = 201, upper = 0.999)
}
\arguments{
  \item{pr}{a profile object}

  \item{npts}{number of points at which to evaluate the
  density}

  \item{upper}{upper bound on cumulative for a cutoff}
}
\value{
  a data frame
}
\description{
  Create an approximating density from a profile object
}


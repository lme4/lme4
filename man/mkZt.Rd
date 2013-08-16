\name{mkZt}
\alias{mkZt}
\title{Generate (transpose of) random effect design matrix}
\usage{
  mkZt(ff, bar, fr)
}
\arguments{
  \item{ff}{the grouping factor}

  \item{bar}{the random effect specification as returned by
  \code{findbars}}

  \item{fr}{a \code{model.frame}}
}
\value{
  a list with entries \code{Zt}, \code{nc} for the number
  of effects per level of the grouping factor, and
  \code{cnms} giving the names of the different random
  effects.
}
\description{
  Generate (transpose of) random effect design matrix
}


\name{glmer.nb}
\alias{glmer.nb}
\title{glmer() for Negative Binomial}
\usage{
  glmer.nb(..., interval = log(th) + c(-3, 3),
    verbose = FALSE)
}
\arguments{
  \item{...}{formula, data, etc: the arguments for
  \code{\link{glmer}(..)} (apart from \code{family}!).}

  \item{interval}{interval in which to start the
  optimization}

  \item{verbose}{logical indicating how much progress
  information should be printed.}
}
\description{
  glmer() for Negative Binomial
}


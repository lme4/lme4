\name{mcmcsamp}
\alias{mcmcsamp}
\title{Create an McMC sample}
\usage{
  mcmcsamp(object, n = 1L, verbose = FALSE, ...)
}
\arguments{
  \item{object}{a fitted model object}

  \item{n}{number of samples to generate.  Defaults to 1
  but that value doesn't make sense.}

  \item{verbose}{should verbose output be given?}

  \item{...}{}
}
\value{
  a Markov chain Monte Carlo sample as a matrix
}
\description{
  Create a Markov chain Monte Carlo sample from the
  posterior distribution of the parameters
}


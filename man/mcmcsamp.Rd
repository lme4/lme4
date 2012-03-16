\name{mcmcsamp}
\alias{mcmcsamp}
\title{Create an MCMC sample}
\usage{
  mcmcsamp(object, n = 1L, verbose = FALSE, ...)
}
\arguments{
  \item{object}{a fitted model object}

  \item{n}{number of samples to generate.  Defaults to 1;
  for real use values of 200-1000 are more typical}

  \item{verbose}{should verbose output be given?}

  \item{...}{additional, optional arguments (not used)}
}
\value{
  a Markov chain Monte Carlo sample as a matrix
}
\description{
  Create a Markov chain Monte Carlo sample from the
  posterior distribution of the parameters
}


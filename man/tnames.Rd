\name{tnames}
\alias{tnames}
\title{Construct names of individual theta/sd:cor components}
\usage{
  tnames(object, diag.only = FALSE, old = TRUE,
    prefix = NULL)
}
\arguments{
  \item{object}{a fixed model}

  \item{diag.only}{include only diagonal elements?}

  \item{old}{(logical) give backward-compatible results?}

  \item{prefix}{a character vector with two elements giving
  the prefix for diagonal (e.g. "sd") and off-diagonal
  (e.g. "cor") elements}
}
\description{
  Construct names of individual theta/sd:cor components
}


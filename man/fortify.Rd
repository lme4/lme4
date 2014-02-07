\name{fortify}
\alias{fortify}
\alias{fortify.merMod}
\title{add information to data based on a fitted model}
\usage{
  fortify.merMod(model, data = getData(model),
    ...)
}
\arguments{
  \item{model}{fitted model}
  \item{data}{original data set, if needed}
  \item{\dots}{additional arguments}
}
\description{
  add information to data based on a fitted model
}
\details{
  \code{fortify} is a function defined in the
  \pkg{ggplot2} package, q.v. for more details.
  \code{fortify} is \emph{not} defined here,
  and \code{fortify.merMod} is defined as a function
  rather than an S3 method, to avoid (1) inducing
  a dependency on \pkg{ggplot2} or (2) masking
  methods from \pkg{ggplot2}.
  This is currently an experimental feature.
}


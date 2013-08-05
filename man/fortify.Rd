\name{fortify}
\alias{fortify}
\alias{fortify.lmerMod}
\title{add information to data based on a fitted model}
\usage{
  fortify(model, data, ...)

  \method{fortify}{lmerMod} (model, data = getData(model),
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
  \code{ggplot2} package, q.v. for more details; the S3
  generic is just defined here to avoid inducing an
  additional \code{Imports:} dependency. This is currently
  an experimental feature.
}


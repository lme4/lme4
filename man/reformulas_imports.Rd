\name{reformulas_imports}
\title{(deprecated) formula-processing functions}
\description{
  Formula-processing functions that have moved to \code{reformulas}
}
\alias{subbars}
\alias{nobars}
\alias{findbars}
\alias{mkReTrms}
\usage{
subbars(...)

nobars(...)

findbars(...)

mkReTrms(...)

}

\arguments{
  \item{...}{arguments to be passed through to functions in the
    \code{reformulas} package}
}
\details{
  All of these functions are now located in the \code{reformulas}
  package. The versions here are 'shims' that throw a warning and
  call the \code{reformulas} version.
}

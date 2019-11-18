\name{fortify}
\alias{fortify}
\alias{fortify.merMod}
\alias{getData}
\alias{getData.merMod}
\title{add information to data based on a fitted model}
\usage{
  fortify.merMod(model, data = getData(model),
    ...)
\method{getData}{merMod}(object)
}
\arguments{
  \item{model}{fitted model}
  \item{object}{fitted model}
  \item{data}{original data set, if needed}
  \item{\dots}{additional arguments}
}
\description{
  \code{fortify} adds information to data based on a fitted model;
  \code{getData} retrieves data as specified in the \code{data}
  argument
}
\details{
  \itemize{
    \item{\code{fortify} is defined in the
    \CRANpkg{ggplot2} package, q.v. for more details.
    \code{fortify} is \emph{not} defined here,
    and \code{fortify.merMod} is defined as a function
    rather than an S3 method, to avoid (1) inducing
    a dependency on \pkg{ggplot2} or (2) masking
    methods from \pkg{ggplot2}.
    This feature is both experimental and semi-deprecated,
    as the help page for \code{fortify} itself says:
    \dQuote{Rather than using this function, I now recommend using the
      \code{broom} package, which implements a much wider range of
      methods. \code{fortify} may be deprecated in the future.}
    The \code{broom.mixed} package is recommended for mixed models
    in general.
  }
  \item{
    \code{getData} is a bare-bones implementation; it relies on
    a \code{data} argument having been specified and the data being
    available in the environment of the formula. Unlike the functions
    in the \code{nlme} package, it does not do anything special with
    \code{na.action} or \code{subset}.
  }
  }
}
\examples{
  fm1 <- lmer(Reaction~Days+(1|Subject),sleepstudy)
  names(fortify.merMod(fm1))
}

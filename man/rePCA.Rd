\name{rePCA}
\alias{rePCA}
\title{PCA of random-effects covariance matrix}
\usage{
rePCA(x)
}
\arguments{
\item{x}{a \code{merMod} object}
}
\value{
a \code{prcomplist} object
}
\description{
PCA of random-effects variance-covariance estimates
}
\details{
  Perform a Principal Components Analysis (PCA) of the random-effects
  variance-covariance estimates from a fitted mixed-effects model.
  This allows the user to detect and diagnose overfitting problems
  in the random effects model (see Bates et al. 2015 for details).
}
\examples{
  fm1 <- lmer(Reaction~Days+(Days|Subject), sleepstudy)
  rePCA(fm1)
}
\author{
  Douglas Bates
}
\seealso{\code{\link{isSingular}}}
\references{
  \itemize{
    \item Douglas Bates, Reinhold Kliegl, Shravan Vasishth, and Harald Baayen.
    Parsimonious Mixed Models. arXiv:1506.04967 [stat], June 2015. arXiv:
    1506.04967.
  }
}


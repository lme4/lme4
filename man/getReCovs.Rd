\name{getReCovs}
\alias{getReCovs}
\title{Extract Fitted Covariance Structures}
\usage{getReCovs(object)}
\arguments{
  \item{object}{a \code{merMod} object (fitted model)}
}
\value{a list of covariance objects}
\description{Returns a list with one covariance object for each
  random effect in the model. If the object is an old-style (pre-version
  2.0.0) fitted model, it automatically converts all terms to
  \code{Covariance.us} objects (except for terms including double-bar
  \sQuote{||} notation, which are converged to \code{Covariance.diag}
  objects).
}
\examples{
  m1 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
  sapply(getReCovs(m1), "class")
}
  


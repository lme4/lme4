\name{getReCovs}
\alias{getReCovs}
\title{Extract Fitted Covariance Structures}
\usage{getReCovs(object)}
\value{a list of covariance objects}
\description{Returns a list with one covariance object for each
  random effect in the model. If the object is an old-style (pre-version
  2.0.0) fitted model, it automatically converts all terms to
  \code{Covariance.us} objects (except for terms including double-bar
  \sQuote{||} notation, which are converged to \code{Covariance.diag}
  objects).
}
\arguments{
  \item{object}{an \R object of class \code{\linkS4class{merMod}}, i.e.,
    as resulting from \code{\link{lmer}()}, or \code{\link{glmer}()},
    etc.}
}
\examples{
  m1 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
  sapply(getReCovs(m1), "class")
}
  


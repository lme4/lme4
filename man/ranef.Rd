\name{ranef}
\docType{genericFunction}
\docType{methods}
\alias{ranef}
\alias{ranef-methods}
\alias{ranef,mer-method}
\title{Extract the modes of the random effects}
\description{
  A generic function to extract the conditional modes of the random
  effects from a fitted model object.  For linear mixed models the
  conditional modes of the random effects are also the conditional
  means.
}
\usage{
ranef(object, \dots)
\S4method{ranef}{mer}(object, postVar = FALSE, drop = FALSE,
      whichel = names(wt), \dots)
}
\arguments{
  \item{object}{an object of a class of fitted models with random
    effects, typically an \code{"\linkS4class{mer}"} object.}
  \item{postVar}{an optional logical argument indicating if the
    conditional variance-covariance matrices, also called the
    \dQuote{posterior variances}, of the random effects should be
    added as an attribute.  Default is \code{FALSE}.}
  \item{drop}{an optional logical argument indicating components of the
    return value that would be data frames with a single column, usually
    a column called \sQuote{\code{(Intercept)}}, should be returned as
    named vectors.}
  \item{whichel}{character vector of names of factors for which to
    return results.}
  \item{\dots}{some methods for this generic function require additional
    arguments.}
}
\value{
  A list of data frames, one for each grouping factor for the random
  effects.  The number of rows in the data frame is the number of levels
  of the grouping factor.  The number of columns is the dimension of the
  random effect associated with each level of the factor.

  If \code{postVar} is \code{TRUE} each of the data frames has an
  attribute called \code{"postVar"} which is a three-dimensional array
  with symmetric faces.

  When \code{drop} is \code{TRUE} any components that would be data
  frames of a single column are converted to named numeric vectors.
}
\details{
  If grouping factor i has k levels and j random effects per level the ith
  component of the list returned by \code{ranef} is a data frame with k
  rows and j columns.  If \code{postVar} is \code{TRUE} the
  \code{"postVar"} attribute is an array of dimension j by j by k.  The
  kth face of this array is a positive definite symmetric j by j
  matrix.  If there is only one grouping factor in the model the
  variance-covariance matrix for the entire random effects vector,
  conditional on the estimates of the model parameters and on the data
  will be block diagonal and this j by j matrix is the kth diagonal block.
  With multiple grouping factors the faces of the \code{"postVar"}
  attributes are still the diagonal blocks of this conditional
  variance-covariance matrix but the matrix itself is no longer block
  diagonal.
}
\note{
  To produce a \dQuote{caterpillar plot} of the random effects apply
  \code{\link[lattice:xyplot]{dotplot}} to the result of a call to \code{ranef}
  with \code{postVar = TRUE}.
}
\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
fm3 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
ranef(fm1)
str(rr1 <- ranef(fm1, postVar = TRUE))
dotplot(rr1,scales = list(x = list(relation = 'free')))[["Subject"]]
##str(ranef(fm2, postVar = TRUE)) ## code not yet written
op <- options(digits = 4)
ranef(fm3, drop = TRUE)
options(op)
}
\keyword{models}
\keyword{methods}

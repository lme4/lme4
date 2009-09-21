% $Id: lmList.Rd,v 1.4 2003/11/01 04:13:18 deepayan Exp $
\name{lmList}
\alias{lmList}
\alias{lmList,formula,data.frame-method}
\title{List of lm Objects with a Common Model}
\usage{
lmList(formula, data, family, subset, weights,
       na.action, offset, pool, \dots)
}
\arguments{
  \item{formula}{a linear formula object of the form
    \code{y ~ x1+...+xn | g}. In the formula object, \code{y} represents
    the response, \code{x1,...,xn} the covariates, and \code{g} the
    grouping factor specifying the partitioning of the data according to
    which different \code{lm} fits should be performed.}
  \item{data}{a data frame in which to interpret the variables named in
    \code{object}.}
  \item{family}{an optional family specification for a generalized
    linear model.}
  \item{weights}{an optional vector of weights to be used
    in the fitting process.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \dQuote{factory-fresh}
    default is \code{\link{na.omit}}.}
  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.}
 \item{pool}{an optional logical value that is preserved as an attribute
   of the returned value.  This will be used as the default for
   \code{pool} in calculations of standard deviations or standard errors
   for summaries.}
 \item{\dots}{optional arguments to be passed to the model-fitting
   function.}
}
\description{
  The \code{data} argument is split according to the levels of the grouping
  factor \code{g} and individual \code{lm} or \code{glm} fits are
  obtained for each \code{data} partition, using the model defined in
  \code{object}.
}
\value{
  an object of class \code{"\linkS4class{lmList}"}, which is a list of
  \code{lm} objects with as many components as the number of groups
  defined by the grouping factor.
}
\seealso{
\code{\link[stats]{lm}}
}
\examples{
(fm1 <- lmList(Reaction ~ Days | Subject, sleepstudy))
}
\keyword{models}


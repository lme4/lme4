\name{ranef}
\alias{ranef}
\alias{ranef.merMod}
\alias{dotplot.ranef.mer}
\alias{qqmath.ranef.mer}
\alias{as.data.frame.ranef.mer}
\title{Extract the modes of the random effects}
\usage{
\method{ranef}{merMod} (object, condVar = TRUE,
      drop = FALSE, whichel = names(ans), postVar = FALSE, ...)
\method{dotplot}{ranef.mer} (x, data, main = TRUE, transf = I, ...)
\method{qqmath}{ranef.mer} (x, data, main = TRUE, ...)
\method{as.data.frame}{ranef.mer} (x, ..., stringsAsFactors = default.stringsAsFactors())
}
\arguments{
  \item{object}{an object of a class of fitted models with
    random effects, typically a
    \code{\linkS4class{merMod}} object.}
  \item{condVar}{a logical argument indicating if
    the conditional variance-covariance matrices of the random effects
    should be added as an attribute.}
  \item{drop}{should components of the return value that would be data frames
    with a single column, usually a column called
    \sQuote{\code{(Intercept)}}, be returned as named
    vectors instead?}
  \item{whichel}{character vector of names of
    grouping factors for which the random effects should be
    returned.}
  \item{postVar}{a (deprecated) synonym for \code{condVar}}
  \item{x}{a random-effects object (of class \code{ranef.mer})
    produced by \code{ranef}}
  \item{main}{include a main title, indicating the grouping factor,
    on each sub-plot?}
  \item{transf}{transformation for random effects: for example,
    \code{exp} for plotting parameters from a (generalized)
    logistic regression on the odds rather than log-odds scale}
  \item{data}{This argument is required by the \code{dotplot}
    and \code{qqmath} generic methods, but is not actually used.}
  \item{stringsAsFactors}{see \code{\link{data.frame}}}
  \item{\dots}{some methods for these generic functions
    require additional arguments.}
}
\value{
  \itemize{
    \item{
      From \code{ranef}:
      An object of class \code{ranef.mer} composed of
      a list of data frames, one for each grouping factor for
      the random effects.  The number of rows in the data frame
      is the number of levels of the grouping factor.  The
      number of columns is the dimension of the random effect
      associated with each level of the factor.

      If \code{condVar} is \code{TRUE} each of the data frames
      has an attribute called \code{"postVar"}.
      \itemize{
	\item If there is
	a single random-effects term for a given grouping factor,
	this attribute is a
	three-dimensional array with symmetric faces; each face
	contains the variance-covariance matrix for a particular
	level of the grouping factor.
	\item If there is more than one random-effects term for
	a given grouping factor (e.g. \code{(1|f) + (0+x|f)}),
	this attribute is a list of arrays as described above,
	one for each term.
      }
      (The name
      of this attribute is a historical artifact,
      and may be changed to \code{condVar} at some point in the future.)

      When \code{drop} is \code{TRUE} any components that would
      be data frames of a single column are converted to named
      numeric vectors.
    }
    \item{
      From \code{as.data.frame}:

      This function converts the random effects to a "long format"
      data frame with columns
      \describe{
	\item{grpvar}{grouping variable}
	\item{term}{random-effects term, e.g. \dQuote{(Intercept)} or
	  \dQuote{Days}}
	\item{grp}{level of the grouping variable (e.g., which Subject)}
	\item{condval}{value of the conditional mean}
	\item{condsd}{conditional standard deviation}
      }
    } % as.data.frame value
  } % itemize
} % value
\description{
  A generic function to extract the conditional modes of
  the random effects from a fitted model object.  For
  linear mixed models the conditional modes of the random
  effects are also the conditional means.
}
\details{
  If grouping factor i has k levels and j random effects
  per level the ith component of the list returned by
  \code{ranef} is a data frame with k rows and j columns.
  If \code{condVar} is \code{TRUE} the \code{"postVar"}
  attribute is an array of dimension j by j by k (or a list
  of such arrays).  The kth
  face of this array is a positive definite symmetric j by
  j matrix.  If there is only one grouping factor in the
  model the variance-covariance matrix for the entire
  random effects vector, conditional on the estimates of
  the model parameters and on the data, will be block
  diagonal; this j by j matrix is the kth diagonal
  block.  With multiple grouping factors the faces of the
  \code{"postVar"} attributes are still the diagonal blocks
  of this conditional variance-covariance matrix but the
  matrix itself is no longer block diagonal.
}
\note{
  To produce a (list of) \dQuote{caterpillar plots} of the random
  effects apply \code{\link[lattice:xyplot]{dotplot}} to
  the result of a call to \code{ranef} with \code{condVar =
    TRUE}; \code{\link[lattice:xyplot]{qqmath}} will generate
  a list of Q-Q plots.
}
\examples{
require(lattice)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
fm3 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
ranef(fm1)
str(rr1 <- ranef(fm1))
dotplot(rr1)  ## default
## specify free scales in order to make Day effects more visible
dotplot(rr1,scales = list(x = list(relation = 'free')))[["Subject"]]
ranef(fm2)
op <- options(digits = 4)
ranef(fm3, drop = TRUE)
options(op)
## as.data.frame() provides RE's and conditional standard deviations:
str(dd <- as.data.frame(rr1))
if (require(ggplot2)) {
    ggplot(dd, aes(y=grp,x=condval)) +
        geom_point() + facet_wrap(~term,scales="free_x") +
        geom_errorbarh(aes(xmin=condval -2*condsd,
                           xmax=condval +2*condsd), height=0)
}
}
\keyword{methods}
\keyword{models}


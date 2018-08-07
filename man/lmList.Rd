\name{lmList}
\alias{lmList}
\alias{plot.lmList}
\title{Fit List of lm or glm Objects with a Common Model}
\description{
  Fit a list of \code{\link{lm}} or \code{\link{glm}} objects with a
  common model for different subgroups of the data.
}
\usage{
lmList(formula, data, family, subset, weights, na.action,
       offset, pool = !isGLM || .hasScale(family2char(family)),
       warn = TRUE, \dots)
}
\arguments{
  \item{formula}{a linear \code{\link{formula}} object of the form
    \code{y ~ x1+...+xn | g}.  In the formula object, \code{y}
    represents the response, \code{x1,...,xn} the covariates,
    and \code{g} the grouping factor specifying the
    partitioning of the data according to which different
    \code{lm} fits should be performed.}
  \item{family}{an optional \code{\link{family}} specification for a
    generalized linear model (\code{\link{glm}}).}
  \item{data}{an optional data frame containing the
    variables named in \code{formula}.  By default the
    variables are taken from the environment from which
    \code{lmer} is called.  See Details.}
  \item{subset}{an optional expression indicating the
    subset of the rows of \code{data} that should be used in
    the fit.  This can be a logical vector, or a numeric
    vector indicating which observation numbers are to be
    included, or a character vector of the row names to be
    included.  All observations are included by default.}
  \item{weights}{an optional vector of \sQuote{prior
      weights} to be used in the fitting process.  Should be
    \code{NULL} or a numeric vector.}
  \item{na.action}{a function that indicates what should
    happen when the data contain \code{NA}s.  The default
    action (\code{na.omit}, inherited from the \sQuote{factory
    fresh} value of \code{getOption("na.action")}) strips any
    observations with any missing values in any variables.}
  \item{offset}{this can be used to specify an \emph{a
      priori} known component to be included in the linear
    predictor during fitting.  This should be \code{NULL} or a
    numeric vector of length equal to the number of cases.
    One or more \code{\link{offset}} terms can be included in
    the formula instead or as well, and if more than one is
    specified their sum is used.  See
    \code{\link{model.offset}}.}
  \item{pool}{logical scalar indicating if the variance estimate should
    pool the residual sums of squares.  By default true if the model has
    a scale parameter (which includes all linear, \code{\link{lmer}()}, ones).}
  \item{warn}{indicating if errors in the single fits should signal a
    \dQuote{summary} \code{\link{warning}}.}
  \item{\dots}{additional, optional arguments to be passed to
    the model function or family evaluation.}
}
\value{an object of \code{\link{class}} \code{\linkS4class{lmList4}} (see
  there, notably for the \code{\link{methods}} defined).
}
\details{
  \itemize{
    \item While \code{data} is optional, the package authors
    \emph{strongly} recommend its use, especially when later applying
    methods such as \code{update} and \code{drop1} to the fitted model
    (\emph{such methods are not guaranteed to work properly if
      \code{data} is omitted}).  If \code{data} is omitted, variables will
    be taken from the environment of \code{formula} (if specified as a
    formula) or from the parent frame (if specified as a character vector).

    \item Since \pkg{lme4} version 1.1-16, if there are errors (see
    \code{\link{stop}}) in the single (\code{lm()} or \code{glm()})
    fits, they are summarized to a warning message which is returned as
    attribute \code{"warnMessage"} and signalled as \code{\link{warning}()}
    when the \code{warn} argument is true.

    In previous \pkg{lme4} versions, a general (different) warning had
    been signalled in this case.
  }
}
\examples{
fm.plm  <- lmList(Reaction ~ Days | Subject, sleepstudy)
coef(fm.plm)
fm.2  <- update(fm.plm, pool = FALSE)
## coefficients are the same, "pooled or unpooled":
stopifnot( all.equal(coef(fm.2), coef(fm.plm)) )

(ci <- confint(fm.plm)) # print and rather *see* :
plot(ci)                # how widely they vary for the individuals
}
\seealso{\code{\linkS4class{lmList4}}}
\keyword{models}
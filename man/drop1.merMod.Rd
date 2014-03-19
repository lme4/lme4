\name{drop1.merMod}
\alias{drop1.merMod}
\title{Drop all possible single fixed-effect terms from a mixed effect model}
\description{
  Drop allowable single terms from the model: see \code{\link{drop1}}
  for details of how the appropriate scope for dropping terms
  is determined.
}
\usage{
\method{drop1}{merMod}(object, scope, scale = 0,
    test = c("none", "Chisq", "user"),
    k = 2, trace = FALSE, sumFun, ...)
}
\arguments{
  \item{object}{a fitted \code{merMod} object.}
  \item{scope}{a formula giving the terms to be considered for adding or
    dropping.}
  \item{scale}{Currently ignored (included for S3 method compatibility)}
  \item{test}{should the results include a test statistic relative to the
    original model?
    The \eqn{\chi^2}{Chisq} test is a likelihood-ratio test,
    which is approximate due to finite-size effects.
  }
  \item{k}{the penalty constant in AIC}
  \item{trace}{print tracing information?}
  \item{sumFun}{a summary function to be used when
    \code{test=="user"}. It must allow arguments \code{scale} and
    \code{k}, but these may be ignored (e.g. specified in \code{dots}).
    The first two arguments must be \code{object}, the full model fit,
    and \code{objectDrop}, a reduced model.  If \code{objectDrop} is missing,
    \code{sumFun} should return a vector of with the appropriate
    length and names (the actual contents are ignored).
  }
  \item{\dots}{other arguments (ignored)}
}
\value{
  An object of class \code{anova} summarizing the differences in fit
  between the models.
}
\details{
  \code{drop1} relies on being able to find the appropriate information
  within the environment of the formula of the original model.  If the
  formula is created in an environment that does not contain the data,
  or other variables passed to the original model (for example, if
  a separate function is called to define the formula), then
  \code{drop1} will fail.  A workaround (see example below) is to
  manually specify an appropriate environment for the formula.
}
\examples{
fm1 <- lmer(Reaction~Days+(Days|Subject),sleepstudy)
## likelihood ratio tests
drop1(fm1,test="Chisq")
## use Kenward-Roger corrected F test, or parametric bootstrap,
## to test the significance of each dropped predictor
if (require(pbkrtest) && packageVersion("pbkrtest")>="0.3.8") {
   KRSumFun <- function(object, objectDrop, ...) {
      krnames <- c("ndf","ddf","Fstat","p.value","F.scaling")
      r <- if (missing(objectDrop)) {
          setNames(rep(NA,length(krnames)),krnames)
      } else {
         krtest <- KRmodcomp(object,objectDrop)
         unlist(krtest$stats[krnames])
      }
      attr(r,"method") <- c("Kenward-Roger via pbkrtest package")
      r
   }
   drop1(fm1,test="user",sumFun=KRSumFun)

   if(lme4:::testLevel() >= 3) { ## takes about 16 sec
     nsim <- 100
     PBSumFun <- function(object, objectDrop, ...) {
	pbnames <- c("stat","p.value")
	r <- if (missing(objectDrop)) {
	    setNames(rep(NA,length(pbnames)),pbnames)
	} else {
	   pbtest <- PBmodcomp(object,objectDrop,nsim=nsim)
	   unlist(pbtest$test[2,pbnames])
	}
	attr(r,"method") <- c("Parametric bootstrap via pbkrtest package")
	r
     }
     system.time(drop1(fm1,test="user",sumFun=PBSumFun))
   }
}
## workaround for creating a formula in a separate environment
createFormula <- function(resp, fixed, rand) {  
    f <- reformulate(c(fixed,rand),response=resp)
    ## use the parent (createModel) environment, not the
    ## environment of this function (which does not contain 'data')
    environment(f) <- parent.frame()
    f
}
createModel <- function(data) {
    mf.final <- createFormula("Reaction", "Days", "(Days|Subject)")
    lmer(mf.final, data=data)
}
drop1(createModel(data=sleepstudy))
}
\keyword{misc}

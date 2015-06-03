\name{refit}
\alias{refit}
\alias{refit.merMod}
\title{Refit a (merMod) Model with a Different Response}
\description{
   Refit a model, possibly after modifying the response vector.  This
   makes use of the model representation and directly goes to the
   optimization.
}
\usage{
refit(object, newresp, \dots)

\method{refit}{merMod}(object, newresp = NULL, rename.response = FALSE,
      maxit = 100, \dots)
}
\arguments{
  \item{object}{a fitted model, usually of class
    \code{\linkS4class{lmerMod}}, to be refit with a new response.}
  \item{newresp}{an (optional) numeric vector providing the new
    response, of the same length as the original response (see
    \code{Details} for information on NA handling).  May also be a data
    frame with a single numeric column, e.g. as produced by
    \code{simulate(object)}.}
  \item{rename.response}{when refitting the model, should the name of
    the response variable in the formula and model frame be replaced with
    the name of \code{newresp}?}
  \item{maxit}{scalar integer, currently only for GLMMs: the maximal
    number of Pwrss update iterations.}
  \item{\dots}{optional additional parameters.  For the \code{merMod}
    method, \code{control}.}
}
\value{
  an object like \code{x}, but fit to a different response vector \eqn{Y}.
}
\seealso{
  \code{\link{update.merMod}} for more flexible and extensive model
  refitting;
  \code{\link{refitML}} for refitting a REML fitted model with maximum
  likelihood (\sQuote{ML}).
}
\details{
  Refit a model, possibly after modifying the response vector.  This
  could be done using \code{\link{update}()}, but the \code{refit()}
  approach should be faster because it bypasses the creation of the
  model representation and goes directly to the optimization step.

  Setting \code{rename.response = TRUE} may be necessary if one
  wants to do further operations (such as \code{update}) on the fitted
  model.  However, the refitted model will still be slightly different
  from the equivalent model fitted via \code{update}; in particular, the
  \code{terms} component is not updated to reflect the new response
  variable, if it has a different name from the original.

  If \code{newresp} has an \code{na.action} attribute, then it is
  assumed that NA values have already been removed from the numeric
  vector; this allows the results of \code{\link{simulate}(object)} to
  be used even if the original response vector contained NA values.
  Otherwise, the length of \code{newresp} must be the same as the
  \emph{original} length of the response.
}
\examples{
## Ex. 1: using refit() to fit each column in a matrix of responses -------
set.seed(101)
Y <- matrix(rnorm(1000),ncol=10)
## combine first column of responses with predictor variables
d <- data.frame(y=Y[,1],x=rnorm(100),f=rep(1:10,10))
## (use check.conv.grad="ignore" to disable convergence checks because we
##  are using a fake example)
## fit first response
fit1 <- lmer(y ~ x+(1|f), data = d,
             control= lmerControl(check.conv.grad="ignore",
                                  check.conv.hess="ignore"))
## combine fit to first response with fits to remaining responses
res <- c(fit1,lapply(as.data.frame(Y[,-1]), refit, object=fit1))

## Ex. 2: refitting simulated data using data that contain NA values ------
sleepstudyNA <- sleepstudy
sleepstudyNA$Reaction[1:3] <- NA
fm0 <- lmer(Reaction ~ Days + (1|Subject), sleepstudyNA)
## the special case of refitting with a single simulation works ...
ss0 <- refit(fm0, simulate(fm0))
## ... but if simulating multiple responses (for efficiency),
## need to use na.action=na.exclude in order to have proper length of data
fm1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudyNA, na.action=na.exclude)
ss <- simulate(fm1, 5)
res2 <- refit(fm1, ss[[5]])
}
%% unit tests --> ../tests/refit.R
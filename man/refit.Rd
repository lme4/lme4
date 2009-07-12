\name{refit}
\alias{refit}
\alias{refit,mer,numeric-method}
\alias{refit,mer,matrix-method}
\title{Re-fit a model to a new response vector}
\description{

  This generic function fits a model to a new response vector.  It
  typically is used as part of a simulation.  For models fit by
  \code{\link{lmer}} or \code{\link{nlmer}} it is much faster to use
  this function to re-fit the model than to begin from the model
  formula.

}
\usage{
refit(object, newresp, \dots)
}
\arguments{
  \item{object}{a fitted model.  Methods are defined for models fit by
    \code{\link{lmer}} and by \code{\link{nlmer}}}.
  \item{newresp}{a new response vector - typically a numeric vector.}
  \item{\dots}{further arguments passed for some methods.}
}
\details{
  For many model-fitting functions in R a large part of the execution
  time is taken up with converting a formula representation and the data
  argument into the numerical structures used to evaluate the parameter
  estimates.  It is common in a simulation to use the same model
  specification and covariate data on many simulated responses.  That
  is, the only thing that changes between model fits is the response
  vector. This generic function allows for the model specification to
  and covariate data to stay the same while the only the response vector
  is changed.
}
\value{
  a fitted model similar to \code{object}
}
%\references{}
%\author{}
%\note{}
%\seealso{}
\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
set.seed(54321)
system.time(simfe <-
        data.frame(t(apply(simulate(fm1, nsim = 100),
                           2, function(y) fixef(refit(fm1, y)))),
                   check.names = FALSE))
qqmath(~ `(Intercept)` + Days, simfe, outer = TRUE,
       ylab = NULL, scales = list(y = list(relation = "free")),
       layout = c(1,2), aspect = 1, type = c("g", "p"))
sapply(simfe, mean)
sapply(simfe, sd)
}
\keyword{models}
\keyword{methods}

\name{predict.merMod}
\alias{predict.merMod}
\title{Predictions from a model at new data values}
\usage{
  \method{predict}{merMod} (object, newdata = NULL,
    newparams=NULL, newX=NULL,
    re.form = NULL, ReForm, REForm, REform,
    terms = NULL,
    type = c("link", "response"), allow.new.levels = FALSE,
    na.action = na.pass, \dots)
}
\arguments{
  \item{object}{a fitted model object}
  \item{newdata}{data frame for which to evaluate
    predictions}
  \item{newparams}{new parameters to use in evaluating predictions,
    specified as in the \code{start} parameter for \code{\link{lmer}}
    or \code{\link{glmer}} -- a list with components \code{theta} and/or
    (for GLMMs) \code{beta}}
  \item{newX}{new design matrix to use in evaluating predictions
    (alternative to \code{newdata})}
  \item{re.form}{formula for random effects to condition on.  If \code{NULL},
    include all random effects; if \code{NA} or \code{~0},
    include no random effects}
  \item{ReForm}{allowed for backward compatibility: \code{re.form} is
    now the preferred argument name}
  \item{REForm}{allowed for backward compatibility: \code{re.form} is
    now the preferred argument name}
  \item{REform}{allowed for backward compatibility: \code{re.form} is
    now the preferred argument name}
  \item{terms}{a \code{\link{terms}} object - not used at
    present}
  \item{type}{character string - either \code{"link"}, the
    default, or \code{"response"} indicating the type of
    prediction object returned}
  \item{allow.new.levels}{(logical) if FALSE (default),
    then any new levels (or NA values) detected in
    \code{newdata} will trigger an error; if TRUE, then the
    prediction will use the unconditional (population-level)
    values for data with previously unobserved levels (or
    NAs)}
  
  \item{na.action}{function determining what should be done
    with missing values for fixed effects in \code{newdata}.
    The default is to predict \code{NA}: see
    \code{\link{na.pass}}.}
  \item{...}{optional additional parameters.  None are used
    at present.}
}
\value{
  a numeric vector of predicted values
}
\description{
  \code{\link{predict}} method for
  \code{\linkS4class{merMod}} objects
}
\note{
  There is no option for computing standard errors of
  predictions because it is difficult to define an
  efficient method that incorporates uncertainty in the
  variance parameters; we recommend \code{\link{bootMer}}
  for this task.
}
\examples{
(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 |herd), cbpp, binomial))
str(p0 <- predict(gm1))            # fitted values
str(p1 <- predict(gm1,ReForm=NA))  # fitted values, unconditional (level-0)
newdata <- with(cbpp, expand.grid(period=unique(period), herd=unique(herd)))
str(p2 <- predict(gm1,newdata))    # new data, all RE
str(p3 <- predict(gm1,newdata,ReForm=NA)) # new data, level-0
str(p4 <- predict(gm1,newdata,ReForm=~(1|herd))) # explicitly specify RE
}


\name{predict.merMod}
\alias{predict.merMod}
\title{Predictions from a model at new data values}
\usage{
  \method{predict}{merMod} (object, newdata = NULL,
    REform = NULL, terms = NULL,
    type = c("link", "response"), allow.new.levels = FALSE,
    na.action = na.pass, ...)
}
\arguments{
  \item{object}{a fitted model object}

  \item{newdata}{data frame for which to evaluate
  predictions}

  \item{REform}{formula for random effects to include.  If
  NULL, include all random effects; if NA, include no
  random effects}

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
str(p1 <- predict(gm1,REform=NA))  # fitted values, unconditional (level-0)
newdata <- with(cbpp, expand.grid(period=unique(period), herd=unique(herd)))
str(p2 <- predict(gm1,newdata))    # new data, all RE
str(p3 <- predict(gm1,newdata,REform=NA)) # new data, level-0
str(p4 <- predict(gm1,newdata,REform=~(1|herd))) # explicitly specify RE
}


\name{predict.merMod}
\alias{predict.merMod}
\title{Predictions from a model at new data values}
\description{
  The \code{\link{predict}} method for \code{\linkS4class{merMod}}
  objects, i.e. results of \code{\link{lmer}()}, \code{\link{glmer}()}, etc.
}
\usage{
\method{predict}{merMod}(object, newdata = NULL, newparams = NULL,
	re.form = NULL, ReForm, REForm, REform,
        random.only=FALSE, terms = NULL,
	type = c("link", "response"), allow.new.levels = FALSE,
	na.action = na.pass, \dots)
}
\arguments{
  \item{object}{a fitted model object}
  \item{newdata}{data frame for which to evaluate
    predictions.}
  \item{newparams}{new parameters to use in evaluating predictions,
    specified as in the \code{start} parameter for \code{\link{lmer}}
    or \code{\link{glmer}} -- a list with components \code{theta} and/or
    (for GLMMs) \code{beta}.}
  \item{re.form}{formula for random effects to condition on.  If \code{NULL},
    include all random effects; if \code{NA} or \code{~0},
    include no random effects.}
  \item{ReForm, REForm, REform}{allowed for backward compatibility: \code{re.form} is
    now the preferred argument name.}
  \item{random.only}{(logical) ignore fixed effects, making predictions
    only using random effects?}
  \item{terms}{a \code{\link{terms}} object - unused at present.}
  \item{type}{character string - either \code{"link"}, the default, or
    \code{"response"} indicating the type of prediction object returned.}
  \item{allow.new.levels}{logical if new levels (or NA values) in
    \code{newdata} are allowed. If FALSE (default), such new values in
    \code{newdata} will trigger an error; if TRUE, then the prediction
    will use the unconditional (population-level) values for data with
    previously unobserved levels (or NAs).}
  \item{na.action}{\code{\link{function}} determining what should be done
    with missing values for fixed effects in \code{newdata}.
    The default is to predict \code{NA}: see \code{\link{na.pass}}.}
  \item{...}{optional additional parameters.  None are used at present.}
}
\value{
  a numeric vector of predicted values
}
\details{
  \itemize{
    \item If any random effects are included in \code{re.form}
    (see below), \code{newdata} \emph{must} contain columns
    corresponding to all of the grouping variables and
    random effects used in the original model, even if not all
    are used in prediction; however, they can be safely set to \code{NA}
    in this case.
    \item There is no option for computing standard errors of
    predictions because it is difficult to define an
    efficient method that incorporates uncertainty in the
    variance parameters; we recommend \code{\link{bootMer}}
    for this task.
  }
}
\examples{
(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 |herd), cbpp, binomial))
str(p0 <- predict(gm1))            # fitted values
str(p1 <- predict(gm1,re.form=NA))  # fitted values, unconditional (level-0)
newdata <- with(cbpp, expand.grid(period=unique(period), herd=unique(herd)))
str(p2 <- predict(gm1,newdata))    # new data, all RE
str(p3 <- predict(gm1,newdata,re.form=NA)) # new data, level-0
str(p4 <- predict(gm1,newdata,re.form= ~(1|herd))) # explicitly specify RE
stopifnot(identical(p2, p4))
\dontshow{%back compatibility - as long as we have it: ReForm, REForm, REform,
suppressWarnings(stopifnot(identical(p1, predict(gm1,ReForm=NA))))
suppressWarnings(stopifnot(identical(p3, predict(gm1,newdata,REForm=NA))))
suppressWarnings(stopifnot(identical(p4, predict(gm1,newdata,REform=~(1|herd)))))

## predict() should work with variable names with spaces [as lm() does]:
dd <- expand.grid(y=1:3, "Animal ID" = 1:9)
fm <- lmer(y ~ 1 + (1 | `Animal ID`),  dd)
summary(fm)
isel <- c(7, 9, 11, 13:17, 20:22)
stopifnot(all.equal(vcov(fm)[1,1], 0.02564102564),
	  all.equal(unname(predict(fm, newdata = dd[isel,])),
		    unname( fitted(fm) [isel])))
}%dont
}


\name{mkReTrms}
\alias{mkReTrms}
\title{Make Random Effect Terms: Create Z, Lambda, Lind, etc.}
\description{
  From the result of \code{\link{findbars}} applied to a model formula
  and the evaluation frame \code{fr}, create the model matrix
  \code{Zt}, etc, associated with the random-effects terms.
}
\usage{
  mkReTrms(bars, fr, drop.unused.levels=TRUE)
}
\arguments{
  \item{bars}{a list of parsed random-effects terms}
  \item{fr}{a model frame in which to evaluate these terms}
  \item{drop.unused.levels}{(logical) drop unused factor levels? (experimental)}
}
\value{
  a \code{\link{list}} with components
  \item{Zt}{transpose of the sparse model matrix for the random effects}
  \item{Lambdat}{transpose of the sparse relative covariance factor}
  \item{Lind}{an integer vector of indices determining the mapping of
    the elements of the \code{theta} vector to the \code{"x"} slot of
    \code{Lambdat}}
  \item{theta}{initial values of the covariance parameters}
  \item{lower}{lower bounds on the covariance parameters}
  \item{flist}{list of grouping factors used in the
    random-effects terms}
  \item{cnms}{a list of column names
    of the random effects according to the grouping factors}
}
\seealso{
  Other utilities: \code{\link{findbars}},
  \code{\link{mkRespMod}}, \code{\link{nlformula}},
  \code{\link{nobars}}, \code{\link{subbars}}
}
\examples{
data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
(bar.f <- findbars(mform)) # list with 3 terms
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FIXME
## Unfinished _FIXME_
## Pfrm <- model.frame(bar.f, data = Pixel)
}
\keyword{utilities}
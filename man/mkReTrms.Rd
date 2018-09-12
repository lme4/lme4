\name{mkReTrms}
\alias{mkReTrms}
\title{Make Random Effect Terms: Create Z, Lambda, Lind, etc.}
\description{
  From the result of \code{\link{findbars}} applied to a model formula
  and the evaluation frame \code{fr}, create the model matrix
  \code{Zt}, etc, associated with the random-effects terms.
}
\usage{
  mkReTrms(bars, fr, drop.unused.levels=TRUE,
                     reorder.terms=TRUE,
                     reorder.vars=FALSE) 
}
\arguments{
  \item{bars}{a list of parsed random-effects terms}
  \item{fr}{a model frame in which to evaluate these terms}
  \item{drop.unused.levels}{(logical) drop unused factor levels?
    (experimental)}
  \item{reorder.terms}{arrange random effects terms in decreasing order
    of number of groups (factor levels)?}
  \item{reorder.vars}{arrange columns of individual random effects terms
    in alphabetical order?}
}
\value{
  a \code{\link{list}} with components
  \item{Zt}{transpose of the sparse model matrix for the random effects}
  \item{theta}{initial values of the covariance parameters}
  \item{Lind}{an integer vector of indices determining the mapping of
    the elements of the \code{theta} vector to the \code{"x"} slot of
    \code{Lambdat}}
  \item{Gp}{a vector indexing the association of
    elements of the conditional mode vector
    with random-effect terms; if \code{nb} is the vector of numbers
    of conditional modes per term (i.e. number of groups times number
    of effects per group), \code{Gp} is \code{c(0,cumsum(nb))}
    (and conversely \code{nb} is \code{diff(Gp)})}
  \item{lower}{lower bounds on the covariance parameters}
  \item{Lambdat}{transpose of the sparse relative covariance factor}
  \item{flist}{list of grouping factors used in the
    random-effects terms}
  \item{cnms}{a list of column names
    of the random effects according to the grouping factors}
  \item{Ztlist}{list of components of the transpose of the
    random-effects model matrix, separated by random-effects term}
}
\seealso{
  Other utilities: \code{\link{findbars}},
  \code{\link{mkRespMod}}, \code{\link{nlformula}},
  \code{\link{nobars}}, \code{\link{subbars}}.
  \code{\link{getME}} can retrieve these components
  from a fitted model, although their values and/or forms
  may be slightly different in the final fitted model from
  their original values as returned from \code{mkReTrms}.
}
\examples{
data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
(bar.f <- findbars(mform)) # list with 3 terms
mf <- model.frame(subbars(mform),data=Pixel)
rt <- mkReTrms(bar.f,mf)
names(rt)
}
%fm1 <- lmer(mform,Pixel)
%rt2 <- getME(fm1,names(rt))
%for (i in seq_along(rt))
%    cat(names(rt)[[i]],
%        isTRUE(all.equal(rt[[i]],rt2[[i]])),
%        "\n")
% ## theta and Lambda components have the same structure, but have
% ## been updated in the fitting; flist is a data frame in reTrms
% ## but converted to a list; Ztlist is quite different (has been
% ## decomposed in getME()
\keyword{utilities}
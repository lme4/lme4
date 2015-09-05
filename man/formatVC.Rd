\name{formatVC}
\title{Format the 'VarCorr' Matrix of Random Effects}
\alias{formatVC}
\description{
  \dQuote{format()} the \code{\link{VarCorr}} matrix of the
  random effects -- for \code{\link{print}()}ing and \code{\link{show}()}ing
}
\usage{
formatVC(varc, digits = max(3, getOption("digits") - 2),
         comp = "Std.Dev.", formatter = format, \dots)
}
\arguments{
  \item{varc}{a \code{\link{VarCorr}}(-like) matrix with
    \code{\link{attributes}}.}
  \item{digits}{the number of significant digits.}
  \item{comp}{character vector of length one
    or two indicating which columns out of "Variance" and
    "Std.Dev." should be shown in the formatted output.}
  \item{formatter}{the \code{\link{function}} to be used
    for formatting the standard deviations and or variances
    (but \emph{not} the correlations which (currently) are
    always formatted as "0.nnn").}
  \item{\dots}{optional arguments for \code{formatter(*)} in addition to
    the first (numeric vector) and \code{digits}.}
}
\value{
  a character matrix of formatted VarCorr entries from
  \code{varc}.
}
%% \details{      ___TODO___
%% }
%% \references{
%% }
\seealso{
  \code{\link{VarCorr}()}.
}
\examples{
data(Orthodont, package="nlme")
fm1 <- lmer(distance ~ age + (age|Subject), data = Orthodont)

vc <- VarCorr(fm1)
## Now print()ing effectively works via  formatVC() :
vc  ## default print method: standard dev and corr
m <- formatVC(vc)
stopifnot(is.matrix(m), is.character(m), ncol(m) == 4)
## this is what you also saw above in 'vc', i.e.,  print(vc):
print(m, quote = FALSE)
}
\keyword{methods}

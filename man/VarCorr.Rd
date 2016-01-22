\name{VarCorr}
\title{Extract Variance and Correlation Components}
\alias{VarCorr}
\alias{VarCorr.merMod}
\alias{as.data.frame.VarCorr.merMod}
\alias{print.VarCorr.merMod}
\usage{
\method{VarCorr}{merMod}(x, sigma=1, \dots)

%% documented because of "surprising arguments":
\method{as.data.frame}{VarCorr.merMod}(x, row.names = NULL,
              optional = FALSE, order = c("cov.last", "lower.tri"), \dots)
\method{print}{VarCorr.merMod}(x, digits = max(3, getOption("digits") - 2),
      comp = "Std.Dev.", formatter = format, \dots)
}
\arguments{
  \item{x}{for \code{VarCorr}: a fitted model object, usually an object inheriting from
    class \code{\linkS4class{merMod}}. For \code{as.data.frame}, a
    \code{VarCorr.merMod} object returned from \code{VarCorr}.}
  \item{sigma}{an optional numeric value used as a multiplier for the
    standard deviations.}
  \item{digits}{an optional integer value specifying the number of digits}
  \item{order}{arrange data frame with variances/standard deviations
    first and covariances/correlations last for each random effects
    term (\code{"cov.last"}), or in the order of the lower triangle of the
    variance-covariance matrix (\code{"lower.tri"})?}
  \item{row.names, optional}{Ignored: necessary for the
    \code{as.data.frame} method.}
  \item{\dots}{Ignored for the \code{as.data.frame} method; passed to
    other \code{\link{print}()} methods for the \code{print()} method.}
  %% print():
  \item{comp}{a \code{\link{character}} vector, specifying the components
    to be printed; simply passed to \code{\link{formatVC}()}.}
  \item{formatter}{a \code{\link{function}} for formatting the numbers;
    simply passed to \code{\link{formatVC}()}.}
}
\description{
  This function calculates the estimated variances,
  standard deviations, and correlations between the
  random-effects terms in a mixed-effects model, of class
  \code{\linkS4class{merMod}} (linear, generalized or
  nonlinear).  The within-group error variance and standard
  deviation are also calculated.
}
\value{
  An object of class \code{VarCorr.merMod}. The internal
  structure of the object is
  a list of matrices, one for each random effects grouping
  term.  For each grouping term, the standard deviations and
  correlation matrices for each grouping term are stored as
  attributes \code{"stddev"} and \code{"correlation"},
  respectively, of the variance-covariance matrix, and the
  residual standard deviation is stored as attribute
  \code{"sc"} (for \code{glmer} fits, this attribute stores
  the scale parameter of the model).

The \code{as.data.frame} method produces a combined data frame with one
row for each variance or covariance parameter (and a row for the
residual error term where applicable) and the following columns:
  \describe{
    \item{grp}{grouping factor}
    \item{var1}{first variable}
    \item{var2}{second variable (\code{NA} for variance parameters)}
    \item{vcov}{variances or covariances}
    \item{sdcor}{standard deviations or correlations}
  }
}
\details{
  The \code{print} method for \code{VarCorr.merMod} objects
  has optional arguments \code{digits} (specify digits of
  precision for printing) and \code{comp}: the latter is
  a character vector with any combination of \code{"Variance"}
  and \code{"Std.Dev."}, to specify whether variances,
  standard deviations, or both should be printed.
}
\author{
  This is modeled after \code{\link[nlme]{VarCorr}} from
  package \pkg{nlme}, by Jose Pinheiro and Douglas Bates.
}
\seealso{
  \code{\link{lmer}}, \code{\link{nlmer}}
}
\examples{
data(Orthodont, package="nlme")
fm1 <- lmer(distance ~ age + (age|Subject), data = Orthodont)
(vc <- VarCorr(fm1))  ## default print method: standard dev and corr
## both variance and std.dev.
print(vc,comp=c("Variance","Std.Dev."),digits=2)
## variance only
print(vc,comp=c("Variance"))
as.data.frame(vc)
as.data.frame(vc,order="lower.tri")
}
\keyword{models}


\name{simulate.formula}
\alias{simulate.formula}
\title{A \code{simulate} Method for \code{formula} objects that dispatches based on the Left-Hand Side}
\description{
This method evaluates the left-hand side (LHS) of the given formula and
dispatches it to an appropriate method based on the result by
setting an nonce class name on the formula.
}
\usage{
\method{simulate}{formula}(object, nsim = 1 , seed = NULL, \dots,
basis, newdata, data)
}
\arguments{
\item{object}{a one- or two-sided \code{\link{formula}}.}
\item{nsim,seed}{number of realisations to simulate and the random
  seed to use; see \code{\link{simulate}}}
\item{...}{additional arguments to methods}
\item{basis}{if given, overrides the LHS of the formula for the
  purposes of dispatching}
\item{newdata,data}{if passed, the \code{object}'s LHS is evaluated in
  this environment; at most one of the two may be passed.}
}
\details{
The dispatching works as follows:

\enumerate{
\item If \code{basis} is not passed, and the formula has an LHS the
   expression on the LHS of the formula in the \code{object} is
   evaluated in the environment \code{newdata} or \code{data} (if given), in
   any case enclosed by the environment of \code{object}. Otherwise,
   \code{basis} is used.
\item The result is set as an attribute \code{".Basis"} on \code{object}. If
   there is no \code{basis} or LHS, it is not set.
\item The class vector of \code{object} has \code{c("formula_lhs_\var{CLASS}",
   "formula_lhs")} prepended to it, where \var{CLASS} is the class
   of the LHS value or \code{basis}. If LHS or \code{basis} has multiple
   classes, they are all prepended; if there is no LHS or \code{basis},
   \code{c("formula_lhs_", "formula_lhs")} is.

\item \code{\link{simulate}} generic is evaluated on the new \code{object}, with all
   arguments passed on, excluding \code{basis}; if \code{newdata} or \code{data}
   are missing, they too are not passed on. The evaluation takes
   place in the parent's environment.
}
A "method" to receive a formula whose LHS evaluates to \var{CLASS}
can therefore be implemented by a function
\code{simulate.formula_lhs_\var{CLASS}()}. This function can expect a
\code{\link{formula}} object, with additional attribute \code{.Basis} giving the
evaluated LHS (so that it does not need to be evaluated again).
}

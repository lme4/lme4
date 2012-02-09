\name{Nelder_Mead}
\alias{Nelder_Mead}
\title{Nelder-Mead optimization}
\usage{
  Nelder_Mead(ff, x0, xst, xt, lower = rep.int(-Inf, n),
    upper = rep.int(Inf, n), control = list())
}
\arguments{
  \item{ff}{a function of one numeric vector argument
  returning a numeric scalar}

  \item{lower}{numeric vector of lower bounds - elements
  may be \code{-Inf}.}

  \item{upper}{numeric vector of upper bounds - elements
  may be \code{Inf}.}

  \item{xst}{numeric vector of initial step sizes to
  establish the simplex - all elements must be non-zero.}

  \item{x0}{numeric vector of starting values for the
  parameters.}

  \item{xt}{numeric vector of tolerances on the
  parameters.}

  \item{control}{a named list of control settings.
  Possible settings are \describe{ \item{iprint}{numeric
  scalar - frequency of printing evaluation information.
  Defaults to 0 indicating no printing.}
  \item{maxfun}{numeric scalar - maximum number of function
  evaluations allowed.} \item{FtolAbs}{numeric scalar -
  absolute tolerance on change in function values}
  \item{FtolRel}{numeric scalar - relative tolerance on
  change in function values} \item{XtolRel}{numeric scalar
  - relative tolerance on change in parameter values}
  \item{MinfMax}{numeric scalar - maximum value of the
  minimum} }}
}
\value{
  a list with 4 components \item{fval}{numeric scalar - the
  minimum function value achieved} \item{par}{numeric
  vector - the value of \code{x} providing the minimum}
  \item{ierr}{integer scalar - error code}
  \item{control}{list - the list of control settings after
  substituting for defaults}
}
\description{
  Nelder-Mead optimization of parameters that may be
  subject to box constraints
}


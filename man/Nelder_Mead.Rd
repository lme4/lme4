\name{Nelder_Mead}
\alias{Nelder_Mead}
\title{Nelder-Mead optimization}
\usage{
  Nelder_Mead(fn, par, lower = rep.int(-Inf, n),
    upper = rep.int(Inf, n), control = list())
}
\arguments{
  \item{fn}{a function of one numeric vector argument
  returning a numeric scalar}

  \item{par}{numeric vector of starting values for the
  parameters.}

  \item{lower}{numeric vector of lower bounds - elements
  may be \code{-Inf}.}

  \item{upper}{numeric vector of upper bounds - elements
  may be \code{Inf}.}

  \item{control}{a named list of control settings.
  Possible settings are \describe{ \item{iprint}{numeric
  scalar - frequency of printing evaluation information.
  Defaults to 0 indicating no printing.}
  \item{maxfun}{numeric scalar - maximum number of function
  evaluations allowed (default:10000).}
  \item{FtolAbs}{numeric scalar - absolute tolerance on
  change in function values (default: 1e-5)}
  \item{FtolRel}{numeric scalar - relative tolerance on
  change in function values (default:1e-15)}
  \item{XtolRel}{numeric scalar - relative tolerance on
  change in parameter values (default: 1e-7)}
  \item{MinfMax}{numeric scalar - maximum value of the
  minimum (default: .Machine$double.xmin)}
  \item{xst}{numeric vector of initial step sizes to
  establish the simplex - all elements must be non-zero
  (default: rep(0.02,length(par)))} \item{xt}{numeric
  vector of tolerances on the parameters (default:
  xst*5e-4)} \item{verbose}{numeric value: 0=no printing,
  1=print every 20 evaluations, 2=print every 10
  evalutions, 3=print every evaluation.  Sets
  \sQuote{iprint}, if specified, but does not override it.}
  }}
}
\value{
  a list with 4 components \item{fval}{numeric scalar - the
  minimum function value achieved} \item{par}{numeric
  vector - the value of \code{x} providing the minimum}
  \item{ierr}{integer scalar - error code (see below)}
  \item{control}{list - the list of control settings after
  substituting for defaults}
}
\description{
  Nelder-Mead optimization of parameters that may be
  subject to box constraints
}
\note{
  Return codes: \describe{ \item{-4}{\code{nm_evals}:
  maximum evaluations reached} \item{-3}{\code{nm_forced}:
  ?} \item{-2}{\code{nm_nofeasible}: cannot generate a
  feasible simplex} \item{-1}{\code{nm_x0notfeasible}:
  initial x is not feasible (?)} }
}


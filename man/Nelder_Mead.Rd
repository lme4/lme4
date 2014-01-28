\name{NelderMead}
\alias{Nelder_Mead}
\title{Nelder-Mead Optimization of Parameters, Possibly (Box) Constrained}
\description{
  Nelder-Mead optimization of parameters, allowing
  optimization subject to box constraints (contrary to the default,
  \code{method = "Nelder-Mead"}, in \R's \code{\link{optim}()}), and
  using reverse communications.
}
\usage{
Nelder_Mead(fn, par, lower = rep.int(-Inf, n), upper = rep.int(Inf, n),
            control = list())
}
\arguments{
  \item{fn}{a \code{\link{function}} of a single numeric vector argument
    returning a numeric scalar.}

  \item{par}{numeric vector of starting values for the parameters.}

  \item{lower}{numeric vector of lower bounds (elements may be \code{-Inf}).}
  \item{upper}{numeric vector of upper bounds (elements may be \code{Inf}).}

  \item{control}{a named list of control settings. Possible settings are
    \describe{
      \item{iprint}{numeric scalar - frequency of printing evaluation information.
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
	(default: rep(0.02,length(par)))}
      \item{xt}{numeric vector of tolerances on the parameters (default: xst*5e-4)}
      \item{verbose}{numeric value: 0=no printing,
	1=print every 20 evaluations, 2=print every 10
	evalutions, 3=print every evaluation.  Sets
	\sQuote{iprint}, if specified, but does not override it.}
      \item{warnOnly}{a logical indicating if non-convergence (codes
	-1,-2,-3) should not \code{\link{stop}(.)}, but rather only call
	\code{\link{warning}} and return a result which might inspected.
	Defaults to \code{FALSE}, i.e., stop on non-convergence.}
    }
  }
}
\value{
  a \code{\link{list}} with components
  \item{fval}{numeric scalar - the minimum function value achieved}
  \item{par}{numeric vector - the value of \code{x} providing the minimum}
  \item{convergence}{integer valued scalar, if not \code{0}, an error code:
    \describe{
      \item{-4}{\code{nm_evals}: maximum evaluations reached}
      \item{-3}{\code{nm_forced}: ?}
      \item{-2}{\code{nm_nofeasible}: cannot generate a feasible simplex}
      \item{-1}{\code{nm_x0notfeasible}: initial x is not feasible (?)}
      \item{0}{successful convergence}
    }
  }
  \item{message}{a string specifying the kind of convergence.}
  \item{control}{the \code{\link{list}} of control settings after
    substituting for defaults.}
  \item{feval}{the number of function evaluations.}
}
\seealso{
  The \code{\linkS4class{NelderMead}} class definition and generator
  function.
}
\examples{
fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
p0 <- c(-1.2, 1)

oo  <- optim(p0, fr) ## also uses Nelder-Mead by default
o.  <- Nelder_Mead(fr, p0)
o.1 <- Nelder_Mead(fr, p0, control=list(verbose=1))# -> some iteration output
stopifnot(identical(o.[1:4], o.1[1:4]),
          all.equal(o.$par, oo$par, tolerance=1e-3))# diff: 0.0003865

%%## but this shows that something "does not work"
o.2 <- Nelder_Mead(fr, p0, control=list(verbose=3, XtolRel=1e-15, FtolAbs= 1e-14))
all.equal(o.2[-5],o.1[-5], tolerance=1e-15)# TRUE, unexpectedly
}
\keyword{classes}


\name{findbars}
\alias{findbars}
\alias{nobars}
\alias{subbars}
\title{Utilities for random-effects terms in formulas}
\description{
  \describe{
    \item{\code{findbars}}{Extracts the pairs of expressions that are
      separated by the vertical bar operator, expanding the slash operator
      as necessary.}
    \item{\code{subbars}}{\dQuote{Sub[stitute] Bars}: Replaces \code{|}
      operators with \code{+} operators.}
    \item{\code{nobars}}{Omits the random-effects expressions.}
  }
}
\usage{
findbars(term)
nobars(term)
subbars(term)
}
\arguments{
  \item{term}{a mixed-effects model formula}
}
\section{Note}{All these functions are called recursively on individual
  terms in the model, which is why the argument is called \code{term}
  and not a name like \code{form}, indicating a formula.}
\value{
  \describe{
    \item{\code{findbars}}{A list of pairs of expressions that were
      separated by vertical bars.}
   \item{\code{nobars}}{The expression with the random effects removed.}
   \item{\code{subbars}}{The expression with all \code{|} operators
     replaced by \code{+} operators.}
   }
}
\details{
  These three functions extend the range of formulas allowed by
  the \code{\link{model.frame}} and \code{\link{model.matrix}} functions
  to those with random-effects terms.

  \code{subbars} is applied to the model formula before it is passed to
  the \code{\link{model.frame}} function.  (Note: It is likely that a
  trivial extension to the \code{model.frame} function would make this
  unnecessary.)
  
  The result of \code{nobars} is passed to \code{\link{model.matrix}} to
  evaluate the fixed-effects model matrix and the result of
  \code{findbars} is used to construct the random-effects model matrix
  and the relative covariance factor.
}
\seealso{
 \code{\link{formula}},
 \code{\link{model.frame}}, \code{\link{model.matrix}}.
}
\examples{
findbars(f1 <- Reaction ~ Days + (Days|Subject))
## => list( Days | Subject )
findbars(y ~ Days + (1|Subject) + (0+Days|Subject))
## => list of length 2:  list ( 1 | Subject ,  0+Days|Subject)

nobars(f1)  ## Reaction ~ Days
subbars(f1) ## Reaction ~ Days + (Days + Subject)
\dontshow{
stopifnot(identical(findbars(f1),
                    list(expression(Days | Subject)[[1]])))
}
}
\keyword{models}
\keyword{utilities}

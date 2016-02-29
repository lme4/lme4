\name{expandDoubleVerts}
\alias{expandDoubleVerts}
\alias{||}
\title{Expand terms with \code{'||'} notation into separate \code{'|'} terms}
\usage{
  expandDoubleVerts(term)
}
\arguments{
  \item{term}{a mixed-model formula}
}
\value{
 the modified term
}
\description{
  From the right hand side of a formula for a mixed-effects model,
  expand terms with the double vertical bar operator
  into separate, independent random effect terms.
}
\note{
  Because \code{||} works at the level of formula parsing, it
  has no way of knowing whether a variable is a factor. It
  just takes the terms within a random-effects term and literally splits them
  into the intercept and separate no-intercept terms,
  e.g. \code{(1+x+y|f)} would be split into \code{(1|f) + (0+x|f) + (0+y|f)}.
  However, \code{||} will fail to break up factors into separate terms;
  the \code{\link{dummy}} function can be useful in this case, although
  it is not as convenient as \code{||}.
}
\examples{
  m <- ~ x + (x || g)
  expandDoubleVerts(m)
  set.seed(101)
  dd <- expand.grid(f=factor(letters[1:3]),g=factor(1:200),rep=1:3)
  dd$y <- simulate(~f + (1|g) + (0+dummy(f,"b")|g) + (0+dummy(f,"c")|g),
          newdata=dd,
          newparams=list(beta=rep(0,3),
                         theta=c(1,2,1),
                         sigma=1),
          family=gaussian)[[1]]
  m1 <- lmer(y~f+(f|g),data=dd)
  VarCorr(m1)
  m2 <- lmer(y~f+(1|g) + (0+dummy(f,"b")|g) + (0+dummy(f,"c")|g),
               data=dd)
  VarCorr(m2)
}
\seealso{
  \code{\link{formula}}, \code{\link{model.frame}},
  \code{\link{model.matrix}}, \code{\link{dummy}}.

  Other utilities: \code{\link{mkRespMod}},
  \code{\link{mkReTrms}}, \code{\link{nlformula}},
  \code{\link{nobars}}, \code{\link{subbars}}
}
\keyword{models}
\keyword{utilities}

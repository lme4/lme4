\name{glmer.nb}
\title{Fitting Negative Binomial GLMMs}
\alias{glmer.nb}
\alias{negative.binomial}% re-exported, needed e.g. in update()
\concept{GLMM}
\description{
  Fits a generalized linear mixed-effects model (GLMM) for the negative
  binomial family, building on \code{\link{glmer}}, and initializing via
  \code{\link[MASS]{theta.ml}} from \pkg{MASS}.
}
\usage{
glmer.nb(\dots, interval = log(th) + c(-3, 3),
         tol = 5e-5, verbose = FALSE, nb.control = NULL,
         initCtrl = list(limit = 20, eps = 2*tol, trace = verbose,
                         theta = NULL))
}
\arguments{
  \item{\dots}{arguments as for \code{glmer(.)} such as \code{formula},
    \code{data}, \code{control}, etc, but \emph{not} \code{family}!}
  \item{interval}{interval in which to start the optimization.  The
    default is symmetric on log scale around the initially estimated theta.}
  \item{tol}{tolerance for the optimization via \code{\link{optimize}}.}
  \item{verbose}{\code{\link{logical}} indicating how much
    progress information should be printed during the optimization.  Use
    \code{verbose = 2} (or larger) to enable \code{verbose=TRUE} in the
    \code{\link{glmer}()} calls.}
  \item{nb.control}{optional \code{\link{list}}, like \code{\link{glmerControl}()},
    used in \code{\link{refit}(*, control = control.nb)} during the
    optimization.}
  \item{initCtrl}{(\emph{\bold{experimental}, do not rely on this}:) a
    \code{\link{list}} with named components as in the default, passed to
    \code{\link[MASS]{theta.ml}} (package \pkg{MASS}) for the initial
    value of the negative binomial parameter \code{theta}.
    May also include a \code{theta} component, in which case the
    initial estimation step is skipped}
}
\value{
  An object of class \code{glmerMod}, for which many
  methods are available (e.g. \code{methods(class="glmerMod")}), see
  \code{\link{glmer}}.
}
\note{For historical reasons, the shape parameter of the negative
  binomial and the random effects parameters in our (G)LMM models are
  both called \code{theta} (\eqn{\theta}), but are unrelated here.

  The negative binomial \eqn{\theta} can be extracted from a fit
  \code{g <- glmer.nb()} by \code{\link{getME}(g, "glmer.nb.theta")}.

  Parts of \code{glmer.nb()} are still experimental and methods are
  still missing or suboptimal.  In particular, there is no inference
  available for the dispersion parameter \eqn{\theta}, yet.

  To fit a negative binomial model with \emph{known} overdispersion
  parameter (e.g. as part of a model comparison exercise, use
  \code{glmer} with the \code{\link[MASS]{negative.binomial}} family from the
  \code{MASS} package, e.g.
  \code{glmer(...,family=MASS::negative.binomial(theta=1.75))}.
  
}
\seealso{
  \code{\link{glmer}}; from package \pkg{MASS},
  \code{\link[MASS]{negative.binomial}} (which we re-export currently) and
  \code{\link[MASS]{theta.ml}}, the latter for initialization of
  optimization.

  The \sQuote{Details} of \code{\link{pnbinom}} for the definition of
  the negative binomial distribution.
}
\examples{
set.seed(101)
dd <- expand.grid(f1 = factor(1:3),
                  f2 = LETTERS[1:2], g=1:9, rep=1:15,
          KEEP.OUT.ATTRS=FALSE)
summary(mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2))))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)
str(dd)
require("MASS")## and use its glm.nb() - as indeed we have zero random effect:
\dontrun{
m.glm <- glm.nb(y ~ f1*f2, data=dd, trace=TRUE)
summary(m.glm)
m.nb <- glmer.nb(y ~ f1*f2 + (1|g), data=dd, verbose=TRUE)
m.nb
## The neg.binomial theta parameter:
getME(m.nb, "glmer.nb.theta")
LL <- logLik(m.nb)
## mixed model has 1 additional parameter (RE variance)
stopifnot(attr(LL,"df")==attr(logLik(m.glm),"df")+1)
plot(m.nb, resid(.) ~ g)# works, as long as data 'dd' is found
}
}
\keyword{models}

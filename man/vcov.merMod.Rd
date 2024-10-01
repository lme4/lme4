\name{vcov.merMod}
\title{Covariance matrix of estimated parameters}
\description{
  Compute the variance-covariance matrix of estimated paramers.
  Optionally also computes correlations, or the full (joint)
  covariance matrix of the fixed-effect coefficients and the conditional
  modes of the random effects.}
\alias{vcov.merMod}
\alias{vcov.summary.merMod}
\usage{
\S3method{vcov}{merMod}(object, correlation = TRUE, sigm = sigma(object),
    use.hessian = NULL, full = FALSE, \dots)
}
\arguments{
  \item{object}{an \R object of class \code{\linkS4class{merMod}}, i.e.,
    as resulting from \code{\link{lmer}()}, or \code{\link{glmer}()},
    etc.}
  \item{correlation}{(logical)
    indicates whether the correlation matrix as well as
    the variance-covariance matrix is desired
  }
  \item{sigm}{the residual standard error; by default \code{\link{sigma}(object)}.}
  \item{use.hessian}{(logical) indicates whether to use the
    finite-difference Hessian of the deviance function to compute
    standard errors of the fixed effects. See \emph{Details}.
  }
  \item{full}{return the 'full' covariance matrix, i.e. the joint
    covariance matrix of the conditional distribution of conditional
    modes (as in \code{getME(., "b")}) and fixed-effect
    parameters. (\code{correlation} and \code{use.hessian} are
    \emph{ignored} in this case.)

    Note that this option may be slow for models with large numbers of
  random-effect levels!
  }
  \item{\dots}{extra arguments for method compatibility (ignored)}
}
\details{
  When \code{use.hessian = FALSE}, the code estimates the
  covariance matrix based on internal information about the inverse of
  the model matrix (see \code{\link{getME}(.,"RX")}). This is exact for
  linear mixed models, but approximate for GLMMs.
  The default is to use the Hessian whenever the
  fixed effect parameters are arguments to the deviance
  function (i.e. for GLMMs with \code{nAGQ>0}), and to use
  \code{\link{getME}(.,"RX")} whenever the fixed effect parameters are
  profiled out (i.e. for GLMMs with \code{nAGQ==0} or LMMs).

  \code{use.hessian=FALSE} is backward-compatible with older versions
  of \code{lme4}, but may give less accurate SE estimates when the
  estimates of the fixed-effect (see \code{\link{getME}(.,"beta")})
  and random-effect (see \code{\link{getME}(.,"theta")}) parameters
  are correlated.

  However, \code{use.hessian=TRUE} is not always more accurate: for some
  numerically unstable fits, the approximation using \code{RX} is
  actually more reliable (because the Hessian has to be computed by
  a finite difference approximation, which is also error-prone): see
  e.g. \url{https://github.com/lme4/lme4/issues/720}

}
    
\examples{
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
(v1 <- vcov(fm1))
v2 <- vcov(fm1, correlation = TRUE)
# extract the hidden 'correlation' entry in @factors
as(v2, "corMatrix")
v3 <- vcov(gm1)
v3X <- vcov(gm1, use.hessian  = FALSE)
all.equal(v3, v3X)
## full correlatiom matrix
cv <- vcov(fm1, full = TRUE)
image(cv, xlab = "", ylab = "",
      scales = list(y = list(labels = rownames(cv)),
                    at = seq(nrow(cv)),
                    x = list(labels = NULL)))
}
\value{a covariance matrix (sparse when \code{full=TRUE})}
\keyword{models}






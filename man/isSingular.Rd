\name{isSingular}
\alias{isSingular}
\alias{boundary}
\title{
  Test fitted model for (near) singularity
}
\description{
  Evaluates whether a fitted mixed model is singular, i.e.
  the parameters are on the boundary of the feasible
  parameter space: variances of one or more linear combinations of effects
  are (close to) zero.
}
\usage{
isSingular(x, tol = 1e-05)
}
\arguments{
  \item{x}{
    a fitted \code{merMod} object (result of \code{lmer} or \code{glmer})
}
  \item{tol}{
    numerical tolerance for detecting singularity
}
}
\details{
  Complex mixed-effect models (i.e., those with a
  large number of variance-covariance parameters) frequently result in
  \emph{singular} fits, i.e. estimated variance-covariance matrices with
  less than full rank. Less technically, this means that some
  "dimensions" of the variance-covariance matrix have been estimated as
  exactly zero. For scalar random effects such as intercept-only models,
  or 2-dimensional random effects such as intercept+slope models,
  singularity is relatively easy to detect because it leads to
  random-effect variance estimates of (nearly) zero, or estimates of
  correlations that are (almost) exactly -1 or 1. However, for more
  complex models (variance-covariance matrices of dimension >=3)
  singularity can be hard to detect; models can often be singular
  without any of their individual variances being close to zero
  or correlations being close to +/-1.

  This function performs a simple test to determine whether any of the
  random effects covariance matrices of a fitted model are singular.
  The \code{\link{rePCA}} method provides more detail about the
  singularity pattern, showing the standard deviations
  of orthogonal variance components and the mapping from
  variance terms in the model to orthogonal components
  (i.e., eigenvector/rotation matrices).
  
  While singular models are statistically well defined (it is
  theoretically sensible for the true maximum likelihood estimate
  to correspond to a singular fit), there are real concerns
  that (1) singular fits correspond to overfitted models that may
  have poor power; (2) chances of numerical problems and mis-convergence
  are higher for singular models (e.g. it may be computationally
  difficult to compute profile confidence intervals for such
  models); (3) standard inferential procedures such as Wald statistics
  and likelihood ratio tests may be inappropriate.
  
  There is not yet consensus about how to deal with singularity, or
  more generally to choose which random-effects specification (from
  a range of choices of varying complexity) to use.  Some proposals
  include:

  \itemize{
    \item avoid fitting overly complex models in the first place,
    i.e. design experiments/restrict models \emph{a priori} such that the
    variance-covariance matrices can be estimated precisely enough
    to avoid singularity (Matuschek et al 2017)
    \item use some form of model selection to choose a model that
    balances predictive accuracy and overfitting/type I error
    (Bates et al 2015, Matuschek et al 2017)
    \item \dQuote{keep it maximal}, i.e. fit the most complex
    model consistent with the experimental design, removing only
    terms required to allow a non-singular fit (Barr et al. 2013),
    or removing further terms based on p-values or AIC
    \item use a partially Bayesian method
    that produces maximum \emph{a posteriori} (MAP)
    estimates using \emph{regularizing} priors to force
    the estimated random-effects variance-covariance matrices
    away from singularity (Chung et al 2013, \pkg{blme} package)
    \item use a fully Bayesian method that both regularizes the
    model via informative priors and gives estimates
    and credible intervals for all parameters that average
    over the uncertainty in the random effects parameters
    (Gelman and Hill 2006, McElreath 2015;
    \pkg{MCMCglmm}, \pkg{rstanarm} and \pkg{brms} packages)
  }  
}
\value{
  a logical value 
}
\seealso{\code{\link{getME}}, \code{\link{rePCA}}}
\references{
  \itemize{
    \item Dale J. Barr, Roger Levy, Christoph Scheepers, and Harry J. Tily. Random effects structure for confirmatory hypothesis testing: Keep it maximal.
    Journal of Memory and Language, 68(3):255–278, April 2013.
    \item Douglas Bates, Reinhold Kliegl, Shravan Vasishth, and Harald Baayen.
    Parsimonious Mixed Models. arXiv:1506.04967 [stat], June 2015. arXiv:
    1506.04967.
    \item Yeojin Chung, Sophia Rabe-Hesketh, Vincent Dorie, Andrew Gelman, and
    Jingchen Liu. A nondegenerate penalized likelihood estimator for variance
    parameters in multilevel models. Psychometrika, pages 1–25, 2013.
    \item Andrew Gelman and Jennifer Hill. Data Analysis Using Regression and
    Multilevel/Hierarchical Models. Cambridge University Press, Cambridge,
    England, 2006.
    \item Hannes Matuschek, Reinhold Kliegl, Shravan Vasishth, Harald Baayen, and
    Douglas Bates. Balancing type I error and power in linear mixed models.
    Journal of Memory and Language, 94:305–315, 2017.
    \item Richard McElreath. Statistical Rethinking: A Bayesian Course with Examples in R and Stan. Chapman and Hall/CRC, Boca Raton, December
    2015
  }
}
\keyword{utilities}

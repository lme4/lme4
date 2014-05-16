\name{isREML}
\alias{isGLMM}
\alias{isLMM}
\alias{isNLMM}
\alias{isREML}
\alias{isGLMM.merMod}
\alias{isLMM.merMod}
\alias{isNLMM.merMod}
\alias{isREML.merMod}
\title{Check characteristics of models}
\usage{
  isREML(x, ...)

  isLMM(x, ...)

  isNLMM(x, ...)

  isGLMM(x, ...)
}
\arguments{
  \item{x}{a fitted model.}

  \item{...}{additional, optional arguments.  (None are
  used in the merMod methods)}
}
\value{
  a logical value
}
\description{
  Check characteristics of models: whether a model fit
  corresponds to a linear (LMM), generalized linear (GLMM),
  or nonlinear (NLMM) mixed model, and whether a linear
  mixed model has been fitted by REML or not
  (\code{isREML(x)} is always \code{FALSE} for GLMMs and
  NLMMs).
}
\details{
  These are generic functions.  At present the only methods
  are for mixed-effects models of class
  \code{\linkS4class{merMod}}.
}
\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
             Orange, start = c(Asym = 200, xmid = 725, scal = 350))

isLMM(fm1)
isGLMM(gm1)
## check all :
is.MM <- function(x) c(LMM = isLMM(x), GLMM= isGLMM(x), NLMM= isNLMM(x))
stopifnot(cbind(is.MM(fm1), is.MM(gm1), is.MM(nm1))
	  == diag(rep(TRUE,3)))
}
\seealso{
  getME
}


\name{mkRespMod}
\alias{mkRespMod}
\title{Create an lmerResp, glmResp or nlsResp instance}
\usage{
  mkRespMod(fr, REML = NULL, family = NULL, nlenv = NULL,
    nlmod = NULL, ...)
}
\arguments{
  \item{fr}{a model frame}

  \item{REML}{logical scalar, value of REML for an lmerResp
  instance}

  \item{family}{the optional glm family (glmResp only)}

  \item{nlenv}{the nonlinear model evaluation environment
  (nlsResp only)}

  \item{nlmod}{the nonlinear model function (nlsResp only)}

  \item{...}{where to look for response information if
  \code{fr} is missing.  Can contain a model response,
  \code{y}, offset, \code{offset}, and weights,
  \code{weights}.}
}
\value{
  an lmerResp or glmResp or nlsResp instance
}
\description{
  Create an lmerResp, glmResp or nlsResp instance
}
\seealso{
  Other utilities: \code{\link{findbars}},
  \code{\link{mkReTrms}}, \code{\link{nlformula}},
  \code{\link{nobars}}, \code{\link{subbars}}
}


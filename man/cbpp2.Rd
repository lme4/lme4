\docType{data}
\name{cbpp2}
\alias{cbpp2}
\title{Contagious bovine pleuropneumonia (extended version)}
\description{
  This dataset is an extension of \code{\link[lme4]{cbpp}}, which describes the 
  serological incidence of CBPP (Contagious bovine pleuropneumonia) in zebu 
  cattle during a follow-up survey implemented in 15 commercial herds located in 
  the Boji district of Ethiopia. There are two extra covariates: \code{treatment} 
  and \code{avg_size}.
}
\format{A data frame with 56 observations on the following 6 variables.
  \describe{
    \item{\code{herd}}{A factor identifying the herd (1 to 15).}
    \item{\code{treatment}}{A factor referring to the control measure used
      to manage CBPP.
      \itemize{
        \item{\code{Complete} = complete isolation or antibiotic treatment,}
        \item{\code{Partial/null} = partial/null isolation and no antibiotic
          treatment,}
        \item{\code{Unknown} = strategy remained.}
      }
    }
    \item{\code{avg_size}}{The average number of animals housed in a della 
      (a temporary paddock used for holding cattle on the farm).}
    \item{\code{incidence}}{The number of new serological cases for a
      given herd and time period.}
    \item{\code{size}}{A numeric vector describing herd size at the
      beginning of a given time period.}
    \item{\code{period}}{A factor with levels \code{1} to \code{4}.}
  }}
\source{
  Lesnoff, M., Laval, G., Bonnet, P., Abdicho, S.,
  Workalemahu, A., Kifle, D., Peyraud, A., Lancelot, R.,
  Thiaucourt, F. (2004) Within-herd spread of contagious
  bovine pleuropneumonia in Ethiopian highlands.
  \emph{Preventive Veterinary Medicine} \bold{64}, 27--40.
}
\details{
  The description here is identical to the \code{cbpp} dataset:
  Contagious bovine pleuropneumonia (CBPP) is a major disease of cattle in 
  Africa, caused by a mycoplasma.  The goal of the survey was to study the
  within-herd spread of CBPP in newly infected herds. Blood samples were 
  quarterly collected from all animals of these herds to determine their CBPP 
  status.  These data were used to compute the serological incidence of CBPP
  (new cases occurring during a given time period).  Some data are missing 
  (lost to follow-up). Serological status was determined using a competitive
  enzyme-linked immuno-sorbent assay (cELISA).
}
\seealso{
  The shorter version, \code{\link[lme4]{cbpp}}.
}
\examples{
## Fitting the model 
gm1 <- glmer(incidence/size ~ period + treatment + avg_size + (1 | herd),
             family = binomial,
             data = cbpp2, weights = size,
             control = glmerControl(optimizer="bobyqa"))
## Adding an observation-level random effect
cbpp2 <- transform(cbpp2,obs=factor(seq(nrow(cbpp2))))
## Herd and observation-level REs (below causes singular fit issues)
gm2 <- update(gm1,.~.+(1|obs)) 
## observation-level REs only (no singular fit issue)
gm3 <- update(gm1,.~.-(1|herd)+(1|obs)) 
}
\keyword{datasets}


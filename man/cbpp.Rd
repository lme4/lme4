\docType{data}
\name{cbpp}
\alias{cbpp}
\title{Contagious bovine pleuropneumonia}
\format{A data frame with 56 observations on the following 4 variables.
  \describe{
    \item{\code{herd}}{A factor identifying the herd (1 to 15).}
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
\description{
  Contagious bovine pleuropneumonia (CBPP) is a major
  disease of cattle in Africa, caused by a mycoplasma.
  This dataset describes the serological incidence of CBPP
  in zebu cattle during a follow-up survey implemented in
  15 commercial herds located in the Boji district of
  Ethiopia.  The goal of the survey was to study the
  within-herd spread of CBPP in newly infected herds. Blood
  samples were quarterly collected from all animals of
  these herds to determine their CBPP status.  These data
  were used to compute the serological incidence of CBPP
  (new cases occurring during a given time period).  Some
  data are missing (lost to follow-up).
}
\details{
  Serological status was determined using a competitive
  enzyme-linked immuno-sorbent assay (cELISA).
}
\examples{
## response as a matrix
(m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             cbpp, binomial, nAGQ=25L))
## response as a vector of probabilities and usage of argument "weights"
m1p <- glmer(incidence / size ~ period + (1 | herd), weights = size,
             cbpp, binomial, nAGQ=25L)
## Confirm that these are equivalent:
stopifnot(all.equal(fixef(m1), fixef(m1p), tol = 1e-5),
          all.equal(ranef(m1), ranef(m1p), tol = 1e-5),
          TRUE)
## Can this section be moved to a test file?  I don't think it belongs in an example. DB
for(m in c(m1, m1p)) {
    cat("-------\\n\\nCall: ",
        paste(format(getCall(m)), collapse="\\n"), "\\n")
    print(logLik(m)); cat("AIC:", AIC(m), "\\n") ; cat("BIC:", BIC(m),"\\n")
}
stopifnot(all.equal(logLik(m1), logLik(m1p), tol = 1e-5),
          all.equal(AIC(m1),    AIC(m1p),    tol = 1e-5),
          all.equal(BIC(m1),    BIC(m1p),    tol = 1e-5))

## GLMM with individual-level variability (accounting for overdispersion)
cbpp$obs <- 1:nrow(cbpp)
(m2 <- glmer(cbind(incidence, size - incidence) ~ period +
    (1 | herd) +  (1|obs),
              family = binomial, data = cbpp))
}
\keyword{datasets}


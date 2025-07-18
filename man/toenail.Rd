\name{toenail}
\alias{toenail}
\docType{data}
\title{
  Toenail Onychomycosis Data from Dermatophyte Infections
}
\description{
  Results of a clinical trial by De Backer et al. (1998) to compare 
  two oral antifungal treatments: itraconazole and terbinafine to treat
  toenail onychomycosis causes by dermatophytes.
}
\usage{data("toenail")}
\format{
  A data frame with 1908 observations on the following 8 variables.
  \describe{
    \item{\code{patientID}}{a unique identifier for each patient in the trial.}
    \item{\code{outcome}}{a factor indicating the severity of the patient's 
      onycholysis (separation of the nail plate from the nail bed), with levels 
      \code{none or mild} and \code{moderate or severe}.}
    \item{\code{treatment}}{a factor indicating the treatment used with levels
      \code{itraconazole} and \code{terbinafine}.}
    \item{\code{month}}{exact timing of visit in months.}
    \item{\code{visit}}{An integer represent the sequential visit number for
      each subject.}
  }
}
\details{
  The data in \code{lme4} was copied from the toenail data found in the 
  \code{HSAUR3} package. From the dataset documentation from the 
  \code{HSAUR3} package: De Backer et al. (1998) describe a clinical trial to 
  compare two competing oral antifungal treatments for toenail infection 
  (dermatophyte onychomycosis). A total of 378 patients were randomly allocated 
  into two treatment groups, one group receiving 250mg per day of terbinafine 
  and the other group 200mg per day of itraconazole. Patients were evaluated at 
  seven visits, intended to be at weeks 0, 4, 8, 12, 24, 36, and 48 for the 
  degree of separation of the nail plate from the nail bed (onycholysis) 
  dichotomized into \code{moderate or severe} and \code{none or mild}. But 
  patients did not always arrive exactly at the scheduled time and the exact 
  time in months that they did attend was recorded. The data is not balanced 
  since not all patients attended for all seven planned visits.
}
\source{
  M. D. Backer and C. D. Vroey and E. Lesaffre and I. Scheys and P. D. Keyser. 
}
\references{
  M. D. Backer and C. D. Vroey and E. Lesaffre and I. Scheys and P. D. Keyser 
  (1998), Twelve weeks of continuous oral therapy for toenail onychomycosis 
  caused by dermatophytes: A double-blind comparative trial of terbinafine 
  250 mg/day versus itraconazole 200 mg/day. \emph{Journal of the American 
  Academy of Dermatology}, \bold{38}, S57--S63; 
  \doi{https://doi.org/10.1016/S0190-9622(98)70486-4}
}
\examples{
toemod <- glmer(outcome ~ time*treatment + (1 | patientID), 
                data = toenail, family = binomial(link = "logit"))
}
\keyword{datasets}

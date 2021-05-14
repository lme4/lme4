\docType{data}
\name{sleepstudy}
\alias{sleepstudy}
\title{Reaction times in a sleep deprivation study}
\format{A data frame with 180 observations on the following 3 variables.
  \describe{
    \item{\code{Reaction}}{Average reaction time (ms)}
    \item{\code{Days}}{Number of days of sleep deprivation}
    \item{\code{Subject}}{Subject number on which the observation was made.}
  }}
\description{
  The average reaction time per day (in milliseconds) for subjects in a sleep
  deprivation study.

  Days 0-1 were adaptation and training (T1/T2), day 2 was baseline (B);
  sleep deprivation started after day 2.

}
\details{
  These data are from the study described in Belenky et al.
  (2003), for the most sleep-deprived group (3 hours time-in-bed)
  and for the first 10 days of the study, up to the recovery period.
  The original study analyzed speed (1/(reaction time)) and treated
  day as a categorical rather than a continuous predictor.
}
\examples{
str(sleepstudy)
require(lattice)
xyplot(Reaction ~ Days | Subject, sleepstudy, type = c("g","p","r"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "Days of sleep deprivation",
       ylab = "Average reaction time (ms)", aspect = "xy")
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, subset=Days>=2))
## independent model
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy, subset=Days>=2))
}
\references{
  Gregory Belenky, Nancy J. Wesensten, David R. Thorne,
  Maria L. Thomas, Helen C. Sing, Daniel P. Redmond,
  Michael B. Russo and Thomas J. Balkin (2003) Patterns of
  performance degradation and restoration during sleep
  restriction and subsequent recovery: a sleep
  dose-response study. \emph{Journal of Sleep Research}
  \bold{12}, 1--12.
}
\keyword{datasets}


\name{grouseticks}
\alias{grouseticks}
\alias{grouseticks_agg}
\docType{data}
\title{
  Data on red grouse ticks from Elston et al. 2001
}
\description{
  Number of ticks on the heads of red grouse chicks
  sampled in the field (\code{grouseticks})
  and an aggregated version (\code{grouseticks_agg}); see original source for more details
}
\usage{data(grouseticks)}
\format{
  \describe{
    \item{\code{INDEX}}{(factor) chick number (observation level)}
    \item{\code{TICKS}}{number of ticks sampled}
    \item{\code{BROOD}}{(factor) brood number}
    \item{\code{HEIGHT}}{height above sea level (meters)}
    \item{\code{YEAR}}{year (-1900)}
    \item{\code{LOCATION}}{(factor) geographic location code}
    \item{\code{cHEIGHT}}{centered height, derived from \code{HEIGHT}}
    \item{\code{meanTICKS}}{mean number of ticks by brood}
    \item{\code{varTICKS}}{variance of number of ticks by brood}
  }
}
\source{
  Robert Moss, via David Elston
}
\details{\code{grouseticks_agg} is just a brood-level
  aggregation of the data}
\references{
  Elston, D. A., R. Moss, T. Boulinier, C. Arrowsmith, and
X. Lambin. 2001. "Analysis of Aggregation, a Worked Example: Numbers of
Ticks on Red Grouse Chicks." Parasitology 122 (05):
563-569. doi:10.1017/S0031182001007740.
\url{http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=82701}.
}
\examples{
data(grouseticks)
## Figure 1a from Elston et al
par(las=1,bty="l")
tvec <- c(0,1,2,5,20,40,80)
pvec <- c(4,1,3)
with(grouseticks_agg,plot(1+meanTICKS~HEIGHT,
                  pch=pvec[factor(YEAR)],
                  log="y",axes=FALSE,
                  xlab="Altitude (m)",
                  ylab="Brood mean ticks"))
axis(side=1)
axis(side=2,at=tvec+1,label=tvec)
box()
abline(v=405,lty=2)
## Figure 1b
with(grouseticks_agg,plot(varTICKS~meanTICKS,
                  pch=4,
                  xlab="Brood mean ticks",
                  ylab="Within-brood variance"))
curve(1*x,from=0,to=70,add=TRUE)
## Model fitting
form <- TICKS~YEAR+HEIGHT+(1|BROOD)+(1|INDEX)+(1|LOCATION)
(full_mod1  <- glmer(form, family="poisson",data=grouseticks))
}
\keyword{datasets}

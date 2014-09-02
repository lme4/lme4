\name{Arabidopsis}
\alias{Arabidopsis}
\docType{data}
\title{
  Arabidopsis clipping/fertilization data
}
\description{
  Data on genetic variation in responses to fertilization
  and simulated herbivory in \emph{Arabidopsis}
}
\usage{data("Arabidopsis")}
\format{
  A data frame with 625 observations on the following 8 variables.
  \describe{
    \item{\code{reg}}{region: a factor with 3 levels \code{NL}
      (Netherlands), \code{SP} (Spain), \code{SW} (Sweden)}
    \item{\code{popu}}{population: a factor with the form \code{n.R}
      representing a population in region \code{R}}
    \item{\code{gen}}{genotype: a factor with 24 (numeric-valued) levels}
    \item{\code{rack}}{a nuisance factor with 2 levels, one for
      each of two greenhouse racks}
    \item{\code{nutrient}}{fertilization treatment/nutrient level (1,
      minimal nutrients or 8, added nutrients)}
    \item{\code{amd}}{simulated herbivory or "clipping" (apical meristem
      damage): \code{unclipped} (baseline) or \code{clipped}}
    \item{\code{status}}{a nuisance factor for germination method
      (\code{Normal}, \code{Petri.Plate}, or \code{Transplant})}
    \item{\code{total.fruits}}{total fruit set per plant (integer)}
  }
}
\source{
  From Josh Banta
}
\references{
  Joshua A. Banta, Martin H. H Stevens, and Massimo Pigliucci (2010)
  A comprehensive test of the 'limiting resources' framework applied to
  plant tolerance to apical meristem damage.
  \emph{Oikos} \bold{119}(2), 359--369;
  \url{http://dx.doi.org/10.1111/j.1600-0706.2009.17726.x}
}
\examples{
data(Arabidopsis)
summary(Arabidopsis[,"total.fruits"])
table(gsub("[0-9].","",levels(Arabidopsis[,"popu"])))
library(lattice)
stripplot(log(total.fruits+1) ~ amd|nutrient, data = Arabidopsis,
          groups = gen,
          strip=strip.custom(strip.names=c(TRUE,TRUE)),
          type=c('p','a'), ## points and panel-average value --
          ## see ?panel.xyplot
          scales=list(x=list(rot=90)),
          main="Panel: nutrient, Color: genotype")
}
\keyword{datasets}

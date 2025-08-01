\name{gopherdat2}
\alias{gopherdat2}
\docType{data}
\title{
  Gopher tortoises shell remains
}
\description{
  Data of fresh Gopher tortoise shell remains \insertCite{ozgul2009upper}{lme4}.
}
\usage{data("gopherdat2")}
\format{
  A data frame with 30 observations on the following 7 variables.
  \describe{
    \item{\code{Site}}{a factor representing the sampling site.  There are 
      \code{10} levels: 
        \itemize{
          \item \code{BS} (Big Shoals State Park).
          \item \code{CB} (Camp Blanding Wildlife Management Area).
          \item \code{Cent} (privately-owned property in central Florida).
          \item \code{CF} (Cecil Field/Branan Field Wildlife and Environmental 
            Area).
          \item \code{FC} (Fort Cooper State Park).
          \item \code{FE} (Flying Eagle Wildlife Management Area).
          \item \code{GH} (Gold Head Branch State Park).
          \item \code{Old} (Perry Oldenburg Wildlife and Environmental Area).
          \item \code{Ord} (Ordway-Swisher Biological Station). 
          \item \code{TE} (Tenoroc Fish Management Area).
        }
      }
    \item{\code{year}}{a numeric vector of the sampling year.}  
    \item{\code{shells}}{a numeric vector of the number of shells found.}
    \item{\code{Area}}{a numeric vector representing site area 
      (units unknown).}
    \item{\code{density}}{a numeric vector representing population density.}
    \item{\code{prev}}{a numeric vector representing the seroprevalence 
      (frequency of antibodies to disease) of \emph{M. agassizi}. }
  }
}
\details{
  Details of the study is described in \insertCite{ozgul2009upper}{lme4} as
  follows:
  "The fieldwork was conducted between 2003 and 2006, during late spring/summer 
  (May–September). Systematic surveys were conducted to locate tortoise burrows 
  and shell remains from deceased tortoises at study sites and consisted of a 
  line of four to eight observers spaced 10 m apart walking parallel transects 
  across the study area."
}
\source{
  \insertCite{ozgul2009upper}{lme4}
}
\references{
  \insertRef{ozgul2009upper}{lme4}
}
\examples{
## Simple model gives a singular fit:
gopher_glmer <- glmer(shells ~ factor(year) + prev + offset(log(Area))
                   + (1|Site), data = gopherdat2, family = "poisson")
## The site-level variance for this model is indeed zero:
VarCorr(gopher_glmer)
## So a Poisson GLM gives the same answer here:
gopher_glm <- glm(shells ~ factor(year) + prev + offset(log(Area)),
                  data = gopherdat2, family = "poisson")
all.equal(fixef(gopher_glmer), coef(gopher_glm))
}
\keyword{datasets}

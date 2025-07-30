\name{gopherdat2}
\alias{gopherdat2}
\docType{data}
\title{
  Gopher tortoises shell remains
}
\description{
  Data of fresh Gopher tortoise shell remains from Ozgul et al.
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
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  From Arpat Ozgul, Madan K. Oli, Benjamin M. Bolker, and Carolina 
  Perez-Heydrich.
}
\references{
  Arpat Ozgul, Madan K. Oli, Benjamin M. Bolker, and Carolina Perez-Heydrich 
  (2009) Upper respiratory tract disease, force of infection, and effects on 
  survival of gopher tortoises. \emph{Ecological Applications} \strong{19}(3), 
  786--798; \doi{10.1890/08-0219.1}
}
\examples{
## Doing this directly will obtain a singularFit issue:
gopher_glmer <- glmer(shells ~ factor(year) + prev + offset(log(Area))
                   + (1|Site), data = gopherdat2, family = "poisson")

## Trying to ignore the singularFit portion:
gopher_glmer <- glmer(
  shells ~ prev + offset(log(Area)) + factor(year) + (1 | Site), 
  data = gopherdat2,  
  family = poisson,                                               
  control = glmerControl(                                        
    optimizer = "bobyqa",                                         
    check.conv.singular = .makeCC(action = "ignore", tol = 1e-4) 
  )
)
## Inspecting for the value of the site-level variance                                          
VarCorr(gopher_glmer)
## This suggests we could have obtained similar results from a Poisson GLM
## Trying to create a simple Poisson GLM:
gopher_glm <- glm(shells ~ factor(year) + prev + offset(log(Area)),
                 data = gopherdat2, family = "poisson")
## Comparing the fixed coefficients, which appear to be identical:
mod_list = list("glm" = gopher_glm, "glmer" = gopher_glmer)
purrr::map_dfr(mod_list, ~broom.mixed::tidy(., effects = "fixed"), 
               .id ="model") |> dplyr::arrange(term)
## As the table suggests, the estimates are similar.
}
\keyword{datasets}

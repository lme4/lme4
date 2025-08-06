\docType{data}
\name{cake}
\alias{cake}
\title{Breakage Angle of Chocolate Cakes}
\format{A data frame with 270 observations on the following 5 variables.
  \describe{
    \item{\code{replicate}}{a factor with levels \code{1} to \code{15}}
    \item{\code{recipe}}{a factor with levels \code{A}, \code{B} and \code{C}}
    \item{\code{temperature}}{an ordered factor with levels \code{175}
      < \code{185} < \code{195} < \code{205} < \code{215} < \code{225}}
    \item{\code{angle}}{a numeric vector giving the angle at which the
      cake broke.}
    \item{\code{temp}}{numeric value of the baking temperature (degrees F).}
  }}
\source{
  Original data were presented in Cook (1938), and reported
  in Cochran and Cox (1957, p. 300).  Also cited in Lee,
  Nelder and Pawitan (2006).
}
\description{
  Data on the breakage angle of chocolate cakes made with
  three different recipes and baked at six different
  temperatures.  This is a split-plot design with the
  recipes being whole-units and the different temperatures
  being applied to sub-units (within replicates). The
  experimental notes suggest that the replicate numbering
  represents temporal ordering.
}
\details{
  The \code{replicate} factor is nested within the
  \code{recipe} factor, and \code{temperature} is nested
  within \code{replicate}.
}
\references{
  \itemize{
    \item \insertRef{cook1938chocolate}{lme4}
    \item \insertRef{cochran1957experimental}{lme4}
    \item \insertRef{lee2006generalized}{lme4}
  }
}
\examples{
str(cake)
## 'temp' is continuous, 'temperature' an ordered factor with 6 levels

(fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake, REML= FALSE))
(fm2 <- lmer(angle ~ recipe + temperature + (1|recipe:replicate), cake, REML= FALSE))
(fm3 <- lmer(angle ~ recipe + temp        + (1|recipe:replicate), cake, REML= FALSE))

## and now "choose" :
anova(fm3, fm2, fm1)
}
\keyword{datasets}


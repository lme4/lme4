\name{culcitalvolume}
\alias{culcitalvolume}
\docType{data}
\title{
 Coral-eating seastar Culcita novaeguineae data (volume loss version)
}
\description{
  This is the same experiment outlined in \code{\link[lme4]{culcitalogreg}} 
  except this data set only contains values where predation occurred, and 
  the original volume as well as volume lost was recorded.
}
\usage{data("culcitalvolume")}
\format{
  A data frame with 50 observations on the following 7 variables.
  \describe{
    \item{\code{ttt}}{represents the combinations of different symbionts,
      which the treatment condition is explicitly written in plain language.}
    \item{\code{volume}}{describes the size of the coral, measured in 
      cm\eqn{^3}.}
    \item{\code{predvolume}}{the amount of volume loss from the coral, measured
      in cm\eqn{^3}.}
    \item{\code{crab}}{describes whether the crab was present in the 
      experiment. \code{n}: not present, \code{y}: present.}
    \item{\code{shrimp}}{describes whether the shrimp was present in the 
      experiment. \code{n}: not present, \code{y}: present.}
    \item{\code{block}}{a numeric variable indicating the experimental block. 
      There are \code{10} blocks in total, each corresponding to a large, 
      octagonal, flow-through seawater tank approximately 0.5m deep and 2m in 
      diameter.}
    \item{\code{ttt2}}{a relabelled version of \code{ttt}.
      \itemize{
          \item \code{1}: no exosymbionts,
          \item \code{2}: pair of \emph{Alpheus lottini} only 
            (Alpheus; shrimp),
          \item \code{3}: pair of \emph{Trapezia serenei} only 
            (Trapezia; crab),
          \item \code{4}: pair of \emph{Alpheus lottini} and pair of 
            \emph{Trapezia serenei} (‘Alpheus and Trapezia’).
        }
      }
  }
}
\seealso{
  The version which contains whether or not predation occurred, 
  \code{\link[lme4]{culcitalogreg}}.
}
\source{
  McKeon, C.S., Stier, A.C., McIlroy, S.E. and Bolker, B.M.
}
\references{
  McKeon, C.S., Stier, A.C., McIlroy, S.E. and Bolker, B.M. Multiple defender 
  effects: synergistic coral defense by mutualist crustaceans. Oecologia 169, 
  1095–1103 (2012); \doi{10.1007/s00442-012-2275-2}.
}
\examples{
## Modifying to create a new response variable
vdata <- transform(culcitalvolume, 
                   propeaten = predvolume/volume,
                   tvol = log(predvolume))
## One-way analysis
(cvm1 <- lmer(tvol ~ ttt2 + (1|block), data = vdata))
(cvm2 <- lmer(propeaten ~ ttt2 + (1|block), data = vdata))
## Two-way analysis 
(cvm3 <- lmer(tvol ~ crab*shrimp + (1|block), data = vdata))
(cvm4 <- lmer(propeaten ~ crab*shrimp + (1|block), data = vdata))
}
\keyword{datasets}

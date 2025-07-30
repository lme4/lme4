\name{salamander}
\alias{salamander}
\docType{data}
\title{
  Mountain dusky salamander mating
}
\description{
  S. Arnold and P. Verrell conducted an experiment at the University of 
  Chicago to study breeding behaviours of mountain dusky salamanders. This
  mating data was used as an example from McCullagh and Nelder's (1989) book 
  on GLMs.
}
\usage{data("salamander")}
\format{
  A data frame with 360 observations on the following 8 variables.
  \describe{
    \item{\code{Season}}{represents the season \code{Fall} \code{Summer} of 
      1986.}
    \item{\code{Experiment}}{experiment number \code{1, 2, 3}.}
    \item{\code{TypeM}}{representing the two types of male salamanders being 
      studied. \code{R} stands for rough butt, and \code{W} stands for white 
      side.}
    \item{\code{TypeF}}{similar to the above except for female salamanders.}
    \item{\code{Cross}}{represents the cross between a female and male type. 
      \itemize{
        \item{\code{RR} = rough butt female crossed with rough butt male},
        \item{\code{RW} = rough butt female crossed with white side male},
        \item{\code{WR} = white side female crossed with rough butt male},
        \item{\code{WW} = white side female crossed with white side male}.
      }
    }
    \item{\code{Male}}{identication number of the male salamander.}
    \item{\code{Female}}{identification number of the male salamander.}
    \item{\code{Mate}}{represents whether mating has occurred. \code{1} = yes,
      \code{0} = no.}
  }
}
\details{
  In this example, every variable is either binary or a factor.  However, 
  most of the variables in this data set are treated like a numerical variable
  (i.e., \code{Experiment, Male, Female, Mate}), so we suggest turning them
  into factor variables before use. As outlined in McCullagh and Nelder (1989): 
  This experiment, conducted by S. Arnold and P. Verrell at the University of 
  Chicago, investigated whether geographically isolated populations of mountain 
  dusky salamanders would interbreed. The goal was to see the difference in 
  mating frquencies between the \code{RW} crosses compared to the \code{WR}
  crosses. Because each salamander was involved in multiple mating trials with 
  different partners, the observations are not independent. Therefore, we use a 
  model (see the \emph{examples} section) that conditions on the specific male 
  and female salamanders in the experiment.
}
\source{
  McCullagh P. and Nelder, J. A.
}
\references{
  McCullagh P. and Nelder, J. A. 1989. \emph{Generalized Linear Models}, 
  Section 14.5, Chapman and Hall/CRC.
}
\examples{
## Making sure Male, Female, and CRoss are treated as factors
salamander$Male <- factor(salamander$Male)
salamander$Female <- factor(salamander$Female)
salamander$Cross <- factor(salamander$Cross)
## Fitting the model described in 14.5.3 from McCullagh and Nelder
sal_mod <- glmer(Mate ~ (1|Female) + (1 | Male) + Cross, data = salamander,
                 family = binomial(link = "logit"))
}
\keyword{datasets}

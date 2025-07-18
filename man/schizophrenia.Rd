\name{schizophrenia}
\alias{schizophrenia}
\docType{data}
\title{
  National Institute of Mental Health Schizophrenia Collaborative Study
}
\description{
  In the NIMH Schizophrenia Collaborative Study, patients were randomized to 
  receive one of four medications, either placebo or one of three different 
  anti-psychotic drugs. The protocol indicated subjects were to then be 
  evaluated at weeks 0, 1, 3, 6 to assess severity of illness; additionally 
  some measurements were made at weeks 2, 4, and 5.
}
\usage{data(schizophrenia)}
\format{
  A data frame with 1603 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{a numeric vector indicating the unique patient 
      identifier.}
    \item{\code{imps79}}{a numeric vector from item 79 on the Inpatient 
      Multidimensional Psychiatric Scale indicating severity of illness.  The 
      scores were applied to be interpreted as follows:
      \itemize{
        \item{\code{1} = normal, not ill at all},
        \item{\code{2} = borderline mentally ill},
        \item{\code{3} = mildly ill},
        \item{\code{4} = moderately ill},
        \item{\code{5} = markedly ill},
        \item{\code{6} = severely ill},
        \item{\code{7} = among the most extremely ill}.
      }
    }
    \item{\code{imps79b}}{a binary version of imps79.}
    \item{\code{imps79o}}{an ordinally scaled version of imps79.}
    \item{\code{int}}{a numeric vector of ones; used in stand-alone package to 
      indicate intercept.}
    \item{\code{TxDrug}}{a numeric vector indicating treatment with drug 
      \code{1} or placebo \code{0}.}
    \item{\code{Week}}{a numeric vector indicating time, in weeks.}
    \item{\code{SqrtWeek}}{the square root of the Week variable.}
    \item{\code{TxSWeek}}{a variable representing the \code{TxDrug * Week} 
      interaction.}
  }
}

\source{
  This data and most of this documentation was taken from the \code{mixor} 
  package, which is no longer from CRAN. The references below mostly reference 
  the deprecated \code{mixor} program.
}
\references{
  Hedeker D. and Gibbons R.D. (1996) A computer program for mixed-effects 
  ordinal regression analysis. \emph{Computer Methods and Programs in 
  Biomedicine} 49, 157-176.
  Hedeker D and Gibbons R.D. (2006) \emph{Longitudinal Data Analysis}, Wiley, 
  Hoboken, New Jesery.
}
\examples{
schmod <- glmer(imps79 ~ TxDrug * Week + (1 | id),
                data = schizophrenia, family = Gamma(link = "log"))
}
\keyword{datasets}

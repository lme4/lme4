\name{culcitalogreg}
\alias{culcitalogreg}
\docType{data}
\title{
 Coral-Eating Seastar Culcita novaeguineae Dataset
}
\description{
  TODO LOTS...
}
\usage{data("culcitalogreg")}
\format{
  A data frame with 80 observations on the following 9 variables.
  \describe{
    \item{\code{block}}{a numeric vector}
    \item{\code{ttt}}{a numeric vector}
    \item{\code{predation}}{a numeric vector}
    \item{\code{ttt.1}}{a character vector}
    \item{\code{crab}}{a character vector}
    \item{\code{shrimp}}{a character vector}
    \item{\code{ttt2}}{a character vector}
    \item{\code{crab2}}{a character vector}
    \item{\code{shrimp2}}{a character vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  McKeon, C.S., Stier, A.C., McIlroy, S.E. et al. Multiple defender effects: 
  synergistic coral defense by mutualist crustaceans. Oecologia 169, 1095–1103 
  (2012); \doi{10.1007/s00442-012-2275-2}.
}
\examples{
cul_mod <- glmer(predation ~ ttt + (1|block), data=culcitalogreg,
                     family = binomial(link = "logit"))
}
\keyword{datasets}

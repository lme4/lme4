\name{Dyestuff}
\alias{Dyestuff}
\alias{Dyestuff2}
\docType{data}
\title{Yield of dyestuff by batch}
\description{
  The \code{Dyestuff} data frame provides the yield of dyestuff
  (Naphthalene Black 12B) from 5 different preparations from each of 6
  different batchs of an intermediate product (H-acid).  The
  \code{Dyestuff2} data were generated data in the same structure but
  with a large residual variance relative to the
  batch variance.
}
\usage{data(Dyestuff)}
\format{
  Data frames, each with 30 observations on the following 2 variables.
  \describe{
    \item{\code{Batch}}{a factor indicating the batch of the
      intermediate product from which the preparation was created.}
    \item{\code{Yield}}{the yield of dyestuff from the preparation
      (grams of standard color).}
  }
}
\details{
 The \code{Dyestuff} data are described in Davies and Goldsmith (1972)
 as coming from
 \dQuote{an investigation to find out how much the
   variation from batch to batch in the quality of an intermediate
   product (H-acid) contributes to the variation in the yield of the
   dyestuff (Naphthalene Black 12B) made from it.  In the experiment six
   samples of the intermediate, representing different batches of works
   manufacture, were obtained, and five preparations of the dyestuff
   were made in the laboratory from each sample.  The equivalent yield
   of each preparation as grams of standard colour was determined by
   dye-trial.}

 The \code{Dyestuff2} data are described in Box and Tiao (1973) as
 illustrating
 \dQuote{
   the case where between-batches mean square is less than the
   within-batches mean square.  These data had to be constructed for
   although examples of this sort undoubtably occur in practice, they
   seem to be rarely published.}
}
\source{
 O.L. Davies and P.L. Goldsmith (eds), \emph{Statistical Methods in
   Research and Production, 4th ed.}, Oliver and Boyd, (1972), section
   6.4

 G.E.P. Box and G.C. Tiao, \emph{Bayesian Inference in Statistical
   Analysis}, Addison-Wesley, (1973), section 5.1.2
}
%\references{}
\examples{
str(Dyestuff)
dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff,
        ylab = "Batch", jitter.y = TRUE, aspect = 0.3,
        type = c("p", "a"))
dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff2,
        ylab = "Batch", jitter.y = TRUE, aspect = 0.3,
        type = c("p", "a"))
(fm1 <- lmer(Yield ~ 1|Batch, Dyestuff))
(fm2 <- lmer(Yield ~ 1|Batch, Dyestuff2))
}
\keyword{datasets}

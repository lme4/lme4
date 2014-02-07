\docType{data}
\name{Dyestuff}
\alias{Dyestuff}
\alias{Dyestuff2}
\title{Yield of dyestuff by batch}
\format{Data frames, each with 30 observations on the following 2 variables.
  \describe{
    \item{\code{Batch}}{a factor indicating the batch of the
      intermediate product from which the preparation was created.}
    \item{\code{Yield}}{the yield of dyestuff from the preparation
      (grams of standard color).}
  }}
\source{
  O.L. Davies and P.L. Goldsmith (eds), \emph{Statistical
  Methods in Research and Production, 4th ed.}, Oliver and
  Boyd, (1972), section 6.4

  G.E.P. Box and G.C. Tiao, \emph{Bayesian Inference in
  Statistical Analysis}, Addison-Wesley, (1973), section
  5.1.2
}
\description{
  The \code{Dyestuff} data frame provides the yield of
  dyestuff (Naphthalene Black 12B) from 5 different
  preparations from each of 6 different batchs of an
  intermediate product (H-acid).  The \code{Dyestuff2} data
  were generated data in the same structure but with a
  large residual variance relative to the batch variance.
}
\details{
  The \code{Dyestuff} data are described in Davies and
  Goldsmith (1972) as coming from \dQuote{an investigation
  to find out how much the variation from batch to batch in
  the quality of an intermediate product (H-acid)
  contributes to the variation in the yield of the dyestuff
  (Naphthalene Black 12B) made from it.  In the experiment
  six samples of the intermediate, representing different
  batches of works manufacture, were obtained, and five
  preparations of the dyestuff were made in the laboratory
  from each sample. The equivalent yield of each
  preparation as grams of standard colour was determined by
  dye-trial.}

  The \code{Dyestuff2} data are described in Box and Tiao
  (1973) as illustrating \dQuote{ the case where
  between-batches mean square is less than the
  within-batches mean square.  These data had to be
  constructed for although examples of this sort
  undoubtably occur in practice, they seem to be rarely
  published.}
}
\examples{
\dontshow{ # useful for the lme4-authors --- development, debugging, etc:
 commandArgs()[-1]
 if(FALSE) ## R environment variables:
 local({ ne <- names(e <- Sys.getenv())
         list(R    = e[grep("^R", ne)],
              "_R" = e[grep("^_R",ne)]) })
 Sys.getenv("R_ENVIRON")
 Sys.getenv("R_PROFILE")
 cat("R_LIBS:\\n"); (RL <- strsplit(Sys.getenv("R_LIBS"), ":")[[1]])
 nRL <- normalizePath(RL)
 cat("and extra(:= not in R_LIBS) .libPaths():\\n")
 .libPaths()[is.na(match(.libPaths(), nRL))]

 structure(Sys.info()[c(4,5,1:3)], class="simple.list") #-> 'nodename' ..
 sessionInfo()
 searchpaths()
 pkgI <- function(pkgname) {
   pd <- tryCatch(packageDescription(pkgname),
                  error=function(e)e, warning=function(w)w)
   if(inherits(pd, "error") || inherits(pd, "warning"))
     cat(sprintf("packageDescription(\\"\%s\\") \%s: \%s\\n",
                 pkgname, class(pd)[2], pd$message))
   else
     cat(sprintf("\%s -- built: \%s\\n\%*s -- dir  : \%s\\n",
                 pkgname, pd$Built, nchar(pkgname), "",
                 dirname(dirname(attr(pd, "file")))))
 }
 pkgI("Matrix")
 pkgI("Rcpp")
 ## 2012-03-12{MM}: fails with --as-cran
 pkgI("RcppEigen")
 pkgI("minqa")
 pkgI("lme4")
}
require(lattice)
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


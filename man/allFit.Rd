\name{allFit}
\alias{allFit}
\title{Refit a fitted model with all available optimizers}
\usage{
allFit(m, meth.tab = NULL, data=NULL,
        verbose = TRUE,
        show.meth.tab = FALSE,
        maxfun = 1e5,
	parallel = c("no", "multicore", "snow"),
	ncpus = getOption("allFit.ncpus", 1L), cl = NULL)
}
\arguments{
  \item{m}{a fitted model}
  \item{meth.tab}{a matrix (or data.frame) with columns
    \describe{
      \item{method}{the name of a specific optimization method to pass
	to the optimizer (leave blank for built-in optimizers)}
      \item{optimizer}{the \code{optimizer} function to use}
    }
  }
  \item{data}{data to be included with result (for later debugging etc.)}
  \item{verbose}{logical: report progress in detail?}
  \item{show.meth.tab}{logical: return table of methods?}
  \item{maxfun}{passed as part of \code{optCtrl} (see \code{\link{lmeControl}})}
  \item{parallel}{The type of parallel operation to be used (if any).
    If missing, the
    default is taken from the option \code{"boot.parallel"} (and if that
    is not set, \code{"no"}).}
  \item{ncpus}{integer: number of processes to be used in parallel operation:
    typically one would choose this to be the number of available CPUs.}
  \item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
    \code{parallel = "snow"}.  If not supplied, a cluster on the
    local machine is created for the duration of the \code{boot} call.}
}
\value{
  a list of fitted \code{merMod} objects (unless \code{show.meth.tab} is
  specified, in which case a data frame of methods is returned)
}
\description{
   Attempt to re-fit a [g]lmer model with a range of optimizers.
The default is to use all known optimizers for R that satisfy the
requirements (do not require explicit gradients, allow
box constraints), in three categories; (i) built-in
(minqa::bobyqa, lme4::Nelder_Mead, nlminbwrap), (ii) wrapped via optimx
(most of optimx's optimizers that allow box constraints require
an explicit gradient function to be specified; the two provided
here are really base R functions that can be accessed via optimx,
(iii) wrapped via nloptr; (iv)
}
\note{
  If you are using \code{parallel="snow"}, you will need to run
  \code{clusterEvalQ(cl,library("lme4"))} before calling
  \code{allFit} to make sure that the
  \code{lme4} package is loaded on all of the workers
}
\seealso{\code{slice},\code{slice2D} from the \code{bbmle} package}
\examples{
library(lme4)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
## show available methods
allFit(show.meth.tab=TRUE) 
gm_all <- allFit(gm1)
ss <- summary(gm_all)
ss$fixef               ## extract fixed effects
ss$llik                ## log-likelihoods
ss$sdcor               ## SDs and correlations
ss$theta               ## Cholesky factors
ss$which.OK            ## which fits worked
}
\keyword{models}

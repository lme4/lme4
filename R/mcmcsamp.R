##' @name pvalues
##' @aliases mcmcsamp
##' @title Getting p-values for fitted models
##' 
##' @description One of the most frequently asked questions about \code{lme4}
##' is "how do I calculate p-values for estimated parameters?"
##' Previous versions of \code{lme4} provided the \code{mcmcsamp}
##' function, which efficiently generated a Markov chain Monte Carlo sample
##' from the posterior distribution of the parameters, assuming
##' flat (scaled likelihood) priors. Due to difficulty in
##' constructing a version of \code{mcmcsamp} that was reliable
##' even in cases where the estimated random effect variances were near
##' zero (e.g. \url{https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q4/003115.html}), \code{mcmcsamp} has been withdrawn (or more precisely,
##' not updated to work with \code{lme4} versions >=1.0.0).
##'
##' Many users, including users of the \code{aovlmer.fnc} function
##' from the \code{languageR} package which relies on \code{mcmcsamp},
##' will be deeply disappointed by this lacuna. Users who need p-values have
##' a variety of options:
##' \itemize{
##' \item likelihood ratio tests via \code{anova} (MC,+)
##' \item profile confidence intervals via \code{\link{profile.merMod}} and \code{\link{confint.merMod}} (CI,+)
##' \item parametric bootstrap confidence intervals and model comparisons via \code{\link{bootMer}} (or \code{PBmodcomp} in the \code{pbkrtest} package) (MC/CI,*,+)
##' \item for random effects, simulation tests via the \code{RLRsim} package (MC,*)
##' \item for fixed effects, F tests via Kenward-Roger approximation using \code{KRmodcomp} from the \code{pbkrtest} package (MC)
##' \item \code{car::Anova} and \code{lmerTest::anova} provide wrappers for \code{pbkrtest}. \code{lmerTest::anova} also provides t tests via the Satterthwaite approximation (P,*)
##' }
##' In the list above, the methods marked \code{MC} provide explicit model comparisons; \code{CI} denotes confidence intervals; and \code{P} denotes parameter-level or sequential tests of all effects in a model. The starred (*) suggestions provide finite-size corrections (important when the number of groups is <50); those marked (+) support GLMMs as well as LMMs.
##' 
##' When all else fails, don't forget to keep p-values in perspective: \url{http://www.phdcomics.com/comics/archive.php?comicid=905}
##' 
if(FALSE) ## C++ code in ../src/mcmcsamp.cpp -- is also  #ifdef 0
# @S3method mcmcsamp merMod
mcmcsamp.merMod <- function(object, n=1L, verbose=FALSE, saveb=FALSE, ...) {
    n <- max(1L, as.integer(n)[1])
    dd <- getME(object, "devcomp")$dims
    ranef <- matrix(numeric(0), nrow = dd[["q"]], ncol = 0)
    if (saveb) ranef <- matrix(, nrow = dd[["q"]], ncol = n)
    sigma <- matrix(unname(sigma(object)), nrow = 1,
                    ncol = (if (dd[["useSc"]]) n else 0))
    ff <- fixef(object)
    fixef <- matrix(ff, nrow=dd[["p"]], ncol=n)
    rownames(fixef) <- names(ff)
    ## FIXME create a copy of the resp and pred modules
    ans <- new("merMCMC",
               Gp = object@Gp,
 #              ST = matrix(.Call(mer_ST_getPars, object), dd[["np"]], n),
               call = object@call,
               dims = object@dims,
               deviance = rep.int(unname(object@deviance[["ML"]]), n),
               fixef = fixef,
               nc = sapply(object@ST, nrow),
               ranef = ranef,
               sigma = sigma)
    .Call(mer_MCMCsamp, ans, object)
}

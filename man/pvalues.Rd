\name{pvalues}
\alias{mcmcsamp}
\alias{pvalues}
\title{Getting p-values for fitted models}
\description{
  One of the most frequently asked questions about
  \code{lme4} is "how do I calculate p-values for estimated
  parameters?" Previous versions of \code{lme4} provided
  the \code{mcmcsamp} function, which efficiently generated
  a Markov chain Monte Carlo sample from the posterior
  distribution of the parameters, assuming flat (scaled
  likelihood) priors. Due to difficulty in constructing a
  version of \code{mcmcsamp} that was reliable even in
  cases where the estimated random effect variances were
  near zero (e.g.
  \url{https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q4/003115.html}),
  \code{mcmcsamp} has been withdrawn (or more precisely,
  not updated to work with \code{lme4} versions >=1.0.0).

  Many users, including users of the \code{aovlmer.fnc} function from
  the \code{languageR} package which relies on \code{mcmcsamp}, will be
  deeply disappointed by this lacuna. Users who need p-values have a
  variety of options. In the list below, the methods marked \code{MC}
  provide explicit model comparisons; \code{CI} denotes confidence
  intervals; and \code{P} denotes parameter-level or sequential tests of
  all effects in a model. The starred (*) suggestions provide
  finite-size corrections (important when the number of groups is <50);
  those marked (+) support GLMMs as well as LMMs.

  \itemize{
    \item likelihood ratio tests via \code{anova} (MC,+)
    \item{profile confidence intervals via \code{\link{profile.merMod}} and
      \code{\link{confint.merMod}} (CI,+)
    }
    \item{parametric bootstrap confidence intervals and model comparisons via
      \code{\link{bootMer}} (or \code{PBmodcomp} in the
      \code{pbkrtest} package) (MC/CI,*,+)
    }
    \item{for random
      effects, simulation tests via the \code{RLRsim} package
      (MC,*)
    }
    \item{for fixed effects, F tests via Kenward-Roger
      approximation using \code{KRmodcomp} from the
      \code{pbkrtest} package (MC)
    }
    \item{\code{car::Anova} and
      \code{lmerTest::anova} provide wrappers for
      \code{pbkrtest}: \code{lmerTest::anova} also provides t tests via the
      Satterthwaite approximation (P,*)
    }
    \item{\code{afex::mixed} is another wrapper for
      \code{pbkrtest} and \code{anova} providing
      "Type 3" tests of all effects (P,*,+)
    }
  }

  \code{arm::sim}, or \code{\link{bootMer}}, can be used
  to compute confidence intervals on predictions.

  When all else fails, don't forget to keep p-values in
  perspective:
  \url{http://www.phdcomics.com/comics/archive.php?comicid=905}
}


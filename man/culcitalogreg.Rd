\name{culcitalogreg}
\alias{culcitalogreg}
\docType{data}
\title{
 Coral-eating seastar Culcita novaeguineae data (binary predation version)
}
\description{
  This data shows the results from a block design experiment to record the 
  interactions of ‘guard’ crustaceans, shrimp (\emph{Alpheus lottini}) and crab 
  (\emph{Trapezia serenei}), that defend their coral host 
  (\emph{Pocillopora cf. meandrina}) from seastar predators 
  (\emph{Culcita novaeguineae}). The coral and their seastar predator were 
  placed in large, octagonal, flow-through seawater tanks measuring 
  approximately 0.5m deep and 2m in diameter. There were 10 tanks in total, 
  each serving as an experimental block. Four treatment groups were tested: 
  (1) no exosymbionts, (2) shrimp only, (3) crab only, and (4) both crab and 
  shrimp present. Each treatment was replicated twice within each block.
}
\usage{data("culcitalogreg")}
\format{
  A data frame with 80 observations on the following 9 variables.
  \describe{
    \item{\code{block}}{a numeric variable indicating the experimental block. 
    There are \code{10} blocks in total, each corresponding to a large, 
    octagonal, flow-through seawater tank approximately 0.5m deep and 2m in 
    diameter.}
    \item{\code{ttt}}{a number which represents the combinations of different 
      symbionts.
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
    \item{\code{predation}}{(binary) indicates whether the seastar predator 
      consumed the coral. \code{0} = did not eat; \code{1} = ate.}
    \item{\code{ttt.1}}{a relabelled version of \code{ttt} that explicitly 
      indicates the corresponding treatment condition in plain language.}
    \item{\code{crab}}{describes whether the crab was present in the 
      experiment. \code{C-}: not present, \code{C+}: present.}
    \item{\code{shrimp}}{describes whether the shrimp was present in the
      experiment. \code{S-}: not present, \code{S+}: present.}
    \item{\code{ttt2}}{a relabelled version of \code{ttt} that uses letters 
      (a, b, c, d) instead of numeric values (1, 2, 3, 4), respectively, to 
      represent the treatment groups.}
    \item{\code{crab2}}{a relabelling of \code{crab} where \code{n} corresponds 
      to \code{C-} (no crab), and \code{y} corresponds to \code{C+} (crab 
      present).}
    \item{\code{shrimp2}}{a relabelling of \code{shrimp} where \code{n} 
      corresponds to \code{S-} (no shrimp), and \code{y} corresponds to 
      \code{S+} (shrimp present).}
  }
}
\details{
  Full details of the experiment, as outlined by McKeon et al., is quoted below:
  "Twenty replicates in total were conducted for each exosymbiont treatment in 
  a temporally blocked design with two replicates in each of ten temporal 
  blocks. Trials were conducted in a large octagonal flow-through seawater tank
  approximately 0.5 m deep and 2 m across. The tank was divided into eight equal 
  sections using plastic screening. Each section was provisioned with a seastar 
  refugium constructed from concrete blocks. We placed Pocillopora colonies 
  into the tank from the field during mid- to late afternoon. To minimize 
  variation driven by search time, a single Culcita was placed directly on top 
  of the coral colony. The following morning, we measured the coral size 
  (length, width, height) and the feeding scars left by the Culcita (length, 
  width, depth)."
}
\seealso{
  The version which only contains predation counts but includes the amount of 
  volume of the coral and the amount of volume lost, 
  \code{\link[lme4]{culcitalvolume}}.
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
cul_mod <- glmer(predation ~ ttt2 + (1|block), data=culcitalogreg,
                     family = binomial(link = "logit"))
}
\keyword{datasets}

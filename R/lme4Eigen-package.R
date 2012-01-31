### roxygen2 documentation for data sets in the package

##' Breakage angle of chocolate cakes
##' 
##' Data on the breakage angle of chocolate cakes made with three different
##' recipes and baked at six different temperatures.  This is a split-plot
##' design with the recipes being whole-units and the different temperatures
##' being applied to sub-units (within replicates). The experimental notes
##' suggest that the replicate numbering represents temporal ordering.
##' 
##' The \code{replicate} factor is nested within the \code{recipe} factor, and
##' \code{temperature} is nested within \code{replicate}.
##' 
##' @name cake
##' @docType data
##' @format A data frame with 270 observations on the following 5 variables.
##'   \describe{
##'     \item{\code{replicate}}{a factor with levels \code{1} to \code{15}}
##'     \item{\code{recipe}}{a factor with levels \code{A}, \code{B} and \code{C}}
##'     \item{\code{temperature}}{an ordered factor with levels \code{175}
##'       < \code{185} < \code{195} < \code{205} < \code{215} < \code{225}}
##'     \item{\code{angle}}{a numeric vector giving the angle at which the
##'       cake broke.}
##'     \item{\code{temp}}{numeric value of the baking temperature (degrees F).}
##'   }
##' @references Cook, F. E. (1938) \emph{Chocolate cake, I. Optimum baking
##'   temperature}. Master's Thesis, Iowa State College.
##' 
##'   Cochran, W. G., and Cox, G. M. (1957) \emph{Experimental designs}, 2nd Ed.
##'   New York, John Wiley \& Sons.
##' 
##'   Lee, Y., Nelder, J. A., and Pawitan, Y. (2006) \emph{Generalized linear
##'   models with random effects. Unified analysis via H-likelihood}. Boca Raton,
##'   Chapman and Hall/CRC.
##' @source Original data were presented in Cook (1938), and reported in Cochran
##'    and Cox (1957, p. 300).  Also cited in Lee, Nelder and Pawitan (2006).
##' @keywords datasets
##' @examples
##' str(cake)
##' ## 'temp' is continuous, 'temperature' an ordered factor with 6 levels
##' 
##' fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake, REML= FALSE)
##' print(fm1, corr=FALSE)
##' fm2 <- lmer(angle ~ recipe + temperature + (1|recipe:replicate), cake, REML= FALSE)
##' print(fm2, corr=FALSE)
##' fm3 <- lmer(angle ~ recipe + temp        + (1|recipe:replicate), cake, REML= FALSE)
##' fm3
##' 
##' ## and now "choose" :
##' anova(fm3, fm2, fm1)
##' 
NULL

##' Contagious bovine pleuropneumonia
##' 
##' Contagious bovine pleuropneumonia (CBPP) is a major disease of cattle in
##' Africa, caused by a mycoplasma.  This dataset describes the serological
##' incidence of CBPP in zebu cattle during a follow-up survey implemented in 15
##' commercial herds located in the Boji district of Ethiopia.  The goal of the
##' survey was to study the within-herd spread of CBPP in newly infected herds.
##' Blood samples were quarterly collected from all animals of these herds to
##' determine their CBPP status.  These data were used to compute the
##' serological incidence of CBPP (new cases occurring during a given time
##' period).  Some data are missing (lost to follow-up).
##' 
##' Serological status was determined using a competitive enzyme-linked
##' immuno-sorbent assay (cELISA).
##' 
##' @name cbpp
##' @docType data
##' @format A data frame with 56 observations on the following 4 variables.
##'   \describe{
##'     \item{\code{herd}}{A factor identifying the herd (1 to 15).}
##'     \item{\code{incidence}}{The number of new serological cases for a
##'       given herd and time period.}
##'     \item{\code{size}}{A numeric vector describing herd size at the
##'       beginning of a given time period.}
##'     \item{\code{period}}{A factor with levels \code{1} to \code{4}.}
##'   }
##' @source Lesnoff, M., Laval, G., Bonnet, P., Abdicho, S., Workalemahu, A.,
##' Kifle, D., Peyraud, A., Lancelot, R., Thiaucourt, F. (2004) Within-herd
##' spread of contagious bovine pleuropneumonia in Ethiopian highlands.
##' \emph{Preventive Veterinary Medicine} \bold{64}, 27--40.
##' @keywords datasets
##' @examples
##' 
##' ## response as a matrix
##' (m1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
##'              cbpp, binomial, nAGQ=25L))
##' dput(unname(fixef(m1)))
##' dput(unname(ranef(m1, drop=TRUE)[[1]]))
##' ## response as a vector of probabilities and usage of argument "weights"
##' m1p <- glmer(incidence / size ~ period + (1 | herd), weights = size,
##'              cbpp, binomial, nAGQ=25L)
##' dput(unname(fixef(m1p)))
##' dput(unname(ranef(m1p, drop=TRUE)[[1]]))
##' ## Confirm that these are equivalent:
##' stopifnot(all.equal(fixef(m1), fixef(m1p), tol = 1e-5),
##'           all.equal(ranef(m1), ranef(m1p), tol = 1e-5),
##'           TRUE)
##' 
##' for(m in c(m1, m1p)) {
##'     cat("-------\n\nCall: ",
##'         paste(format(getCall(m)), collapse="\n"), "\n")
##'     print(logLik(m)); cat("AIC:", AIC(m), "\n") ; cat("BIC:", BIC(m),"\n")
##' }
##' stopifnot(all.equal(logLik(m1), logLik(m1p), tol = 1e-5),
##'           all.equal(AIC(m1),    AIC(m1p),    tol = 1e-5),
##'           all.equal(BIC(m1),    BIC(m1p),    tol = 1e-5))
##' 
##' ## GLMM with individual-level variability (accounting for overdispersion)
##' cbpp$obs <- 1:nrow(cbpp)
##' (m2 <- glmer(cbind(incidence, size - incidence) ~ period +
##'     (1 | herd) +  (1|obs),
##'               family = binomial, data = cbpp))
##' 
##' 
NULL

##' Yield of dyestuff by batch
##' 
##' The \code{Dyestuff} data frame provides the yield of dyestuff (Naphthalene
##' Black 12B) from 5 different preparations from each of 6 different batchs of
##' an intermediate product (H-acid).  The \code{Dyestuff2} data were generated
##' data in the same structure but with a large residual variance relative to
##' the batch variance.
##' 
##' The \code{Dyestuff} data are described in Davies and Goldsmith (1972) as
##' coming from \dQuote{an investigation to find out how much the variation from
##' batch to batch in the quality of an intermediate product (H-acid)
##' contributes to the variation in the yield of the dyestuff (Naphthalene Black
##' 12B) made from it.  In the experiment six samples of the intermediate,
##' representing different batches of works manufacture, were obtained, and five
##' preparations of the dyestuff were made in the laboratory from each sample.
##' The equivalent yield of each preparation as grams of standard colour was
##' determined by dye-trial.}
##' 
##' The \code{Dyestuff2} data are described in Box and Tiao (1973) as
##' illustrating \dQuote{ the case where between-batches mean square is less
##' than the within-batches mean square.  These data had to be constructed for
##' although examples of this sort undoubtably occur in practice, they seem to
##' be rarely published.}
##' 
##' @name Dyestuff
##' @aliases Dyestuff Dyestuff2
##' @docType data
##' @format Data frames, each with 30 observations on the following 2 variables.
##'   \describe{
##'     \item{\code{Batch}}{a factor indicating the batch of the
##'       intermediate product from which the preparation was created.}
##'     \item{\code{Yield}}{the yield of dyestuff from the preparation
##'       (grams of standard color).}
##'   }
##' @source O.L. Davies and P.L. Goldsmith (eds), \emph{Statistical Methods in
##' Research and Production, 4th ed.}, Oliver and Boyd, (1972), section 6.4
##' 
##' G.E.P. Box and G.C. Tiao, \emph{Bayesian Inference in Statistical Analysis},
##' Addison-Wesley, (1973), section 5.1.2
##' @keywords datasets
##' @examples
##' 
##' \dontshow{ # useful for the lme4-authors --- development, debugging, etc:
##'  commandArgs()[-1]
##'  if(FALSE) ## R environment variables:
##'  local({ ne <- names(e <- Sys.getenv())
##'          list(R    = e[grep("^R", ne)],
##'               "_R" = e[grep("^_R",ne)]) })
##'  Sys.getenv("R_ENVIRON")
##'  Sys.getenv("R_PROFILE")
##'  cat("R_LIBS:\n"); (RL <- strsplit(Sys.getenv("R_LIBS"), ":")[[1]])
##'  nRL <- normalizePath(RL)
##'  cat("and extra .libPaths():\n")
##'  .libPaths()[is.na(match(.libPaths(), nRL))]
##' 
##'  sessionInfo()
##'  pkgI <- function(pkgname) {
##'    pd <- packageDescription(pkgname)
##'    cat(sprintf("%s -- built: %s\n%*s -- dir  : %s\n",
##'                pkgname, pd$Built, nchar(pkgname), "",
##'                dirname(dirname(attr(pd, "file")))))
##'  }
##'  pkgI("Matrix")
##'  pkgI("Rcpp")
##'  pkgI("RcppEigen")
##'  pkgI("minqa")
##'  pkgI("lme4Eigen")
##' }
##' str(Dyestuff)
##' dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff,
##'         ylab = "Batch", jitter.y = TRUE, aspect = 0.3,
##'         type = c("p", "a"))
##' dotplot(reorder(Batch, Yield) ~ Yield, Dyestuff2,
##'         ylab = "Batch", jitter.y = TRUE, aspect = 0.3,
##'         type = c("p", "a"))
##' (fm1 <- lmer(Yield ~ 1|Batch, Dyestuff))
##' (fm2 <- lmer(Yield ~ 1|Batch, Dyestuff2))
##' 
NULL

##' Sparse Gauss-Hermite quadrature grids
##'
##' \code{GQN} contains the non-redundant quadrature nodes and weights for
##' integration of a scalar function of a \code{d}-dimensional argument with
##' respect to the density function of the \code{d}-dimensional Gaussian
##' density function.  These are stored in a list of lists.  The outer list
##' is indexed by the dimension, \code{d}, in the range of 1 to 20.  The inner
##' list is indexed by \code{k}, the order of the quadrature.
##'
##' @note These are only the non-redundant nodes.  To regenerate the whole
##' array of nodes, all possible permutations of axes and all possible
##' combinations of \eqn{\pm 1}{+/- 1} must be applied to the axes.
##' The function \code{\link{GQdk}} reproduces the entire array of nodes.
##' @seealso \code{\link{GQdk}}
##' @name GQN
##' @docType data
##' @format A list of lists.
##' @examples
##' GQN[[3]][[5]]
##'
NULL

##' University Lecture/Instructor Evaluations by Students at ETH
##' 
##' University lecture evaluations by students at ETH Zurich, anonymized for
##' privacy protection.  This is an interesting \dQuote{medium} sized example of
##' a \emph{partially} nested mixed effect model.
##' 
##' The main goal of the survey is to find \dQuote{the best liked prof},
##' according to the lectures given.  Statistical analysis of such data has been
##' the basis for a (student) jury selecting the final winners.
##' 
##' The present data set has been anonymized and slightly simplified on purpose.
##' 
##' @name InstEval
##' @docType data
##' @format A data frame with 73421 observations on the following 7 variables.
##'   \describe{
##'     \item{\code{s}}{a factor with levels \code{1:2972} denoting
##'       individual students.}
##'     \item{\code{d}}{a factor with 1128 levels from \code{1:2160}, denoting
##'       individual professors or lecturers.}% ("d": \dQuote{Dozierende} in German)
##'     \item{\code{studage}}{an ordered factor with levels \code{2} <
##'       \code{4} < \code{6} < \code{8}, denoting student's \dQuote{age}
##'       measured in the \emph{semester} number the student has been enrolled.}
##'     \item{\code{lectage}}{an ordered factor with 6 levels, \code{1} <
##'       \code{2} < ... < \code{6}, measuring how many semesters back the
##'       lecture rated had taken place.}
##'     \item{\code{service}}{a binary factor with levels \code{0} and
##'       \code{1}; a lecture is a \dQuote{service}, if held for a
##'       different department than the lecturer's main one.}
##'     \item{\code{dept}}{a factor with 14 levels from \code{1:15}, using a
##'       random code for the department of the lecture.}
##' 
##'     \item{\code{y}}{a numeric vector of \emph{ratings} of lectures by
##'       the students, using the discrete scale \code{1:5}, with meanings
##'       of \sQuote{poor} to \sQuote{very good}.}
##'   }
##'   Each observation is one student's rating for a specific lecture
##'   (of one lecturer, during one semester in the past).
##' @keywords datasets
##' @examples
##' 
##' str(InstEval)
##' 
##' head(InstEval, 16)
##' xtabs(~ service + dept, InstEval)
##' 
NULL

##' Paste strength by batch and cask
##' 
##' Strength of a chemical paste product; its quality depending on the delivery
##' batch, and the cask within the delivery.
##' 
##' The data are described in Davies and Goldsmith (1972) as coming from
##' \dQuote{ deliveries of a chemical paste product contained in casks where, in
##' addition to sampling and testing errors, there are variations in quality
##' between deliveries \dots{} As a routine, three casks selected at random from
##' each delivery were sampled and the samples were kept for reference. \dots{}
##' Ten of the delivery batches were sampled at random and two analytical tests
##' carried out on each of the 30 samples}.
##' 
##' @name Pastes
##' @docType data
##' @format A data frame with 60 observations on the following 4 variables.
##'   \describe{
##'     \item{\code{strength}}{paste strength.}
##'     \item{\code{batch}}{delivery batch from which the sample was
##'       sample.  A factor with 10 levels: \sQuote{A} to \sQuote{J}.}
##'     \item{\code{cask}}{cask within the delivery batch from which the
##'       sample was chosen.  A factor with 3 levels: \sQuote{a} to
##'       \sQuote{c}.}
##'     \item{\code{sample}}{the sample of paste whose strength was assayed,
##'       two assays per sample. A factor with 30 levels: \sQuote{A:a} to
##'       \sQuote{J:c}.}
##'   }
##' @source O.L. Davies and P.L. Goldsmith (eds), \emph{Statistical Methods in
##' Research and Production, 4th ed.}, Oliver and Boyd, (1972), section 6.5
##' @keywords datasets
##' @examples
##' str(Pastes)
##' dotplot(cask ~ strength | reorder(batch, strength), Pastes,
##'         strip = FALSE, strip.left = TRUE, layout = c(1, 10),
##'         ylab = "Cask within batch",
##'         xlab = "Paste strength", jitter.y = TRUE)
##' ## Modifying the factors to enhance the plot
##' Pastes <- within(Pastes, batch <- reorder(batch, strength))
##' Pastes <- within(Pastes, sample <- reorder(reorder(sample, strength),
##'           as.numeric(batch)))
##' dotplot(sample ~ strength | batch, Pastes,
##'         strip = FALSE, strip.left = TRUE, layout = c(1, 10),
##'         scales = list(y = list(relation = "free")),
##'         ylab = "Sample within batch",
##'         xlab = "Paste strength", jitter.y = TRUE)
##' ## Four equivalent models differing only in specification
##' (fm1 <- lmer(strength ~ (1|batch) + (1|sample), Pastes))
##' (fm2 <- lmer(strength ~ (1|batch/cask), Pastes))
##' (fm3 <- lmer(strength ~ (1|batch) + (1|batch:cask), Pastes))
##' (fm4 <- lmer(strength ~ (1|batch/sample), Pastes))
##' ## fm4 results in redundant labels on the sample:batch interaction
##' head(ranef(fm4)[[1]])
##' ## compare to fm1
##' head(ranef(fm1)[[1]])
##' ## This model is different and NOT appropriate for these data
##' (fm5 <- lmer(strength ~ (1|batch) + (1|cask), Pastes))
##' 
##' L <- getME(fm1, "L")
##' Matrix::image(L, sub = "Structure of random effects interaction in pastes model")
##' 
NULL

##' Variation in penicillin testing
##' 
##' Six samples of penicillin were tested using the \emph{B. subtilis} plate
##' method on each of 24 plates.  The response is the diameter (mm) of the zone
##' of inhibition of growth of the organism.
##' 
##' The data are described in Davies and Goldsmith (1972) as coming from an
##' investigation to \dQuote{assess the variability between samples of
##' penicillin by the \emph{B. subtilis} method.  I this test method a
##' bulk-innoculated nutrient agar medium is poured into a Petri dish of
##' approximately 90 mm. diameter, known as a plate.  When the medium has set,
##' six small hollow cylinders or pots (about 4 mm. in diameter) are cemented
##' onto the surface at equally spaced intervals.  A few drops of the penicillin
##' solutions to be compared are placed in the respective cylinders, and the
##' whole plate is placed in an incubator for a given time.  Penicillin diffuses
##' from the pots into the agar, and this produces a clear circular zone of
##' inhibition of growth of the organisms, which can be readily measured.  The
##' diameter of the zone is related in a known way to the concentration of
##' penicillin in the solution.}
##' 
##' @name Penicillin
##' @docType data
##' @format A data frame with 144 observations on the following 3 variables.
##'   \describe{
##'     \item{\code{diameter}}{diameter (mm) of the zone of inhibition of
##'       the growth of the organism.}
##'     \item{\code{plate}}{assay plate.  A factor with levels \sQuote{a} to
##'       \sQuote{x}.}
##'     \item{\code{sample}}{penicillin sample.  A factor with levels
##'       \sQuote{A} to \sQuote{F}.}
##'   }
##' @source O.L. Davies and P.L. Goldsmith (eds), \emph{Statistical Methods in
##' Research and Production, 4th ed.}, Oliver and Boyd, (1972), section 6.6
##' @keywords datasets
##' @examples
##' 
##' str(Penicillin)
##' dotplot(reorder(plate, diameter) ~ diameter, Penicillin, groups = sample,
##'         ylab = "Plate", xlab = "Diameter of growth inhibition zone (mm)",
##'         type = c("p", "a"), auto.key = list(columns = 3, lines = TRUE,
##'         title = "Penicillin sample"))
##' (fm1 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin))
##' 
##' L <- getME(fm1, "L")
##' Matrix::image(L, main = "L",
##'               sub = "Penicillin: Structure of random effects interaction")
##' 
NULL

##' Reaction times in a sleep deprivation study
##' 
##' The average reaction time per day for subjects in a sleep deprivation study.
##' On day 0 the subjects had their normal amount of sleep.  Starting that night
##' they were restricted to 3 hours of sleep per night.  The observations
##' represent the average reaction time on a series of tests given each day to
##' each subject.
##' 
##' These data are from the study described in Belenky et al. (2003), for the
##' sleep-deprived group and for the first 10 days of the study, up to the
##' recovery period.
##' 
##' @name sleepstudy
##' @docType data
##' @format A data frame with 180 observations on the following 3 variables.
##'   \describe{
##'     \item{\code{Reaction}}{Average reaction time (ms)}
##'     \item{\code{Days}}{Number of days of sleep deprivation}
##'     \item{\code{Subject}}{Subject number on which the observation was made.}
##'   }
##' @references Gregory Belenky, Nancy J. Wesensten, David R. Thorne, Maria L.
##' Thomas, Helen C. Sing, Daniel P. Redmond, Michael B. Russo and Thomas J.
##' Balkin (2003) Patterns of performance degradation and restoration during
##' sleep restriction and subsequent recovery: a sleep dose-response study.
##' \emph{Journal of Sleep Research} \bold{12}, 1--12.
##' @keywords datasets
##' @examples
##' 
##' str(sleepstudy)
##' xyplot(Reaction ~ Days | Subject, sleepstudy, type = c("g","p","r"),
##'        index = function(x,y) coef(lm(y ~ x))[1],
##'        xlab = "Days of sleep deprivation",
##'        ylab = "Average reaction time (ms)", aspect = "xy")
##' (fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
##' (fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
##' 
NULL


##' Verbal Aggression item responses
##' 
##' These are the item responses to a questionaire on verbal aggression.  These
##' data are used throughout De Boeck and Wilson, \emph{Explanatory Item
##' Response Models} (Springer, 2004) to illustrate various forms of item
##' response models.
##' 
##' 
##' @name VerbAgg
##' @docType data
##' @format A data frame with 7584 observations on the following 13 variables.
##'   \describe{
##'     \item{\code{Anger}}{the subject's Trait Anger score as measured on
##'       the State-Trait Anger Expression Inventory (STAXI)}
##'     \item{\code{Gender}}{the subject's gender - a factor with levels
##'       \code{M} and \code{F}}
##'     \item{\code{item}}{the item on the questionaire, as a factor}
##'     \item{\code{resp}}{the subject's response to the item - an ordered
##'       factor with levels \code{no} < \code{perhaps} < \code{yes}}
##'     \item{\code{id}}{the subject identifier, as a factor}
##'     \item{\code{btype}}{behavior type - a factor with levels
##'       \code{curse}, \code{scold} and \code{shout}}
##'     \item{\code{situ}}{situation type - a factor with levels
##'       \code{other} and \code{self} indicating other-to-blame and self-to-blame}
##'     \item{\code{mode}}{behavior mode - a factor with levels \code{want}
##'       and \code{do}}
##'     \item{\code{r2}}{dichotomous version of the response - a factor with
##'       levels \code{N} and \code{Y}}
##'   }
##' @references De Boeck and Wilson (2004), \emph{Explanatory Item Response
##' Models}, Springer.
##' @source \url{http://bear.soe.berkeley.edu/EIRM/}
##' @keywords datasets
##' @examples
##' 
##' str(VerbAgg)
##' ## Show how  r2 := h(resp) is defined:
##' with(VerbAgg, stopifnot( identical(r2, {
##'      r <- factor(resp, ordered=FALSE); levels(r) <- c("N","Y","Y"); r})))
##' 
##' xtabs(~ item + resp, VerbAgg)
##' xtabs(~ btype + resp, VerbAgg)
##' round(100 * ftable(prop.table(xtabs(~ situ + mode + resp, VerbAgg), 1:2), 1))
##' person <- unique(subset(VerbAgg, select = c(id, Gender, Anger)))
##' if (require(lattice)) {   # is this necessary when the package depends on lattice?
##'     densityplot(~ Anger, person, groups = Gender, auto.key = list(columns = 2),
##'                 xlab = "Trait Anger score (STAXI)")
##' }
##' 
##' \dontrun{## takes about 15 sec
##' print(fmVA <- glmer(r2 ~ (Anger + Gender + btype + situ)^2 +
##'                    (1|id) + (1|item), family = binomial, data =
##'                    VerbAgg), corr=FALSE)
##' }
##'                        ## much faster but less accurate
##' print(fmVA0 <- glmer(r2 ~ (Anger + Gender + btype + situ)^2 +
##'                     (1|id) + (1|item), family = binomial, data =
##'                     VerbAgg, nAGQ=0L), corr=FALSE)
##' 
NULL




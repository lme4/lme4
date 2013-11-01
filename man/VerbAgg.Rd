\docType{data}
\name{VerbAgg}
\alias{VerbAgg}
\title{Verbal Aggression item responses}
\format{A data frame with 7584 observations on the following 13 variables.
  \describe{
    \item{\code{Anger}}{the subject's Trait Anger score as measured on
      the State-Trait Anger Expression Inventory (STAXI)}
    \item{\code{Gender}}{the subject's gender - a factor with levels
      \code{M} and \code{F}}
    \item{\code{item}}{the item on the questionaire, as a factor}
    \item{\code{resp}}{the subject's response to the item - an ordered
      factor with levels \code{no} < \code{perhaps} < \code{yes}}
    \item{\code{id}}{the subject identifier, as a factor}
    \item{\code{btype}}{behavior type - a factor with levels
      \code{curse}, \code{scold} and \code{shout}}
    \item{\code{situ}}{situation type - a factor with levels
      \code{other} and \code{self} indicating other-to-blame and self-to-blame}
    \item{\code{mode}}{behavior mode - a factor with levels \code{want}
      and \code{do}}
    \item{\code{r2}}{dichotomous version of the response - a factor with
      levels \code{N} and \code{Y}}
  }}
\source{
  \url{http://bear.soe.berkeley.edu/EIRM/}
}
\description{
  These are the item responses to a questionaire on verbal
  aggression.  These data are used throughout De Boeck and
  Wilson, \emph{Explanatory Item Response Models}
  (Springer, 2004) to illustrate various forms of item
  response models.
}
\examples{
str(VerbAgg)
## Show how  r2 := h(resp) is defined:
with(VerbAgg, stopifnot( identical(r2, {
     r <- factor(resp, ordered=FALSE); levels(r) <- c("N","Y","Y"); r})))

xtabs(~ item + resp, VerbAgg)
xtabs(~ btype + resp, VerbAgg)
round(100 * ftable(prop.table(xtabs(~ situ + mode + resp, VerbAgg), 1:2), 1))
person <- unique(subset(VerbAgg, select = c(id, Gender, Anger)))
require(lattice)
densityplot(~ Anger, person, groups = Gender, auto.key = list(columns = 2),
            xlab = "Trait Anger score (STAXI)")

if(lme4:::testLevel() >= 3) { ## takes about 15 sec
print(fmVA <- glmer(r2 ~ (Anger + Gender + btype + situ)^2 +
 		   (1|id) + (1|item), family = binomial, data =
		   VerbAgg), corr=FALSE)
}
                       ## much faster but less accurate
print(fmVA0 <- glmer(r2 ~ (Anger + Gender + btype + situ)^2 +
                    (1|id) + (1|item), family = binomial, data =
                    VerbAgg, nAGQ=0L), corr=FALSE)
}
\references{
  De Boeck and Wilson (2004), \emph{Explanatory Item
  Response Models}, Springer.
}
\keyword{datasets}


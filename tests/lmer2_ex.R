stopifnot(suppressPackageStartupMessages(require(lme4)))

## Using simple generated data -- fully balanced here, unbalanced later
set.seed(1)
dat <- within(data.frame(lagoon = factor(rep(1:4, each = 25)),
                         habitat = factor(rep(1:20, each = 5))),
          {   ## a simple  lagoon effect but no random effect
              y <- round(10*rnorm(100, m = 10*as.numeric(lagoon)))
              ## Here, *with* an RE, sigma_a = 100
              RE <- rep(round(rnorm(nlevels(habitat), sd = 100)), each = 5)
              y2 <- y + RE
          })

## FIXME:   want  lmer(* , sparseX = TRUE )  {as in lme4a}
if (FALSE) {                            # need to adapt to new structure

##' <description>
##'
##' <details>
##' @title Comparing the different versions of lmer() for same data & model
##' @param form
##' @param data
##' @param verbose
##' @return
chkLmers <- function(form, data, verbose = FALSE,
                     tol = 200e-7) # had tol = 7e-7 working ..
{
#    m   <- lmer1(form, data = data)  # ok, and more clear
#    m.  <- lmer1(form, data = data, sparseX = TRUE, verbose = verbose)
    m2  <- lmer (form, data = data, verbose = verbose) # lmem-dense
    m2. <- lmer (form, data = data, sparseX = TRUE, verbose = verbose)
    ##
    Eq <- function(x,y) all.equal(x,y, tolerance = tol)
    stopifnot(## Compare  sparse & dense of the new class results
              identical(slotNames(m2), slotNames(m2.))
              ,
              identical(slotNames(m2@fe), slotNames(m2.@fe))
              ,
              Eq(m2@resp, m2.@resp)
              ,
              Eq(m2@re, m2.@re)
              ,
              Eq(m2@fe@coef, m2.@fe@coef)
              ,
              ## and now compare with the "old" (class 'mer')
#              Eq(unname(fixef(m)), m2@fe@beta)
#              ,
#              Eq(unname(fixef(m.)), m2.@fe@beta)
#              ,
              ## to do
              ## all.equal(ranef(m)), m2@re)
              ## all.equal(ranef(m.)), m2.@re)
              TRUE)
    invisible(list(#m=m, m.=m.,
                   m2 = m2, m2. = m2.))
}

chk1 <- chkLmers(y  ~ 0+lagoon + (1|habitat), data = dat, verbose = TRUE)
chk2 <- chkLmers(y2 ~ 0+lagoon + (1|habitat), data = dat, verbose = TRUE)
chk1$m2  ## show( lmer() ) -- sigma_a == 0
chk2$m2. ## show( lmer( <sparseX>) ) --

n <- nrow(dat)
for(i in 1:20) {
    iOut <- sort(sample(n, 1+rpois(1, 3), replace=FALSE))
    cat(i,":  w/o ", paste(iOut, collapse=", ")," ")
    chkLmers(y  ~ 0+lagoon + (1|habitat), data = dat[- iOut,])
    chkLmers(y2 ~   lagoon + (1|habitat), data = dat[- iOut,])
    cat("\n")
}

## One (rare) example where the default tolerance is not sufficient:
dat. <- dat[- c(14, 34, 66, 67, 71, 88),]
try( chkLmers(y ~ 0+lagoon + (1|habitat), data = dat.) )
## Error: Eq(unname(fixef(m)), m2@fe@beta) is not TRUE
##
## but higher tolerance works:
chkLmers(y ~ 0+lagoon + (1|habitat), data = dat., tol = 2e-4, verbose=TRUE)

}
proc.time()
sessionInfo()

## example posted by Stéphane Laurent
## exercises bug where Nelder-Mead min objective function value was >0
set.seed(666)
sims <- function(I, J, sigmab0, sigmaw0){
    Mu <- rnorm(I, mean=0, sd=sigmab0)
    y <- c(sapply(Mu, function(mu) rnorm(J, mu, sigmaw0)))
    data.frame(y=y, group=gl(I,J))
}

I <- 3  # number of groups
J <- 8  # number of repeats per group
sigmab0 <- 0.15  # between standard deviation
sigmaw0 <- 0.15  # within standard deviation

dat <- sims(I, J, sigmab0, sigmaw0)

library(lme4)
fm3 <- lmer(y ~ (1|group), data=dat)
stopifnot(all.equal(unname(unlist(VarCorr(fm3))),
		    if(fm3@optinfo$optimizer == "Nelder_Mead")
		    0.029662844057 else
		    0.02966269809, tolerance = 1e-7))

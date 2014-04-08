library(lme4)
set.seed(101)
                                        # simulation formula
form <- y ~ x + (x|g)
                                        # simulation parameters
tmpPar <- lme4:::mkParsTemplate(form)
tmpPar$beta[] <- c(5, -1)
tmpPar$sigma <- 1
                                        # simulation design
tmpDat <- lme4:::mkDataTemplate(form, nGrp = 100, nPerGrp = 20, rfunc = rnorm)
tmpDat$y <- abs(tmpDat$y)
                                        # 100 simulations
sims <- simulate(form, nsim = 100, newdata = tmpDat, family = Gamma, newparams = tmpPar)

                                        # fit models to simulations and
                                        # catch cases with warnings
                                        # (takes some time)
options(warn = 2)
mods <- lapply(1:100, function(i) {
    simDat <- within(tmpDat, y <- sims[,i])
    m <- try(glmer(form, data = simDat, family = Gamma))
    return(m)
})
options(warn = 0)
wsim <- which(sapply(mods, inherits, "try-error"))

                                        # about a third of the models had warnings
length(wsim)/100

                                        # examine first model with a warning
simDat <- within(tmpDat, y <- sims[,wsim[1]])
simDat0 <- simDat[!is.nan(simDat$y),]
(m <- glmer(form, data = simDat0, family = Gamma))

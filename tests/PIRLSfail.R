## example from HSAUR2 package
 if(FALSE) {
     ## (commented out to avoid R CMD check warning about undeclared dependency
     ## library("multcomp")
     trees513 <- subset(trees513, !species %in%
                        c("fir", "ash/maple/elm/lime", "softwood (other)"))
     trees513$species <- trees513$species[,drop = TRUE]
     levels(trees513$species)[nlevels(trees513$species)] <- "hardwood"
     save("trees513", file="trees513.RData")
 }

library("lme4")
load("trees513.RData")
dfun <- glmer(damage ~ species - 1 + (1 | lattice / plot),
              data = trees513, family = binomial, devFunOnly = TRUE)
ls.str(environment(dfun))# and you can investigate ...

mmod <- glmer(damage ~ species - 1 + (1 | lattice / plot),
              data = trees513, family = binomial())


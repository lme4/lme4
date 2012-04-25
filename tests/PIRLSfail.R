## example from HSAUR2 package
## library("multcomp")
## trees513 <- subset(trees513, !species %in% c("fir", "ash/maple/elm/lime", "softwood (other)"))
## trees513$species <- trees513$species[,drop = TRUE]
## levels(trees513$species)[nlevels(trees513$species)] <- "hardwood"
## save("trees513",file="trees513.RData")
library("lme4")
load("trees513.RData")
## FIXME: try() to pass checks
mmod <- try(lmer(damage ~ species - 1 + (1 | lattice / plot),
              data = trees513, family = binomial()))


library(readr)
library(plyr)
library(dplyr)
dataset <- read_csv("glmernb_foraging_subset.csv") %>%
  mutate(
    name = factor(name),
    fieldid = factor(rank(origarea)),
    logarea = log(origarea)
  )

library(ggplot2); theme_set(theme_bw())
ggplot(dataset,aes(logarea,field_count,colour=habitat))+stat_sum(alpha=0.8)+
    scale_colour_brewer(palette="Dark2")+geom_line(aes(group=name),alpha=0.5)

library(lme4)
library(glmmADMB)
library(glmmTMB)
## devtools::install_github("glmmTMB/glmmTMB",sub="glmmTMB")
library(INLA)
## install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
source(system.file("utils", "allFit.R", package="lme4"))

fn <- "glmernb_foraging_models.rda"
if (file.exists(fn)) {
    load(fn)
} else {
    form <- field_count ~ logarea + habitat + (1 | name) + (1 | fieldid)
    ## fit everything ...
    m.lme4 <- glmer.nb(form, data = dataset, verbose=TRUE)
    m.lme4.nlopt <- glmer.nb(form, data = dataset, verbose=TRUE,
                             control=glmerControl(optimizer="nloptwrap"))
    m.lme4.all <- allFit(m.lme4)  ## does this work on glmer.nb objects?  will be slow ...
    m.admb <- glmmadmb(form, data = dataset, family = "nbinom")
    m.tmb <- glmmTMB(form, data = dataset, family = "nbinom2")
    m.inla <- inla(field_count ~ f(name, model = "iid") + f(fieldid, model = "iid") + logarea + habitat, data = dataset, family = "nbinomial")
    save(list=ls(pattern="m\\."),file=fn)
}

summary(m.lme4)$coefficients  ## bogus SE
summary(m.admb)$coefficients[,1:2]
summary(m.tmb)$coefficients[[1]][,1:2]
m.inla$summary.fixed[,1:2]

library(broom)
library(dotwhisker)
names(m.lme4.all)  <- paste0("lme4_all_",names(m.lme4.all))
modList <- c(list(admb=m.admb,tmb=m.tmb,lme4_basic=m.lme4,lme4_nlopt=m.lme4.nlopt),
             unclass(m.lme4.all))
br <- ldply(modList,tidy) %>% rename(model=.id)
dwplot(br)


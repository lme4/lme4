load("randcrabdata.RData")
## randdata0: simulated data, in form suitable for plotting
## randdata: simulated data, in form suitable for analysis

fr  ## alive/dead formula
fr2 ## proportion alive formula (use with weights=initial.snail.density)

library(ggplot2)
library(grid)
zmargin <- opts(panel.margin=unit(0,"lines"))
theme_update(theme_bw())
g1 <- ggplot(randdata0,aes(x=snail.size,y=surv,colour=snail.size,fill=snail.size))+
    geom_hline(yintercept=1,colour="black")+
    stat_sum(aes(size=factor(..n..)),alpha=0.6)+
    facet_grid(.~ttt)+zmargin+
    geom_boxplot(fill=NA,outlier.colour=NULL,outlier.shape=3)+  ## set outliers to same colour as points
    ## (hard to see which are outliers, but it doesn't really matter in this case)
    scale_size_discrete("# obs",range=c(2,5))


library(lme4)
if (packageVersion("lme4")>"0.999375-42") {
    ## using development lme4 ...
    try(glmer1 <- glmer(fr2,weights=initial.snail.density,family ="binomial", data=randdata))
    ## pwrssUpdate did not converge
    try(glmer1B <- glmer(fr,family ="binomial", data=randdata))
    if (require("lme4.0")) {
        detach("package:lme4")
        ## prop/weights formulation
        glmer1 <- glmer(fr2,weights=initial.snail.density,family ="binomial", data=randdata)
        ## alive/dead formulation
        glmer1B <- glmer(fr,family ="binomial", data=randdata)
        coef(glmer1B)
        detach("package:lme4.0")
    }
    if (require("glmmADMB")) {
        ## prop/weights formulation
        glmer1B <- glmmadmb(fr,family ="binomial", data=randdata)
    }
} else {
    ## CRAN version of lme4
    glmer1 <- glmer(fr2,weights=initial.snail.density,family ="binomial", data=randdata)
    glmer1B <- glmer(fr,family ="binomial", data=randdata)
}
    
fixef(glmer1)
fixef(glmer1B)
## note in this case that variances are not converging to zero
VarCorr(glmer1)

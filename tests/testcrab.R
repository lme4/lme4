library("lme4")
L <- load(system.file("testdata","crabs_randdata2.Rda",package="lme4"))
## randdata0: simulated data, in form suitable for plotting
## randdata: simulated data, in form suitable for analysis

## fr  ## alive/dead formula
## fr2 ## proportion alive formula (use with weights=initial.snail.density)

## FIXME: there are still bigger differences than I'd like between the approaches
## (mostly in the random-effects correlation).  It's not clear who's right;
## lme4 thinks its parameters are better, but ??  Could be explored further.

if (FALSE) {
## library(ggplot2)  ## commented to avoid triggering Suggests: requirement
library(grid)
zmargin <- theme(panel.margin=unit(0,"lines"))
theme_set(theme_bw())
g1 <- ggplot(randdata0,aes(x=snail.size,y=surv,colour=snail.size,fill=snail.size))+
    geom_hline(yintercept=1,colour="black")+
    stat_sum(aes(size=factor(..n..)),alpha=0.6)+
    facet_grid(.~ttt)+zmargin+
    geom_boxplot(fill=NA,outlier.colour=NULL,outlier.shape=3)+  ## set outliers to same colour as points
    ## (hard to see which are outliers, but it doesn't really matter in this case)
    scale_size_discrete("# obs",range=c(2,5))
}

t1 <- system.time(glmer1 <- glmer(fr2,weights=initial.snail.density,
                                  family ="binomial", data=randdata))
t1B <- system.time(glmer1B <- glmer(fr,family ="binomial", data=randdata))

res1 <- c(fixef(glmer1),c(VarCorr(glmer1)$plot))
res1B <- c(fixef(glmer1B),c(VarCorr(glmer1B)$plot))
p1 <- unlist(getME(glmer1,c("theta","beta")))
stopifnot(all.equal(res1,res1B))

dfun <- update(glmer1,devFunOnly=TRUE)
stopifnot(all.equal(dfun(p1),c(-2*logLik(glmer1))))
##
## library(lme4.0)  ## version 0.999999.2 results
## t1_lme4.0 <- system.time(glmer1X <-
##                           glmer(fr2,weights=initial.snail.density,
##                                 family ="binomial", data=randdata))
## dput(c(fixef(glmer1X),c(VarCorr(glmer1X)$plot)))
## p1X <- c(getME(glmer1X,"theta"),getME(glmer1X,"beta"))
p1X <- c(0.681301656652347, -1.14775239687404, 0.436143018123226, 
2.77730476938968, 0.609023583738824, -1.60055813739844, 2.0324468778545, 
0.624173873057839, -1.7908793509579, -2.44540201631615, -1.42365990002708, 
-2.26780929006268, 0.700928084600075, -1.26220238391029, 0.369024582097804, 
3.44325347343035, 2.26400391093108)
stopifnot(all.equal(unname(p1),p1X,tolerance=0.03))
dfun(p1X)
dfun(p1)
## ~ 1.8 seconds elapsed time
lme4.0_res <- 
    structure(c(2.77730476938968, 0.609023583738824, -1.60055813739844, 
                2.0324468778545, 0.624173873057839, -1.7908793509579, -2.44540201631615, 
                -1.42365990002708, -2.26780929006268, 0.700928084600075, -1.26220238391029, 
                0.369024582097804, 3.44325347343035, 2.26400391093108, 0.464171947357232, 
                -0.532754465140956, -0.532754465140956, 0.801690946568518),
          .Names = c("(Intercept)", 
          "crab.speciesS", "crab.speciesW", "crab.sizeS", "crab.sizeM", 
          "snail.sizeS", "crab.speciesS:crab.sizeS", "crab.speciesS:crab.sizeM", 
          "crab.speciesS:snail.sizeS", "crab.speciesW:snail.sizeS",
          "crab.sizeS:snail.sizeS", "crab.sizeM:snail.sizeS",
          "crab.speciesS:crab.sizeS:snail.sizeS", 
          "crab.speciesS:crab.sizeM:snail.sizeS", "", "", "", ""))

stopifnot(all.equal(res1,lme4.0_res,tolerance=0.015))


## library("glmmADMB")
## prop/weights formulation: ~ 7 seconds
## t1_glmmadmb <- system.time(glmer1B <- glmmadmb(fr,family ="binomial",
##                                               corStruct="full",data=randdata))
## dput(c(fixef(glmer1B),c(VarCorr(glmer1B)$plot)))
glmmADMB_res <- structure(c(2.7773101267224, 0.609026276823218, -1.60055704634712, 
2.03244174458562, 0.624171008585953, -1.79088398816641, -2.44540300134182, 
-1.42366043619683, -2.26780858382505, 0.700927141726545, -1.26219964572264, 
0.369029052442189, 3.44326297908383, 2.26403738918967, 0.46417, 
-0.53253, -0.53253, 0.80169), .Names = c("(Intercept)", "crab.speciesS", 
"crab.speciesW", "crab.sizeS", "crab.sizeM", "snail.sizeS", "crab.speciesS:crab.sizeS", 
"crab.speciesS:crab.sizeM", "crab.speciesS:snail.sizeS", "crab.speciesW:snail.sizeS", 
"crab.sizeS:snail.sizeS", "crab.sizeM:snail.sizeS", "crab.speciesS:crab.sizeS:snail.sizeS", 
"crab.speciesS:crab.sizeM:snail.sizeS", "", "", "", ""))

stopifnot(all.equal(res1B,glmmADMB_res,tolerance=0.015))

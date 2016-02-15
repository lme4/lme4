## original code for reading/aggregating:

## tickdata <- read.table("Elston2001_tickdata.txt",header=TRUE,
##   colClasses=c("factor","numeric","factor","numeric","factor","factor"))

## tickdata <- transform(tickdata,cHEIGHT=scale(HEIGHT,scale=FALSE))
## for (i in names(tickdata)) {
##   if (is.factor(tickdata[[i]])) {
##     tickdata[[i]] <- factor(tickdata[[i]],levels=sort(as.numeric(levels(tickdata[[i]]))))
##   }
## }
## summary(tickdata)
## grouseticks <- tickdata

## library(reshape)
## meantick <- rename(aggregate(TICKS~BROOD,data=tickdata,FUN=mean),
##                    c(TICKS="meanTICKS"))
## vartick <- rename(aggregate(TICKS~BROOD,data=tickdata,FUN=var),
##                    c(TICKS="varTICKS"))
## uniqtick <- unique(subset(tickdata,select=-c(INDEX,TICKS)))
## grouseticks_agg <- Reduce(merge,list(meantick,vartick,uniqtick))

## save("grouseticks","grouseticks_agg",file="grouseticks.rda")

library("lme4")
data(grouseticks)
do.plots <- FALSE
form <- TICKS~YEAR+HEIGHT+(1|BROOD)+(1|INDEX)+(1|LOCATION)

## fit with lme4
## library(lme4)
## t1 <- system.time(full_mod1  <- glmer(form, family="poisson",data=grouseticks))
## c1 <- c(fixef(full_mod1),unlist(VarCorr(full_mod1)), logLik=logLik(full_mod1),time=t1["elapsed"])
## allcoefs1 <- c(unlist(full_mod1@ST),fixef(full_mod1))
## detach("package:lme4")

## lme4 summary results:
t1 <- structure(c(1.288, 0.048, 1.36, 0, 0), class = "proc_time",
                .Names = c("user.self",
                  "sys.self", "elapsed", "user.child", "sys.child"))

c1 <- structure(c(11.3559080756861, 1.1804105508475, -0.978704335712111,
                  -0.0237607330254979, 0.293232458048324, 0.562551624933584,
                  0.279548178949372,
                  -424.771990224991, 1.36),
                .Names = c("(Intercept)", "YEAR96",
                  "YEAR97", "HEIGHT", "INDEX", "BROOD", "LOCATION",
                  "logLik", "time.elapsed"
                  ))
allcoefs1 <- structure(c(0.541509425632023, 0.750034415832756,
                         0.528723159081737,
                         11.3559080756861, 1.1804105508475,
                         -0.978704335712111, -0.0237607330254979
                         ),
                       .Names = c("", "", "", "(Intercept)",
                         "YEAR96", "YEAR97",  "HEIGHT"))

sfun <- function(x,time) {
    res <- c(fixef(x),unlist(VarCorr(x)),logLik=logLik(x),time=NA)
    if (!missing(time)) res["time"] <- time["elapsed"]
    res
}

allcoefs <- function(x) c(getME(x,"theta"),getME(x,"beta"))
    
if (lme4:::testLevel() > 1) {
    t2 <- system.time(full_mod2  <- glmer(form, family="poisson",
                                          data=grouseticks,
                control=glmerControl(optCtrl=list(maxfun=60000))))
    ## back-compatible
    t3 <- system.time(full_mod3  <- glmer(form, family="poisson",
                                          data=grouseticks,
                control=glmerControl(optimizer=c("bobyqa","Nelder_Mead"))))
    
    c2 <- sfun(full_mod2)
    c3 <- sfun(full_mod3)
    ## refit

    ## FIXME: with bobyqa, *fails* with non-pos-def Cholmod ...
    full_mod4 <- try(refit(full_mod2,grouseticks$TICKS))

    full_mod5 <- refit(full_mod3,grouseticks$TICKS)
    stopifnot(all.equal(c2,c3,tolerance=1e-6))
    stopifnot(all.equal(c3,sfun(full_mod5),tolerance=1e-6))
    stopifnot(all.equal(unname(allcoefs(full_mod3)),unname(allcoefs1),
              tolerance=1e-5))
}


## deviance function
## FIXME: does compDev do _anything_ any more?
mm <- glmer(form, family="poisson",data=grouseticks,devFunOnly=TRUE)
mm2 <- glmer(form, family="poisson",data=grouseticks,
             devFunOnly=TRUE,control=glmerControl(compDev=TRUE))
stopifnot(all.equal(mm(allcoefs1),mm2(allcoefs1)))


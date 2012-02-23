tickdata <- read.table("Elston2001_tickdata.txt",header=TRUE,
  colClasses=c("factor","numeric","factor","numeric","factor","factor"))

tickdata <- transform(tickdata,cHEIGHT=scale(HEIGHT,scale=FALSE))

form <- TICKS~YEAR+HEIGHT+(1|BROOD)+(1|INDEX)+(1|LOCATION)
## fit with lme4
library(lme4)
t1 <- system.time(full_mod1  <- glmer(form, family="poisson",data=tickdata))
c1 <- c(fixef(full_mod1),unlist(VarCorr(full_mod1)), logLik=logLik(full_mod1),time=t1["elapsed"])
allcoefs1 <- c(unlist(full_mod1@ST),fixef(full_mod1))
detach("package:lme4")

## fit with lme4Eigen
library(lme4Eigen)
t2 <- system.time(full_mod2  <- glmer(form, family="poisson",data=tickdata))
c2 <- c(fixef(full_mod2),unlist(VarCorr(full_mod2)), logLik=logLik(full_mod2),time=t2["elapsed"])

allcoefs <- function(x) c(getME(x,"theta"),getME(x,"beta"))

## deviance function
mm <- glmer(form, family="poisson",data=tickdata,devFunOnly=TRUE,tolPwrss=1e-12,verbose=4,
            compDev=FALSE)
mm(allcoefs1)
## works with compDev=FALSE, fails with compDev=TRUE

## refit
full_mod3 <- refit(full_mod2,tickdata$TICKS)

fn <- "elston_fits.RData"
## FIXME::: some possibility of differing results? 1780.
## what changed ???
if (!file.exists(fn)) {
  tvec <- seq(-13,-7,by=0.1)
  dmat <- matrix(nrow=length(tvec),ncol=9,  
                 dimnames=list(NULL,c("deviance","time_elapsed",
                   paste("theta",1:3,sep=""),paste("beta",1:4,sep=""))))
  for (i in seq_along(tvec)) {
    tt <- system.time(gg <- glmer(form,family="poisson",data=tickdata,tolPwrss=10^tvec[i]))["elapsed"]
    cat(i,tvec[i],deviance(gg),"\n")
    dmat[i,] <- c(deviance(gg),tt,allcoefs(gg))
  }
  dmat <- data.frame(logtol=tvec,dmat)
  save("dmat",file=fn)
} else load(fn)

newdev <- apply(dmat[,-(1:3)],1,mm)
library(ggplot2)
library(reshape)
qplot(logtol,value,data=melt(dmat,id.var="logtol"),geom=c("line"))+
  facet_wrap(~variable,scale="free")+
  geom_line(data=data.frame(logtol=dmat$logtol,value=newdev,variable="deviance"),
            colour="blue")+
  geom_hline(data=data.frame(variable="deviance",val=mm(allcoefs1)),
             aes(yintercept=val),colour="red")+theme_bw()+
  geom_vline(data=data.frame(variable="deviance",val=-10),
             aes(xintercept=val),colour="gray")

detach("package:lme4Eigen")
cbind(c1,c2)


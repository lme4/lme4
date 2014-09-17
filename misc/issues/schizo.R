## data from http://www.uic.edu/~hedeker/SCHIZX1.DAT.txt
## inspiration from http://www.uic.edu/~hedeker/long.html

fn <- "schizx1.dat"
if (!file.exists(fn)) {
    download.file("http://www.uic.edu/~hedeker/SCHIZX1.DAT.txt",dest=fn)
}
## almost works: get extra ^Z, need to strip it ...
dd0 <- read.table("schizx1.dat")
names(dd0) <- c("id","imps79","imps79b","imps79o","int","tx",
               "week","sweek","txswk")
apply(dd0,2,function(x) sum(x==-9,na.rm=TRUE))

## dangerous to use na.strings==c("-9","NA") in general?
## what if there are legitimate -9 values in *some* columns?
na.vals <- function(x,crit) {
    x[crit] <- NA
    return(x)
}
library(plyr)  ## for mutate()
library(MASS)  ## for glmmPQL
library(lme4)
library(ggplot2); theme_set(theme_bw())
dd <- mutate(dd0,
                imps79=na.vals(imps79,imps79<0),
                imps79b=na.vals(imps79b,imps79b<0),
                imps79b=factor(imps79b,labels=c("le mild","ge moderate")),
                tx=factor(tx,labels=c("placebo","drug")))

## sweek is square-root of week ...
ggplot(dd,aes(sweek,imps79b,colour=tx,group=id))+geom_line()+
    stat_sum(aes(group=tx))+
    geom_smooth(method="glm",family="binomial",se=FALSE)

ggplot(dd,aes(sweek,imps79b,colour=tx,group=id))+geom_line()+
    facet_wrap(~id)

with(dd,table(tx,sweek))
with(na.omit(dd),table(tx,sweek))

## reduce this to sequences: most data are observed at only 4 time periods
dd$rweek <- round(dd$sweek^2)
dd2 <- na.omit(subset(dd,rweek %in% c(0,1,3,6)))
nrow(na.omit(dd))-nrow(dd2)  ## lose only 34 cases (out of 1600)
dd3 <- ddply(dd2,"id",summarise,
             tx=tx[1],
             trans=paste(imps79b,collapse=""))
ttab <- with(dd3,table(trans,tx))
ttab[rowSums(ttab)>10,]

m0 <- glm(imps79b~tx*sweek,dd,family=binomial)
m1 <- MASS::glmmPQL(imps79b~tx*sweek,random=~1|id,dd,family=binomial)
m2 <- MASS::glmmPQL(imps79b~tx*sweek,random=~1+sweek|id,dd,family=binomial)
m3 <- glmer(imps79b~tx*sweek+(1|id),dd,family=binomial)
## extremely slow: lots of 'step-halving' action
m4 <- glmer(imps79b~tx*sweek+(1+sweek|id),dd,family=binomial,
            verbose=100)
## eventually fails to converge with maxgrad=180.4 (!!)
library(MCMCglmm)
m5 <- MCMCglmm(imps79b~tx*sweek,
               random=~id,
               data=na.omit(dd),
               family="categorical",verbose=FALSE)
m6 <- MCMCglmm(imps79b~tx*sweek,
               random=~id,
               data=na.omit(dd),
               family="categorical",verbose=FALSE)
m6 <- MCMCglmm(imps79b~tx*sweek,
               random=~us(sweek):id,
               data=na.omit(dd),
               family="categorical",verbose=FALSE)

cbind(glm=coef(m0),glmmPQL_int=fixef(m1),
      glmmPQL_slope=fixef(m2),lme4_int=fixef(m3),
      lme4_slope=fixef(m4))

library(coefplot2)
coefplot2(list(glm=m0,PQL_int=m1,PQL_slope=m2,
               Lapl_int=m3,Lapl_slope=m4),
          legend=TRUE)
## needs more work to combine m5
library(RColorBrewer)
cvec <- brewer.pal(7,"Dark2")
coefplot2(list(glm=m0,PQL_int=m1,PQL_slope=m2,
               Lapl_int=m3,Lapl_slope=m4,
               MCMC_int=m5,MCMC_slope=m6),
          col=cvec,
          legend=TRUE)

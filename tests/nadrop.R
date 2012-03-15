library(lme4)
d <- data.frame(x=runif(100),f=factor(rep(1:10,10)))
set.seed(101)
u <- rnorm(10)
d <- transform(d,y=rnorm(100,1+2*x+u[f],0.2))
d0 <- d
d[c(3,5,7),"x"] <- NA

## 'omit' and 'exclude' are the only choices under which
##  we will see NA values in the results
fm0 <- lmer(y~x+(1|f),data=d0)
fm1 <- lmer(y~x+(1|f),data=d)
fm2 <- lmer(y~x+(1|f),data=d,na.action="na.exclude")
try(fm3 <- lmer(y~x+(1|f),data=d,na.action="na.pass"))

refit(fm0)
refit(fm1)
refit(fm2)

refit(fm0,runif(100))
refit(fm1,runif(100))
refit(fm2,runif(100))

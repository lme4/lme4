## private data: request access from Alex Whitworth for debugging purposes
load("PXX_BenBolkerExDat.Rdata")  

library(lme4)
## source("./F01_predict_funcs.R")

## this is the problematic statement:
predict(mod1, dat1, allow.new.levels=TRUE)
## -> Error in X %*% fixef(object) : non-conformable arguments

(form1 <- formula(mod1))
## log(physical + 2000) ~ 0 + geography:format2 + geography:yr_since_rel + 
##     geography:xmas + (1 | geography/title)

## count levels ...
with(model.frame(mod1),{
         nf2 <- length(unique(format2))
         ng <- length(levels(geography))
         ng*nf2 + 2*ng  ## 40
     })
## the formula implies 40 parameters (ng*nf2 + 2*ng) but X only has
## 36 due to missing combinations.  BUT ...
## model matrix does *not* have any dropped columns recorded (why not???)
attr(getME(mod1,"X"),"col.dropped") ## NULL

dat0 <- model.frame(mod1,fixed.only=TRUE)

## there are three combinations missing here:
i1 <- with(dat0,levels(interaction(geography,format2,drop=TRUE)))
length(i1) ## 21
i2 <- with(dat0,levels(interaction(geography,format2)))
length(i2) ## 24
setdiff(i2,i1)
## [1] "UMG_FRA.SuperAudio" "UMG_GER.SuperAudio" "JAPAN.Vinyl"

## recreate logic from inside predict.merMod, for convenience
RHS <- formula(mod1,fixed.only=TRUE)[-2]
object <- mod1
newdata <- dat1
X0 <- getME(object, "X")
Terms <- terms(object,fixed.only=TRUE)
mf <- model.frame(object, fixed.only=TRUE)
isFac <- vapply(mf, is.factor, FUN.VALUE=TRUE)
isFac[attr(Terms,"response")] <- FALSE
orig_levs <- if (length(isFac)==0) NULL else lapply(mf[isFac],levels)
mfnew <- model.frame(delete.response(Terms),
                     newdata,
                     xlev=orig_levs)
X <- model.matrix(RHS, data=mfnew,
                  contrasts.arg=attr(X0,"contrasts"))
c0 <- colnames(X0)  ## original: 36
c1 <- colnames(X)  ## new : 40
setdiff(c1,c0)
## these are the missing values from above, plus one additional
## geography*xmas interaction level
##
## [1] "geographyUMG_FRA:format2SuperAudio" "geographyUMG_GER:format2SuperAudio"
## [3] "geographyJAPAN:format2Vinyl"        "geographyLT:xmasTRUE"

## this could fix it, but we shouldn't have to?
X <- X[,colnames(X0)]
       
## attempt at a reproducible example
dd <- expand.grid(f=LETTERS[1:3],g=letters[1:3],h=1:3,rep=1:10)
dd <- subset(dd,!(f=="A" & g=="a"))
set.seed(101)
dd$y <- rnorm(nrow(dd))

m1 <- lmer(y~f:g+(1|h),dd)
newdd <- subset(dd,rep=1:5)
predict(m1,newdd)  ## works fine without hack because X has col.dropped
                   ## attribute ...

## refit with extracted data
dat2 <- model.frame(mod1)
names(dat2)[1] <- "y"  ## response has an awkward name
form2 <- y ~ 0 + geography:format2 + geography:yr_since_rel + 
    geography:xmas + (1 | geography/title)
m2 <- lmer(form2,data=dat2)
## fixed-effect model matrix is rank deficient
##     so dropping 4 columns / coefficients
attr(getME(m2,"X"),"col.dropped")
predict(m2,newdata=dat1,allow.new.levels=TRUE) ## works

## an attempt at a reproducible example of @bachlaw's issue
set.seed(101)
## construct data frame that has one level with only a single obs.
## so slope can't be estimated for that level (induce rank-deficiency)
d <- data.frame(y=rnorm(60),x=rnorm(60),f=factor(rep(1:3,c(1,29,30))),
                g=factor(rep(1:6,each=10)))
m1 <- lm(y~f*x,data=d)
coef(m1)
## f3:x is NA

library("lme4")
m2 <- lmer(y~f*x+(1|g),d)
predict(m2,re.form=~1|g)
predict(m2)


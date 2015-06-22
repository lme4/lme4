## https://github.com/lme4/lme4/issues/59
library(lme4)
dat <- read.csv(system.file("testdata","dat20101314.csv",package="lme4"))
NMcopy <- lme4:::Nelder_Mead

cc <- capture.output(lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
                          control=
                          lmerControl("NMcopy", 
                                      optCtrl= list(iprint=20))))

## check that printing goes through step 140 twice and up to 240 once
## findStep <- function(str,n) sum(grepl(paste0("\\(NM\\) ",n,": "),str))
cc <- paste(cc,collapse="")
countStep <- function(str,n) {
    length(gregexpr(paste0("\\(NM\\) ",n,": "),str)[[1]])
}

stopifnot(countStep(cc,140)==2 && countStep(cc,240)==1)

## testStr <-
## "(NM) 20: f = -53.3709 at 0.706667 0.813333  1.46444(NM) 40: f = -147.132 at     0     0 19.18(NM) 60: f = -147.159 at       0       0 17.4275(NM) 80: f = -147.159 at       0       0 17.5615(NM) 100: f = -147.159 at    0       0 17.5754(NM) 120: f = -147.159 at       0       0 17.5769(NM) 140: f = -147.159 at       0       0 17.5768(NM) 20: f = -165.55 at 0.0933333  0.573333   17.3168(NM) 40: f = -173.704 at 0.23799  1.4697 16.9728(NM) 60: f = -173.849 at 0.449634  1.39998  16.9452(NM) 80: f = -174.421 at 0.52329 1.69123 18.1534(NM) 100: f = -176.747 at 0.762043  1.88271  32.8993(NM) 120: f = -176.839 at 0.751206  1.75371 37.2128(NM) 140: f = -176.853 at 0.706425   1.7307  35.7528(NM) 160: f = -176.853 at 0.710803  1.73476  35.7032(NM) 180: f = -176.853 at 0.710159  1.73449  35.6699(NM) 200: f = -176.853 at 0.710271  1.73461 35.6689(NM) 220: f = -176.853 at 0.710259   1.7346  35.6684(NM) 240: f = -176.853 at 0.710257  1.73459  35.6685Linear mixed model fit by REML ['lmerMod']"

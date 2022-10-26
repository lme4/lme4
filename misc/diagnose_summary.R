f <- function(L) {
    load(sprintf("test-summary_testlevel_%d.rda", L))
    return(mget(ls()))
}

r1 <- f(1)
r100 <- f(100)
waldo::compare(r1$m1, r100$m1)
summary(r1$m1)
summary(r100$m1)

## cc1 doesn't include 

waldo::compare(r1$ss, r100$ss)
waldo::compare(r1$m1, r100$m1, ignore_attr = TRUE)

setdiff(names(r100$ss$loadedOnly), names(r1$ss$loadedOnly))  ## MASS, boot, splines
setdiff(names(r100$ss$otherPkgs), names(r1$ss$otherPkgs))  ## MASS, boot, splines
## what is messing things up? stats4, gamm4?

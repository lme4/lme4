## utility for processing times from win-builder logs etc.

tmpdir <- "1yqHv23qSIGU"
tmpdir <- "JJVaXb6jq9AR"
head_url <- sprintf("https://win-builder.r-project.org/%s/",tmpdir)

## extract times listed in square brackets
get_brackets <- function(r) as.numeric(gsub("^.*\\[([0-9]+)s\\].*","\\1",r))
get_testname <- function(r) gsub("^ *Running '([^']+)'.*$","\\1",r)

r0 <- readLines(paste0(head_url,"00check.log"))
r <- r0[grep("[",r0,fixed=TRUE)]
r2 <- r[grep("Running",r,fixed=TRUE)]
r <- r[grep("Running",r,invert=TRUE)]


r3 <- r2[order(-get_brackets(r2),get_testname(r2))]
View(r3)
n <- get_brackets(r)
r[order(-n)]
print(full_time <- sum(n)/60)
print(test_time <- sum(get_brackets(r2))/60)

## full time: 664 seconds
## need to drop another 90 seconds or so?
## offset (9)  + hatvalues (4) + optimizer (8) + predsim (8)
##  + priorWeightsModComp (8)  + testcolonizer (4) + varcorr(4) +
##  drop (6s) + drop1contrasts (6)
## 
## old
## skipped stuff
skipped <- c("glmerControlPass",
             "glmmWeights",
             "boundary",
             "glmmExt",
             "bootMer",
             "priorWeights",
             "respiratory",
             "testcrab",
             "simulate",
             "lmer-1",
             "confint",
             "glmer-1",
             "throw",
             "agridat_gotway")

r3 <- r2[grep(paste(skipped,collapse="|"),r2,invert=TRUE)]
new_test_time <- sum(get_brackets(r3))/60
test_time-new_test_time

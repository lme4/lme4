library(lme4)
library(reshape2)
library(dplyr)

##' Construct a lookup table for the random effects of a \code{merMod}
##' object
##'
##' @param object a \code{merMod} object
##' @return a data frame containing one row per random effect, and
##' with four columns:
##'     
##' (1) index -- the index of the random effects, which may be used to
##' subset the columns of Z, the rows/columns of Lambda, and the
##' rows/columns of the conditional variance matrix,
##' 
##' (2) level -- the level of the grouping factor of the random
##' effect,
##' 
##' (3) variable -- the explanatory variable associated with the
##' random effect, and
##' 
##' (4) group -- the grouping factor of the random effect
mkReIndex <- function(object) {
                                        # get information about the
                                        # dimensions of the random
                                        # effects
    rp <- rePos$new(object)
                                        # list with the levels for
                                        # each grouping factor (one
                                        # character-valued list
                                        # element per factor)
    levs <- lapply(rp$flist, levels)
                                        # names of the variables
                                        # associated with each random
                                        # effect
    variableNames <- 
        mapply(rep, rp$cnms,
               times = rp$nlevs[attr(rp$flist, "assign")],
               SIMPLIFY = FALSE) %>%
        unlist(use.names = FALSE)
                                        # construct the output data
                                        # frame
    rep %>%
        mapply(levs[attr(rp$flist, "assign")], # levels associated with each RE term
               each = rp$ncols,                # num vars associated with each RE term
               SIMPLIFY = FALSE) %>%
        melt() %>% 
        setNames(c("level", "group")) %>%
        mutate(variable = variableNames) %>%  
        mutate(index = row_number()) %>%
        select(index, level, variable, group)
}

mkRanefTable <- function(object, saveCondVarMatrix = FALSE) {
    ## q table
    re <- getME(object, "b")
    cv <- lme4:::condVar(object)
    out <- 
        object %>%
            mkReIndex() %>%
            mutate(ranef = as.numeric(re), condVar = as.numeric(diag(cv))) %>%
            as_data_frame
    attr(out, "condVarMatrix") <- cv
    return(out)
}

mkFixefTable <- function(object) {
    ## p table
    object %>%
        fixef() %>%
        melt(value.name = "fixef") %>%
        add_rownames("variable") %>%
        mutate(fixefVar = diag(vcov(object)))
}

mkRanVCTable <- function(object) {
    VarCorr(object) %>%
        setNames(names(getME(object, "cnms"))) %>%
        as.data.frame() %>%
        setNames(c("group", "variable", "variableCor", "ranefVar", "ranefSdCor")) %>%
        as_data_frame
}

data("sleepstudy")
fm1 <- lmer(Reaction ~ Days + (0 + Days | Subject) + (1 | Subject), data = sleepstudy)
fm2 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)

## random effect tables
##
## the index column indexes ...
##   (1) the elements of the random effects vector
##   (2) the columns of the Z matrix
##   (3) the rows and columns of the Lambda matrix and the condVar matrix
mkRanefTable(fm1)
mkRanefTable(fm2)

## can be used to return only part of the condVar matrix,
## say Subject 331 ...
(ii1 <-
    fm1 %>%
    mkRanefTable() %>%
    filter(grepl("331", level)))
(ii2 <-
    fm2 %>%
    mkRanefTable() %>%
    filter(grepl("331", level)))

cv1 <- lme4:::condVar(fm1)[ii1$index, ii1$index]
dimnames(cv1) <- list(ii1$variable, ii1$variable)
cv2 <- lme4:::condVar(fm2)[ii2$index, ii2$index]
dimnames(cv2) <- list(ii2$variable, ii2$variable)
cv1
cv2

## can tack on fixed effect or variance-covariance information
fm2 %>%
    mkRanefTable() %>%
    left_join(fm2 %>%
              mkRanVCTable %>%
              filter(is.na(variableCor)) %>%
              select(-variableCor, -ranefSdCor),
              by = c("group", "variable")) %>%
    left_join(mkFixefTable(fm2),
              by = "variable") %>%
    mutate(totef = fixef + ranef) %>%
    mutate(totefVar = fixefVar + condVar)


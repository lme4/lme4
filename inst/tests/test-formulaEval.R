library("testthat")
library("lme4")

context("data= argument and formula evaluation")

test_that("glmerFormX", {
    set.seed(101)

    n <- 50
    x <- rbinom(n, 1, 1/2)
    y <- rnorm(n)
    z <- rnorm(n)
    r <- sample(1:5, size=n, replace=TRUE)
    d <- data.frame(x,y,z,r)

    F <- "z"
    rF <- "(1|r)"
    modStr <- (paste("x ~", "y +", F, "+", rF))
    modForm <- as.formula(modStr)

    expect_that(m_data.3 <- glmer( modStr , data=d, family="binomial"), is_a("glmerMod"))
    expect_error(drop1(m_data.3),"'data' not found")
    expect_that(m_data.4 <- glmer( "x ~ y + z + (1|r)" , data=d, family="binomial"), is_a("glmerMod"))
    expect_error(drop1(m_data.4),"'data' not found")
})

test_that("glmerForm", {
    set.seed(101)

    n <- 50
    x <- rbinom(n, 1, 1/2)
    y <- rnorm(n)
    z <- rnorm(n)
    r <- sample(1:5, size=n, replace=TRUE)
    d <- data.frame(x,y,z,r)

    F <- "z"
    rF <- "(1|r)"
    modStr <- (paste("x ~", "y +", F, "+", rF))
    modForm <- as.formula(modStr)

    ## formulas have environments associated, but character vectors don't
    ## data argument not specified:
    ## should work, but documentation warns against it
    expect_that(m_nodata.0 <- glmer( x ~ y + z + (1|r) , family="binomial"), is_a("glmerMod"))
    expect_that(m_nodata.1 <- glmer( as.formula(modStr) , family="binomial"), is_a("glmerMod"))
    expect_that(m_nodata.2 <- glmer( modForm , family="binomial"), is_a("glmerMod"))
    expect_that(m_nodata.3 <- glmer( modStr , family="binomial"), is_a("glmerMod"))
    expect_that(m_nodata.4 <- glmer( "x ~ y + z + (1|r)" , family="binomial"), is_a("glmerMod"))

    ## apply drop1 to all of these ...
    m_nodata_List <- list(m_nodata.0,m_nodata.1,m_nodata.2,m_nodata.3,m_nodata.4)
    d_nodata_List <- lapply(m_nodata_List,drop1)

    rm(list=c("x","y","z","r"))

    ## data argument specified
    expect_that(m_data.0 <- glmer( x ~ y + z + (1|r) , data=d, family="binomial"), is_a("glmerMod"))
    expect_that(m_data.1 <- glmer( as.formula(modStr) , data=d, family="binomial"), is_a("glmerMod"))
    expect_that(m_data.2 <- glmer( modForm , data=d, family="binomial"), is_a("glmerMod"))
    expect_that(m_data.3 <- glmer( modStr , data=d, family="binomial"), is_a("glmerMod"))
    expect_that(m_data.4 <- glmer( "x ~ y + z + (1|r)" , data=d, family="binomial"), is_a("glmerMod"))

    ff <- function() {
        set.seed(101)
        n <- 50
        x <- rbinom(n, 1, 1/2)
        y <- rnorm(n)
        z <- rnorm(n)
        r <- sample(1:5, size=n, replace=TRUE)
        d2 <- data.frame(x,y,z,r)
        glmer( x ~ y + z + (1|r), data=d2, family="binomial")
    }
    m_data.5 <- ff()

    ff2 <- function() {
        set.seed(101)
        n <- 50
        x <- rbinom(n, 1, 1/2)
        y <- rnorm(n)
        z <- rnorm(n)
        r <- sample(1:5, size=n, replace=TRUE)
        glmer( x ~ y + z + (1|r), family="binomial")
    }
    m_data.6 <- ff2()
   

    m_data_List <- list(m_data.0,m_data.1,m_data.2,m_data.3,m_data.4,m_data.5,m_data.6)
    badNums <- 4:5
    d_data_List <- lapply(m_data_List[-badNums],drop1)

    ## these do NOT fail if there is a variable 'd' living in the global environment --
    ## they DO fail in the testthat context
    expect_error(drop1(m_data.3),"'data' not found")
    expect_error(drop1(m_data.4),"'data' not found")
    
    ## expect_error(lapply(m_data_List[4],drop1))
    ## expect_error(lapply(m_data_List[5],drop1))
    ## d_data_List <- lapply(m_data_List,drop1,evalhack="parent")  ## fails on element 1
    ## d_data_List <- lapply(m_data_List,drop1,evalhack="formulaenv")  ## fails on element 4
    ## d_data_List <- lapply(m_data_List,drop1,evalhack="nulldata")  ## succeeds
    ## drop1(m_data.5,evalhack="parent") ## 'd2' not found
    ## drop1(m_data.5,evalhack="nulldata") ## 'x' not found (d2 is in environment ...)
    ## should we try to make update smarter ... ??

    ## test equivalence of (i vs i+1) for all models, all drop1() results
    for (i in 1:(length(m_nodata_List)-1)) {
        expect_equivalent(m_nodata_List[[i]],m_nodata_List[[i+1]])
        expect_equivalent(d_nodata_List[[i]],d_nodata_List[[i+1]])
    }

    expect_equivalent(m_nodata_List[[1]],m_data_List[[1]])
    expect_equivalent(d_nodata_List[[1]],d_data_List[[1]])

    for (i in 1:(length(m_data_List)-1)) {
        expect_equivalent(m_data_List[[i]],m_data_List[[i+1]])
    }
    ## allow for dropped 'bad' vals
    for (i in 1:(length(d_data_List)-1)) {
        expect_equivalent(d_data_List[[i]],d_data_List[[i+1]])
    }

})


test_that("lmerForm", {

    set.seed(101)

    x <- rnorm(10)
    y <- rnorm(10)
    z <- rnorm(10)
    r <- sample(1:3, size=10, replace=TRUE)
    d <- data.frame(x,y,z,r)

    ## example from Joehanes Roeby
    m2 <- suppressWarnings(lmer(x ~ y + z + (1|r), data=d))
    ff <- function() {
        m1 <- suppressWarnings(lmer(x ~ y + z + (1|r), data=d))
        return(anova(m1))
    }

    ff1 <- Reaction ~ Days + (Days|Subject)
    fm1 <- lmer(ff1, sleepstudy)
    fun <- function () {
        ff1 <- Reaction ~ Days + (Days|Subject)
        fm1 <- suppressWarnings(lmer(ff1, sleepstudy))
        return (anova(fm1))
    }
    anova(m2)
    ff()
    expect_equal(anova(m2),ff())
    anova(fm1)
    fun()
    expect_equal(anova(fm1),fun())
})

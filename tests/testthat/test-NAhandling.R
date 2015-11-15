stopifnot(require("testthat"), require("lme4"))

context("NA (and Inf) handling")

## Modified sleepstudy data :
sleepst.a <- sleepstudy
rownames(sleepst.a) <- paste0("a", rownames(sleepstudy))
sleepstudyNA  <- within(sleepst.a,  Reaction[1:3] <- NA)
sleepstudyNA2 <- within(sleepst.a,  Days[1:3] <- NA)
sleepInf      <- within(sleepstudy, Reaction[Reaction > 400] <- Inf)

## Modified  cake  data :
cakeNA <- rbind(cake, tail(cake,1))
cakeNA[nrow(cakeNA), "angle"] <- NA
## Create new data frame with some NAs in fixed effect
cakeNA.X <- within(cake, temp[1:5] <- NA)
## NA values in random effects -- should get treated
cakeNA.Z <- within(cake, replicate[1:5] <- NA)

test_that("naming", {
    ## baseline model
    fm1 <- lmer(Reaction~Days+(Days|Subject), sleepst.a)
    ## default: na.omit
    fm2 <- update(fm1, data=sleepstudyNA,
                  control=lmerControl(check.conv.grad="ignore"))
    expect_equal(head(names(fitted(fm1))), paste0("a",1:6))
    expect_equal(head(names(fitted(fm2))), paste0("a",4:9))
    expect_equal(names(predict(fm2)), names(fitted(fm2)))
    expect_equal(length(p1 <- predict(fm2)), 177)
    ## predict with na.exclude -> has 3 NA's, but otherwise identical:
    expect_equal(length(p2 <- predict(fm2, na.action=na.exclude)), 180)
    expect_identical(p1, p2[!is.na(p2)])
    expect_equal(length((s1 <- simulate(fm1,1))[[1]]),180)
    expect_equal(length((s2 <- simulate(fm2,1))[[1]]),177)
    expect_equal(head(rownames(s1)),paste0("a",1:6))
    expect_equal(head(rownames(s2)),paste0("a",4:9))

    ## test simulation
    expect_is(attr(simulate(fm2),"na.action"),"omit")
    expect_is(refit(fm2,simulate(fm2)),"merMod")
    expect_equal(fixef(fm2),
                 fixef(refit(fm2, sleepstudyNA$Reaction)), tolerance = 1e-5)
    fm2ex <- update(fm2, na.action=na.exclude)
    expect_equal(nrow(ss2 <- simulate(fm2ex)),180)
    expect_is(refit(fm2,ss2[[1]]),"merMod")
    ## issue #197, 18 new subjects; some with NA in y
    d2 <- sleepstudyNA[c(1:180, 1:180),]
    d2[,"Subject"] <- factor(rep(1:36, each=10))
    d2[d2$Subject == 19, "Reaction"] <- NA
    expect_equal(dim( simulate(fm1, newdata=d2, allow.new.levels=TRUE) ), c(360,1))

    ## na.pass (pretty messed up)
    expect_error(update(fm1,data=sleepstudyNA,
                  control=lmerControl(check.conv.grad="ignore"),
                  na.action=na.pass),
                 "NA/NaN/Inf in 'y'")
    sleepstudyNA2 <- within(sleepst.a, Days[1:3] <- NA)
    expect_error(fm4 <- update(fm1, data = sleepstudyNA2,
                               control=lmerControl(check.conv.grad="ignore"),
                               na.action=na.pass),"NA in Z")
    expect_is(suppressWarnings(confint(fm2,method="boot",nsim=3,
                                       quiet=TRUE)),"matrix")
    expect_error(update(fm1, data = sleepstudyNA2,
                        control = lmerControl(check.conv.grad="ignore"),
                        na.action = na.pass),
                 "NA in Z")
    expect_is(suppressWarnings(
                  ci2 <- confint(fm2, method="boot", nsim=3, quiet=TRUE)), "matrix")
})

test_that("other_NA", {
    expect_error(lmer(Reaction ~ Days + (Days | Subject), sleepInf), "\\<Inf\\>")

    fm0 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
    ## NA's in response :
    fm1 <- update(fm0, data = cakeNA)
    expect_true(all.equal(  fixef(fm0),  fixef(fm1)))
    expect_true(all.equal(VarCorr(fm0),VarCorr(fm1)))
    expect_true(all.equal(  ranef(fm0),  ranef(fm1)))

    fm1_omit <-  update(fm1, na.action = na.omit)
    fm1_excl <-  update(fm1, na.action = na.exclude)
    expect_error(update(fm1, na.action = na.pass), "NA/NaN")
    expect_error(update(fm1, na.action = na.fail), "missing values in object")
    fm1_omit@call <- fm1@call ## <- just for comparing:
    expect_equal(fm1, fm1_omit)
    expect_equal(length(fitted(fm1_omit)), 270)
    expect_equal(length(fitted(fm1_excl)), 271)
    expect_true(is.na(tail(predict(fm1_excl),1)))

    ## test predict.lm
    d <- data.frame(x = 1:10, y = c(rnorm(9),NA))
    lm1 <- lm(y~x, data=d, na.action=na.exclude)
    expect_is(predict(lm1), "numeric")
    expect_equal(1, sum(is.na(predict(lm1, newdata = data.frame(x=c(1:4,NA))))))

    ## Triq examples ...
    m.lmer <- lmer (angle ~ temp + (1 | recipe) + (1 | replicate), data=cake)
    ## NAs in fixed effect
    p1_pass <- predict(m.lmer, newdata=cakeNA.X, re.form=NA,
                       na.action=na.pass)
    expect_true(length(p1_pass)==nrow(cakeNA.X))
    expect_true(all(is.na(p1_pass[1:5])))
    p1_omit <- predict(m.lmer, newdata=cakeNA.X, re.form=NA,
                       na.action=na.omit)
    p1_exclude <- predict(m.lmer, newdata=cakeNA.X, re.form=NA,
                          na.action=na.exclude)
    expect_true(length(p1_omit)==nrow(na.omit(cakeNA.X)))
    expect_true(length(p1_exclude)==nrow(cakeNA.X))
    expect_true(all.equal(c(na.omit(p1_exclude)),p1_omit))
    expect_error(predict(m.lmer, newdata=cakeNA.X, re.form=NA, na.action=na.fail),
                 "missing values in object")

    ## now try it with re.form==NULL
    p2_pass <- predict(m.lmer, newdata=cakeNA.X, re.form=NULL,
                       na.action=na.pass)
    expect_true(length(p2_pass)==nrow(cakeNA.X))
    expect_true(all(is.na(p2_pass[1:5])))
    p2_omit <- predict(m.lmer, newdata=cakeNA.X, re.form=NULL,
                       na.action=na.omit)
    p2_exclude <- predict(m.lmer, newdata=cakeNA.X, re.form=NULL,
                          na.action=na.exclude)
    expect_true(length(p2_omit)==nrow(na.omit(cakeNA.X)))
    expect_true(all.equal(c(na.omit(p2_exclude)),p2_omit))
    expect_error(predict(m.lmer, newdata=cakeNA.X, re.form=NULL, na.action=na.fail),
                 "missing values in object")

    ## experiment with NA values in random effects -- should get treated
    expect_error(predict(m.lmer, newdata=cakeNA.Z, re.form=NULL),
                "NAs are not allowed in prediction data")
    p4 <- predict(m.lmer, newdata=cakeNA.Z, re.form=NULL,
                  allow.new.levels=TRUE)
    p4B <- predict(m.lmer, newdata=cakeNA.Z, re.form=~1|recipe,
                   allow.new.levels=TRUE)
    expect_true(all.equal(p4[1:5],p4B[1:5]))
    p4C <- predict(m.lmer, newdata=cakeNA.Z, re.form=NA)

    d <- data.frame(x=runif(100),f=factor(rep(1:10,10)))
    set.seed(101)
    u <- rnorm(10)
    d <- transform(d,y=rnorm(100,1+2*x+u[f],0.2))
    d0 <- d
    d[c(3,5,7),"x"] <- NA

    ## 'omit' and 'exclude' are the only choices under which
    ##  we will see NA values in the results
    fm0 <- lmer(y~x+(1|f), data=d0)
    ## no 'na.action' attribute because no NAs in this data set
    expect_equal(attr(model.frame(fm0),"na.action"),NULL)
    fm1 <- update(fm0, data=d)
    ## no NAs in predict or residuals because na.omit
    expect_false(any(is.na(predict(fm1))))
    expect_false(any(is.na(residuals(fm1))))
    fm2 <- update(fm1,na.action="na.exclude")
    ## no NAs in predict or residuals because na.omit
    nNA <- sum(is.na(d$x))
    expect_equal(sum(is.na(predict(fm2))),nNA)
    expect_equal(sum(is.na(residuals(fm2))),nNA)
    expect_error(fm3 <- lmer(y~x+(1|f), data=d, na.action="na.pass"),
                 "Error in qr.default")
    expect_is(refit(fm0),"merMod")
    expect_is(refit(fm1),"merMod")
    expect_is(refit(fm2),"merMod")
})


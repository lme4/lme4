library(lme4)
library(influence.ME)
library(broom.mixed)
library(Matrix)
library(parallel)

## simulate artificial data
simfun <- function(theta = 0, seed = 101, n = 100, nf = 2, sigma = 5) {
    if (!is.null(seed)) set.seed(seed)
    dd <- data.frame(x = rnorm(n),
                     f = factor(rep(1:nf, each = n/nf)))
    dd$y <- suppressMessages(simulate(~ x + (1|f),
                     newdata = dd,
                     newparams = list(beta = c(2,2),
                                      theta = theta,
                                      sigma = sigma),
                     family = gaussian)[[1]])
    return(dd)
}

## compute all influence measures
all_influence <- function(m, n_correct = TRUE, data, debug = FALSE) {

    if (debug) cat("hat values\n")
    i1 <- lme4:::cooks.distance.merMod(m)
    i2 <- stats::cooks.distance(m)  ## calls lme4::cooks.distance
    i3 <- broom.mixed::augment(m)$.cooksd

    ## code:
    ## p <- Matrix::rankMatrix(getME(model, "X"))
    ## hat <- hatvalues(model)
    ## dispersion <- sigma(model)^2
    ## res <- residuals(model, type = "pearson")
    ## res <- (res/(1 - hat))^2 * hat/(dispersion * p)
    ## res[is.infinite(res)] <- NaN
    ## res

    ## evaluation is tricky, need 'data' arg

    if (debug) cat("brute force (influence.merMod)\n")
    n <- lme4:::influence.merMod(m, data = data)
    i4 <- cooks.distance(n)
    if (debug) cat("brute force (influence.ME)\n")
    p <- try(influence.ME::influence(m, obs = TRUE), silent = TRUE)
    ## fails under parLapply - problem with update/evaluation scope
    if (is(p, "try-error")) {
        i5 <- rep(NA_real_, length(i4))
    } else {
        i5 <- c(cooks.distance(p))
    }
    
    ret <- cbind(lme4_hat=i1, stat_hat=i2, broom_hat=i3,
          lme4_brute=i4, iME_brute = i5)

    ## scaling to make hatvalues results match hatvalues(lm)
    if (n_correct) {
        n <- nobs(m)
        nc <- n - length(fixef(m))
        ret[,1:3] <- ret[,1:3]*nc/n
    }

    return(ret)
}

## scaled root mean square error
srmse <- function(x,y) {
    sqrt(mean(((x-y)/(x+y)/2)^2))
}

## sRMSE comparison of all columns
srmse_mat <- function(x) {
    n <- ncol(x)
    x_srmse <- matrix(NA, nrow = n, ncol =n,
                      dimnames = list(colnames(x), colnames(x)))
    for (i in 1:n) {
        for (j in 1:n) {
            x_srmse[i,j] <- srmse(x[,i], x[,j])
        }
    }
    as(zapsmall(x_srmse, 10), "sparseMatrix")
}

## find discrepancies between brute-force and hat-values
get_discrepmat <- function(mList) {
    discrep <- matrix(NA, ncol = 2, nrow = length(tvec),
                      dimnames = list(NULL, c("lme4_brute", "iME_brute")))
    for (i in seq_along(mList)) {
        mm <- srmse_mat(mList[[i]])
        discrep[i,"lme4_brute"] <- mm["lme4_brute", "lme4_hat"]
        discrep[i, "iME_brute"] <- mm["iME_brute", "lme4_hat"]
    }
    discrep
}

###
m <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)

dd <- simfun()
fm2 <- lmer(y~x + (1|f), dd, REML = FALSE)
fm2L <- lm(y~x , dd)

## influence measures for a range of values
tvec <- seq(0, 2, by = 0.05)
mList <- suppressMessages(
    lapply(tvec,
           function(t) {
               dd <- simfun(theta = t)
               m <- lmer(y~x + (1|f), data = dd)
               all_influence(m, data = dd)
           })
)

## now for a larger simulated data set
cl <- makeCluster(20, outfile = "")
clusterExport(cl, c("simfun", "all_influence"))
    invisible(clusterEvalQ(cl, sapply(FUN = library,
                                      c("broom.mixed", "lme4", "influence.ME",
                                        "Matrix"),
                                      character.only = TRUE)))

mList_big <- suppressMessages(
    parLapply(cl = cl,
              tvec,
              function(t) {
                  cat(t,"\n")
                  dd <- simfun(theta = t, n = 1000)
                  m <- lmer(y~x + (1|f), data = dd)
                  all_influence(m, data = dd)
              })
)
stopCluster(cl)

## summarize
aa <- all_influence(m, data = sleepstudy)
aa2 <- cbind(suppressMessages(all_influence(fm2, data = dd)),
             lm_hat = cooks.distance(fm2L))

discrep <- get_discrepmat(mList)
discrep_big <- get_discrepmat(mList_big)
colnames(discrep_big) <- paste0(colnames(discrep_big), "_n=1000")
discrep_all <- cbind(discrep, discrep_big)[,-4] ## last column is all NA

png("cooks_distance.png")
par(las=1, bty = "n")
cvec <- c(1, 2, 4)
matplot(tvec, discrep_all, type = "b", col = cvec, pch = cvec,
        xlab = "random-effects SD",
        ylab = "scaled RMSE from hatvalues result")
legend("bottomright", col = cvec, lty = cvec, pch = cvec,
       legend = colnames(discrep_all))
dev.off()

pairs(aa2, gap = 0,
      panel = function(x,y, ...) {
          points(x,y)
          abline(a=0, b=1, col = 2)
      })


print(srmse_mat(aa2), digits=2)
print(srmse_mat(aa), digits=2)

MASS::fractions(sigma(fm2L)^2/sigma(fm2)^2)  ## 50/49



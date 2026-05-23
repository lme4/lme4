get_est <- function(mod, model_name) {
  v <- as.data.frame(VarCorr(mod))
  vnames <- with(v,
                 ifelse(is.na(var2),
                        paste0("sd_", var1, "|", grp),
                        paste0("cor_", var1, ".", var2, "|", grp)))
  est <- c(fixef(mod), setNames(v$sdcor, vnames))
  data.frame(model=model_name, var=names(est), est=est, row.names=NULL)
}

ncores <- min(10, parallel::detectCores() - 1)

b_cifun <- function(x, nsim = 501) {
  cat(deparse(getCall(x)$formula), "\n")
  confint(x, method = "boot", seed = 101, nsim = nsim, signames = FALSE,
          parallel = "multicore", ncpus = ncores)
}

p_fun <- function(x) {
  cat(deparse(getCall(x)$formula), "\n")
  profile(x, signames = FALSE,
          devtol = 1e-2,
          maxpts = 250,
          parallel = "multicore", ncpus = ncores)
}


combfun1 <- function(x, model, model_id) {
   df <- data.frame(var = rownames(x), x, check.names = FALSE)
   names(df)[names(df) == "2.5 %"]  <- "lwr"
   names(df)[names(df) == "97.5 %"] <- "upr"
   rownames(df) <- NULL
   df$model <- model
   df$model_id <- model_id
   df
}

combfun2 <- function(x, type, df_name, df_est) {
  df <- do.call("rbind", Map(combfun1, x, df_name$mdesc, df_name$mnames))
  ## if necessary, hack correlation names (ugh); different parts of the pipeline
  ## put the term names in a different order
  ## this assumes that this is the **ONLY** model we are considering with a correlation ...
  corvals <- grep("cor", df$var)
  if (length(corvals) > 0 && !any(grepl("urban", df$var[corvals]))) stop("unexpected cor name pattern in combfun2")
  if (length(corvals) > 0) {
    df$var[corvals] <- "cor_(Intercept).urban1|district"
  }
  df <- merge(df, df_est, by.x = c("model_id", "var"), by.y = c("model", "var"))
  df$model <- factor(df$model, levels = df_name$mdesc)
  df$type <- type
  df <- df[order(df$model),]
  df
}


re_wald_ci_hack <- function(x, q = 0.95) {
  lev <- qnorm((q+1)/2)
  th <- unname(getME(x, "theta"))
  nth <- length(th)
  H <- x@optinfo$derivs$Hessian
  v <- try(solve(H)[1:nth, 1:nth, drop = FALSE], silent = TRUE)
  if (inherits(v, "try-error")) return (data.frame(lwr = rep(NA_real_, nth), upr = rep(NA_real_, nth)))
  ## in this special case, this means one or two
  ##  scalar REs, each of which correspond to a SD
  if (nth <= 2) {
    vals <- th
    s <- sqrt(diag(v))
  } else {
    ## delta method hack.
    phi <- th[2]^2 + th[3]^2
    vals <- c(sd1 = th[1], sd2 = sqrt(phi), cor = th[2]/sqrt(phi))
    s <- rep(NA_real_, 3)
    s[1] <- sqrt(v[1,1])
    cvec2 <- c(0, th[2]/sqrt(phi), th[3]/sqrt(phi))
    s[2] <- sqrt(drop(t(cvec2) %*% (v %*% cvec2)))
    cvec3 <- c(0, (1-th[2]^2/phi)/sqrt(phi), -th[2]*th[3]*phi^(-3/2))
    s[3] <- sqrt(drop(t(cvec3) %*% (v %*% cvec3)))
  }
  return(data.frame(lwr = vals - lev*s, upr = vals + lev*s))
}

wald_cifun <- function(x) {
  cc <- confint(x, method = "Wald", signames = FALSE)
  nth <- length(getME(x, "theta"))
  cc[1:nth, ] <- as.matrix(re_wald_ci_hack(x))
  cc
}

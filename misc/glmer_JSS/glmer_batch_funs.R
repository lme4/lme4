get_est <- function(mod, model_name) {
  v <- as.data.frame(VarCorr(mod))
  vnames <- with(v,
                 ifelse(is.na(var2),
                        paste0("sd_", var1, "|", grp),
                        paste0("cor_", var1, ".", var2, "|", grp)))
  est <- c(fixef(mod), setNames(v$sdcor, vnames))
  data.frame(model=model_name, var=names(est), est=est, row.names=NULL)
}

b_cifun <- function(x) {
  cat(deparse(getCall(x)$formula), "\n")
  confint(x, method = "boot", seed = 101, nsim = 501, signames = FALSE,
          parallel = "multicore", ncpus = 10)
}

p_fun <- function(x) {
  cat(deparse(getCall(x)$formula), "\n")
  profile(x, signames = FALSE,
          devtol = 1e-2,
          maxpts = 250,
          parallel = "multicore", ncpus = 10)
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
  df <- merge(df, df_est, by.x = c("model_id", "var"), by.y = c("model", "var"))
  df$model <- factor(df$model, levels = df_name$mdesc)
  df$type <- type
  df <- df[order(df$model),]
  df
}

library("lme4")

cbpp2 <- read.csv("cbpp2.csv")
cbpp2 <- transform(
  cbpp2,
  period = factor(period),
  treatment = factor(
    treatment,
    levels = c("Partial/null", "Complete", "Unknown")
  ),
  obs = factor(seq(nrow(cbpp2))),
  avg_size = avg_size / sd(avg_size, na.rm = TRUE)
)

gm_herd <- glmer(
  incidence / size ~ period + treatment + avg_size + (1 | herd),
  family = binomial,
  data = cbpp2,
  weights = size,
  control = glmerControl(optimizer = "bobyqa")
)

gm_herdobs <- update(gm_herd, . ~ . + (1 | obs)) ## herd and observation-level REs
gm_obs <- update(gm_herd, . ~ . - (1 | herd) + (1 | obs)) ## observation-level REs only

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

mnames <- c("herd", "obs", "herdobs")
mod_list <- mget(sprintf("gm_%s", mnames))
cbpp_confint_boot <- lapply(mod_list, b_cifun)
names(cbpp_confint_boot) <- mnames

cbpp_prof <- lapply(mod_list, p_fun)

if (FALSE) {
  ## experimenting/exploring
  profile.gm1 <- cbpp_confint_prof$herd
  lattice::xyplot(profile.gm1)

  library(ggplot2)
  p2 <- as.data.frame(profile.gm1)
  ggplot(p2, aes(.focal, .zeta)) +
    geom_point() +
    geom_line() +
    facet_wrap(~.par, scale = "free")
}

cbpp_confint_prof <- lapply(cbpp_prof, confint)

cbpp_confint_wald <- lapply(mod_list, function(x) confint(x, method = "Wald", signames = FALSE))

get_est <- function(mod, model_name) {
  v <- as.data.frame(VarCorr(mod))
  vnames <- paste0("sd_", v$var1, "|", v$grp)
  est <- c(fixef(mod), setNames(v$sdcor, vnames))
  data.frame(model=model_name, var=names(est), est=est, row.names=NULL)
}

cbpp_est <- do.call("rbind",
                     Map(get_est, mod_list, mnames))
rownames(cbpp_est) <- NULL

combfun <- function(
save(
  list = c("cbpp_confint_prof", "cbpp_confint_boot", "cbpp_confint_wald", "cbpp_prof", "cbpp_est"),
  file = "CBPP_batch.rda"
)
sessionInfo()

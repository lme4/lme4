
library(glmmTMB)

# Model 1; taken from github.com/lme4/lme4/issues/962
set.seed(101)
form <- out_bin_1 ~ ftime + trt + ar1(ftime + 0 | id_cluster) +
  (1 | id_individual)
sfun <- function(n) factor(sample(1:n, replace=TRUE, size = 1465))
dd <- data.frame(id_cluster = sfun(7), id_individual = sfun(486),
                 ftime = sfun(4), trt = sfun(2))
np <- list(beta = c(0, rep(0.1, 3), 1), theta = rep(0, 3))
dd$out_bin_1 <- simulate_new(form[-2], newdata = dd,
                             newparams = np,
                             family = binomial)[[1]]

dd$out_gauss <- simulate_new(form[-2], newdata = dd,
                             newparams = c(np, list(betadisp=0)),
                             family = gaussian)[[1]]

glmmTMB_1 <- glmmTMB(form, data = dd, family = binomial)

glmmTMB_glmer_ar1_comparison <- list(
  glmer_ar1_dd = dd,
  glmmTMB_ar1_binomial = glmmTMB_1
)

saveRDS(glmmTMB_glmer_ar1_comparison, "inst/testdata/glmmTMB_glmer_ar1_comparison.rds")

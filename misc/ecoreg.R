download.file("https://github.com/bbolker/mmd_utils/raw/refs/heads/master/ecoreg.rds", dest = "ecoreg.rds")
ee <- readRDS("ecoreg.rds")

form_full <- mmamm_log ~
    (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc)^2 +
    (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc | biome)

form_diag <- mmamm_log ~
  (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc)^2 +
  diag(NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc | biome)

form_cs <- mmamm_log ~
  (NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc)^2 +
  cs(NPP_log_sc + Feat_log_sc + NPP_cv_sc + Feat_cv_sc | biome)

form_list <- namedList(form_full, form_diag, form_cs)
lmer_list <- lapply(form_list, \(f) try(lmer(f, ee, REML = FALSE)))

library(glmmTMB)
glmmTMB_list <- lapply(form_list, \(f) try(glmmTMB(f, ee, REML = FALSE)))

lapply(glmmTMB_list, VarCorr)

lapply(lmer_list, \(m) try(VarCorr(m)))
try(VarCorr(lmer_list[[1]]))  ## this is broken and wasn't before?


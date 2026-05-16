## Helper functions for glmmTMB batch fitting, parallel to glmer_batch_funs.R.
## Key differences from lme4:
##   - fixef(mod) returns a list; use fixef(mod)$cond for conditional fixed effects
##   - VarCorr(mod) returns a list; use VarCorr(mod)$cond for conditional RE
##   - confint() row names have a "cond." prefix for fixed effects — normalize_tmb_names() strips it
##   - Wald CIs for RE parameters are labeled "Std.Dev.X|grp" — normalize_tmb_names() renames them
##   - profile/uniroot RE rows ("theta_i|grp.j") are on the internal log(sd) scale —
##     back-transformed to sdcor scale by backtransform_theta() for scalar RE groups;
##     multi-RE groups (Cholesky off-diagonals) are not back-transformed and will be absent from plots

get_est_tmb <- function(mod, model_name) {
  vc <- VarCorr(mod)$cond
  rows <- lapply(names(vc), function(grp) {
    m   <- vc[[grp]]
    sds <- attr(m, "stddev")
    nv  <- length(sds)
    out <- data.frame(grp = grp, var1 = names(sds), var2 = NA_character_,
                      sdcor = unname(sds), stringsAsFactors = FALSE)
    if (nv > 1) {
      co  <- attr(m, "correlation")
      idx <- which(upper.tri(co), arr.ind = TRUE)
      out <- rbind(out, data.frame(
        grp   = grp,
        var1  = rownames(co)[idx[, 1]],
        var2  = colnames(co)[idx[, 2]],
        sdcor = co[idx],
        stringsAsFactors = FALSE
      ))
    }
    out
  })
  v      <- do.call(rbind, rows)
  vnames <- with(v, ifelse(is.na(var2),
                           paste0("sd_", var1, "|", grp),
                           paste0("cor_", var1, ".", var2, "|", grp)))
  est <- c(fixef(mod)$cond, setNames(v$sdcor, vnames))
  data.frame(model = model_name, var = names(est), est = est, row.names = NULL)
}

## Build a mapping from glmmTMB theta parameter names ("theta_i|grp.j") to
## sdcor names ("sd_X|grp"), for scalar RE groups only (nv == 1).
## glmmTMB numbers thetas globally across groups in formula order; each scalar
## group consumes one theta. Multi-RE groups (nv > 1) are skipped — their
## Cholesky off-diagonal elements are not simply log(sd) and cannot be
## back-transformed without full model-structure awareness.
build_theta_map <- function(mod) {
  vc         <- VarCorr(mod)$cond
  result     <- character(0)
  global_idx <- 1L
  for (grp in names(vc)) {
    m       <- vc[[grp]]
    sds     <- attr(m, "stddev")
    nv      <- length(sds)
    n_theta <- nv * (nv + 1L) / 2L
    if (nv == 1L) {
      theta_name          <- paste0("theta_", global_idx, "|", grp, ".1")
      result[theta_name]  <- paste0("sd_", names(sds)[1], "|", grp)
    }
    global_idx <- global_idx + n_theta
  }
  result
}

## Back-transform scalar-RE theta rows from log(sd) to sd scale using theta_map.
## Only CI bound columns ("2.5 %" / "97.5 %") are exponentiated; other columns
## (e.g. "Estimate" in Wald output) are left unchanged.
backtransform_theta <- function(ci, theta_map) {
  rn      <- rownames(ci)
  ci_cols <- colnames(ci) %in% c("2.5 %", "97.5 %")
  for (i in seq_len(nrow(ci))) {
    if (rn[i] %in% names(theta_map)) {
      ci[i, ci_cols] <- exp(ci[i, ci_cols])
      rn[i]          <- theta_map[rn[i]]
    }
  }
  rownames(ci) <- rn
  ci
}

## Normalize glmmTMB confint row names to match get_est_tmb conventions:
##   - strip "cond." prefix from fixed effects
##   - rename Wald RE rows: "Std.Dev.X|grp" -> "sd_X|grp", "Cor.X.Y|grp" -> "cor_X.Y|grp"
## Called after backtransform_theta(), so scalar-RE theta rows have already been
## renamed to "sd_X|grp"; any remaining "theta_i|grp.j" rows (multi-RE off-diagonals)
## are left as-is and will be dropped by the merge in combfun2.
normalize_tmb_names <- function(x) {
  rn <- rownames(x)
  rn <- sub("^cond\\.", "", rn)
  rn <- sub("^Std\\.Dev\\.", "sd_", rn)
  rn <- sub("^Cor\\.",      "cor_", rn)
  rownames(x) <- rn
  x
}

## glmmTMB does not support method = "boot" in confint(); uniroot is the
## closest substitute — it finds CIs by root-finding on the profile likelihood.
## Run sequentially (no parallel=) so tryCatch can catch errors from child calls.
u_cifun_tmb <- function(x) {
  cat(deparse(getCall(x)$formula), "\n")
  tryCatch({
    ci        <- confint(x, method = "uniroot")
    theta_map <- build_theta_map(x)
    ci        <- backtransform_theta(ci, theta_map)
    normalize_tmb_names(ci)
  }, error = function(e) {
    warning("uniroot CI failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}

## Note: glmmTMB profile CIs are on the internal (log/Cholesky) scale for RE
## parameters, unlike lme4's sdcor scale — direct RE comparison requires
## back-transformation. Fixed-effect CIs are on the same scale as lme4.
## Returns NULL (with a warning) when profiling fails (e.g. near-zero RE variance).
p_cifun_tmb <- function(x) {
  cat(deparse(getCall(x)$formula), "\n")
  tryCatch({
    ci        <- confint(x, method = "profile")
    theta_map <- build_theta_map(x)
    ci        <- backtransform_theta(ci, theta_map)
    normalize_tmb_names(ci)
  }, error = function(e) {
    ## Full profile failed (typically a near-singular RE variance parameter).
    ## Fall back to profiling fixed effects only, which are numerically stable.
    ## No theta back-transformation needed since beta_ rows have no theta names.
    message("Full profile CI failed (", conditionMessage(e),
            "); falling back to fixed-effects-only profile.")
    tryCatch({
      ci <- confint(x, method = "profile", parm = "beta_")
      normalize_tmb_names(ci)
    }, error = function(e2) {
      warning("Fixed-effects profile CI also failed: ", conditionMessage(e2), call. = FALSE)
      NULL
    })
  })
}

w_cifun_tmb <- function(x) {
  tryCatch({
    normalize_tmb_names(confint(x, method = "Wald"))
  }, error = function(e) {
    warning("Wald CI failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}

combfun1 <- function(x, model, model_id) {
  df <- data.frame(var = rownames(x), x, check.names = FALSE)
  names(df)[names(df) == "2.5 %"]  <- "lwr"
  names(df)[names(df) == "97.5 %"] <- "upr"
  rownames(df) <- NULL
  df$model <- model
  df$model_id <- model_id
  df[, c("var", "lwr", "upr", "model", "model_id")]
}

combfun2 <- function(x, type, df_name, df_est) {
  ok <- !sapply(x, is.null)
  if (!any(ok)) return(NULL)
  df <- do.call("rbind", Map(combfun1, x[ok], df_name$mdesc[ok], df_name$mnames[ok]))
  df <- merge(df, df_est, by.x = c("model_id", "var"), by.y = c("model", "var"))
  df$model <- factor(df$model, levels = df_name$mdesc)
  df$type <- type
  df <- df[order(df$model),]
  df
}

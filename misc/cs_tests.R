# VarCorr Workflow Test
library(glmmTMB)
library(lme4)
data("sleepstudy", package = "lme4") ## redundant, but ...

# Models
fit_us <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)
fit_cs_glmm <- glmmTMB(Reaction ~ Days + cs(Days | Subject), sleepstudy, REML = FALSE)
fit_cs <- lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy, REML = FALSE)
fit_cs_hom <- lmer(Reaction ~ Days + cs(Days | Subject, hom = TRUE), sleepstudy, REML = FALSE)
fit_cs_glmm_hom <- glmmTMB(Reaction ~ Days + homcs(Days | Subject), sleepstudy, REML = FALSE)

##############################

# Manual mkVarCorr_for_structure Reconstruction Test


# STEP 1: GET MODEL COMPONENTS
fit_cs <- lmer(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)

# Set the model to test
model <- fit_cs

structure_info <- lme4:::extract_structure_info(model)
cs_structure <- structure_info$structures[[1]]
theta_values <- getME(model, "theta")
sigma_value <- sigma(model)
term_cnms <- getME(model, "cnms")[[1]]

print("Input values:")
print(cs_structure)
print(theta_values)
print(sigma_value)
print(term_cnms)

# STEP 2: MANUAL STRUCTURE OBJECT CREATION
manual_structure <- new("HeterogeneousCSCovariance", dimension = 2L)
print(manual_structure)

# STEP 3: MANUAL PARAMETER ASSIGNMENT
d <- 2L
n_v_params <- d

if (n_v_params <= length(theta_values)) {
    manual_vparams <- theta_values[seq_len(n_v_params)]

    if (length(theta_values) > n_v_params) {
        manual_cparams <- theta_values[-seq_len(n_v_params)]
    } else {
        manual_cparams <- numeric(0)
    }
} else {
    manual_vparams <- theta_values
    manual_cparams <- numeric(0)
}

print(manual_vparams)
print(manual_cparams)

manual_structure@vparameters <- manual_vparams
manual_structure@cparameters <- manual_cparams

comparison_structure <- new("HeterogeneousCSCovariance", dimension = 2L)
comparison_structure <- lme4:::set_parameters(comparison_structure, theta_values)

print(all.equal(manual_structure@vparameters, comparison_structure@vparameters) && 
      all.equal(manual_structure@cparameters, comparison_structure@cparameters))

# STEP 4: MANUAL CORRELATION TRANSFORMATION
theta_corr <- manual_cparams[1]
print(theta_corr)

a <- 1/(d-1)
manual_rho <- plogis(theta_corr) * (1 + a) - a
print(manual_rho)

actual_rho <- lme4:::cs_theta_to_rho(theta_corr, d)
print(actual_rho)
print(all.equal(manual_rho, actual_rho))

target_glmmTMB_rho <- 0.081
required_theta_for_glmmTMB <- lme4:::cs_rho_to_theta(target_glmmTMB_rho, d)
print(required_theta_for_glmmTMB)
print(manual_rho)
print(theta_corr - required_theta_for_glmmTMB)
print(required_theta_for_glmmTMB/theta_corr)

# STEP 5: MANUAL COVARIANCE MATRIX CONSTRUCTION
manual_R <- matrix(manual_rho, nrow = d, ncol = d)
diag(manual_R) <- 1.0
print(manual_R)

manual_D <- diag(manual_vparams)
print(manual_D)

manual_Sigma <- manual_D %*% manual_R %*% manual_D
print(manual_Sigma)
print(diag(manual_Sigma))

actual_cov_matrix <- lme4:::compute_covariance_matrix(manual_structure)
print(manual_Sigma)
print(actual_cov_matrix)
print(all.equal(as.matrix(manual_Sigma), as.matrix(actual_cov_matrix)))

# STEP 6: MANUAL CORRELATION MATRIX EXTRACTION
manual_stddevs_from_cov <- sqrt(diag(manual_Sigma))
print(manual_stddevs_from_cov)

manual_D_inv <- diag(1 / manual_stddevs_from_cov)
manual_R_extracted <- manual_D_inv %*% manual_Sigma %*% manual_D_inv
print(manual_R_extracted)

actual_cor_matrix <- lme4:::compute_correlation_matrix(manual_structure)
print(manual_R_extracted)
print(actual_cor_matrix)
print(all.equal(as.matrix(manual_R_extracted), as.matrix(actual_cor_matrix)))

# STEP 7: MANUAL SIGMA SCALING AND FINAL RESULT CONSTRUCTION
print(sigma_value)
print(sigma_value^2)
manual_scaled_cov <- sigma_value^2 * manual_Sigma
print(manual_scaled_cov)

manual_final_stddevs <- sqrt(diag(manual_scaled_cov))
print(manual_final_stddevs)

dimnames(manual_scaled_cov) <- list(term_cnms, term_cnms)
dimnames(manual_R_extracted) <- list(term_cnms, term_cnms)
names(manual_final_stddevs) <- term_cnms

manual_final_result <- structure(manual_scaled_cov,
                                stddev = manual_final_stddevs,
                                correlation = as.matrix(manual_R_extracted))

print(manual_final_result)
print(attr(manual_final_result, "stddev"))
print(attr(manual_final_result, "correlation"))

# STEP 8: THETA SPLITTING AND PARAMETER SIZE VERIFICATION
param_sizes <- structure_info$param_sizes
print(param_sizes)

theta_splits <- lme4:::split_theta_by_structure(theta_values, param_sizes)
actual_theta_block <- theta_splits[[1]]
print(actual_theta_block)
print(theta_values)
print(all.equal(actual_theta_block, theta_values))

manual_test_structure <- new("HeterogeneousCSCovariance", dimension = 2L)
manual_test_structure <- lme4:::set_parameters(manual_test_structure, theta_values)
print(manual_test_structure@vparameters)
print(manual_test_structure@cparameters)

actual_test_structure <- new("HeterogeneousCSCovariance", dimension = 2L)
actual_test_structure <- lme4:::set_parameters(actual_test_structure, actual_theta_block)
print(actual_test_structure@vparameters)
print(actual_test_structure@cparameters)
print(all.equal(manual_test_structure@vparameters, actual_test_structure@vparameters) &&
      all.equal(manual_test_structure@cparameters, actual_test_structure@cparameters))

# STEP 9: COMPARISON
method_result <- lme4:::mkVarCorr_for_structure(cs_structure, actual_theta_block, sigma_value, term_cnms)
print(attr(manual_final_result, "stddev"))
print(attr(method_result, "stddev"))
print(all.equal(attr(manual_final_result, "stddev"), attr(method_result, "stddev")))

print(attr(manual_final_result, "correlation"))
print(attr(method_result, "correlation"))
print(all.equal(attr(manual_final_result, "correlation"), attr(method_result, "correlation")))

full_varcorr_result <- VarCorr(model)
subject_component <- full_varcorr_result$Subject
print(attr(manual_final_result, "stddev"))
print(attr(subject_component, "stddev"))
print(all.equal(attr(manual_final_result, "stddev"), attr(subject_component, "stddev")))

print(attr(manual_final_result, "correlation"))
print(attr(subject_component, "correlation"))
print(all.equal(attr(manual_final_result, "correlation"), attr(subject_component, "correlation")))

# STEP 10: PARAMETER FLOW DIAGNOSTIC SUMMARY
print(theta_values)
print(manual_vparams)
print(manual_cparams)
print(manual_cparams[1])
print(manual_rho)
print(diag(manual_Sigma))
print(sigma_value^2)
print(diag(manual_scaled_cov))
print(manual_final_stddevs)
print(manual_R_extracted[1, 2])

###########################################


# Test mkVarCorr_for_structure Pipeline


# Set the model to test
model <- fit_cs

# STEP 1: EXTRACT COMPONENTS
structure_info <- lme4:::extract_structure_info(model)
cs_structure <- structure_info$structures[[1]]
theta_values <- getME(model, "theta")
sigma_value <- sigma(model)
cnms_values <- getME(model, "cnms")

print(cs_structure)
print(theta_values)
print(sigma_value)
print(cnms_values[[1]])

# STEP 2: TEST mkVarCorr_for_structure
theta_block <- theta_values  
term_cnms <- cnms_values[[1]]

print(theta_block)
print(sigma_value)
print(term_cnms)

varcorr_result <- lme4:::mkVarCorr_for_structure(cs_structure, theta_block, sigma_value, term_cnms)
print(varcorr_result)
print(attr(varcorr_result, "stddev"))
print(attr(varcorr_result, "correlation"))

# STEP 3: TEST mkVarCorr
nc <- lengths(cnms_values)
nms <- names(getME(model, "flist"))[attr(getME(model, "flist"), "assign")]

print(sigma_value)
print(nc)
print(theta_values)
print(nms)

mkvarcorr_result <- mkVarCorr(sc = sigma_value,
                              cnms = cnms_values,
                              nc = nc,
                              theta = theta_values,
                              nms = nms,
                              structure_info = structure_info)

print(mkvarcorr_result)

# STEP 4: TEST VarCorr.merMod
varcorr_full_result <- VarCorr(model)
print(varcorr_full_result)

# STEP 5: EXTRACT AND COMPARE VALUES
extract_values <- function(vc_result) {
  if (is.list(vc_result) && length(vc_result) > 0) {
    first_term <- vc_result[[1]]
    stddev <- attr(first_term, "stddev")
    corr <- attr(first_term, "correlation")
    corr_val <- if (is.matrix(corr) && nrow(corr) >= 2) corr[2,1] else NA
  } else {
    stddev <- attr(vc_result, "stddev")
    corr <- attr(vc_result, "correlation")
    corr_val <- if (is.matrix(corr) && nrow(corr) >= 2) corr[2,1] else NA
  }
  
  print(stddev)
  print(corr_val)
  return(list(stddev = stddev, correlation = corr_val))
}

result1 <- extract_values(varcorr_result)
result2 <- extract_values(mkvarcorr_result)
result3 <- extract_values(varcorr_full_result)

# STEP 6: COMPARE WITH glmmTMB
glmm_varcorr <- VarCorr(fit_cs_glmm)
print(glmm_varcorr)

# STEP 7: OPTIMIZER COMPARISON
glmm_theta <- getME(fit_cs_glmm, "theta")
print(glmm_theta)
print(logLik(fit_cs_glmm))

lme4_theta_default <- getME(model, "theta")
print(lme4_theta_default)
print(logLik(model))

# STEP 8: USE glmmTMB fitted theta values as start values for lme4
lmod <- lFormula(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)
original_start <- getME(model, "theta")
print(original_start)

devfun <- do.call(mkLmerDevfun, lmod)
dev_default <- devfun(original_start)
dev_glmm <- devfun(glmm_theta)
print(dev_glmm)

opt_custom <- optimizeLmer(devfun, start = glmm_theta)
print(opt_custom$par)
print(opt_custom$fval)
 
fit <- mkMerMod(environment(devfun), opt_custom, lmod$reTrms, fr = lmod$fr, mc = match.call())
print(logLik(fit))
print(VarCorr(fit))
 
# STEP 9: ALTERNATIVE OPTIMIZER
fit_alt <- lmer(Reaction ~ Days + cs(Days | Subject, hom = FALSE), 
                sleepstudy, REML = FALSE,
                start = glmm_theta,  
                control = lmerControl(optimizer = "bobyqa"))
print(fit_alt)

















# Test mkVarCorr_for_structure Pipeline
# ====================================

library(lme4)
library(glmmTMB)
data(sleepstudy)

# Fit models
fit_cs <- lmer(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)
fit_cs_glmm <- glmmTMB(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)

# Set the model to test
model <- fit_cs

# STEP 1: EXTRACT COMPONENTS
structure_info <- lme4:::extract_structure_info(model)
cs_structure <- structure_info$structures[[1]]
theta_values <- getME(model, "theta")
sigma_value <- sigma(model)
cnms_values <- getME(model, "cnms")

print(cs_structure)
print(theta_values)
print(sigma_value)
print(cnms_values[[1]])

# STEP 2: TEST mkVarCorr_for_structure
theta_block <- theta_values  
term_cnms <- cnms_values[[1]]

print(theta_block)
print(sigma_value)
print(term_cnms)

varcorr_result <- lme4:::mkVarCorr_for_structure(cs_structure, theta_block, sigma_value, term_cnms)
print(varcorr_result)
print(attr(varcorr_result, "stddev"))
print(attr(varcorr_result, "correlation"))

# STEP 3: TEST mkVarCorr
nc <- lengths(cnms_values)
nms <- names(getME(model, "flist"))[attr(getME(model, "flist"), "assign")]

print(sigma_value)
print(nc)
print(theta_values)
print(nms)

mkvarcorr_result <- mkVarCorr(sc = sigma_value,
                              cnms = cnms_values,
                              nc = nc,
                              theta = theta_values,
                              nms = nms,
                              structure_info = structure_info)

print(mkvarcorr_result)

# STEP 4: TEST VarCorr.merMod
varcorr_full_result <- VarCorr(model)
print(varcorr_full_result)

# STEP 5: EXTRACT AND COMPARE VALUES
extract_values <- function(vc_result) {
  if (is.list(vc_result) && length(vc_result) > 0) {
    first_term <- vc_result[[1]]
    stddev <- attr(first_term, "stddev")
    corr <- attr(first_term, "correlation")
    corr_val <- if (is.matrix(corr) && nrow(corr) >= 2) corr[2,1] else NA
  } else {
    stddev <- attr(vc_result, "stddev")
    corr <- attr(vc_result, "correlation")
    corr_val <- if (is.matrix(corr) && nrow(corr) >= 2) corr[2,1] else NA
  }
  
  print(stddev)
  print(corr_val)
  return(list(stddev = stddev, correlation = corr_val))
}

result1 <- extract_values(varcorr_result)
result2 <- extract_values(mkvarcorr_result)
result3 <- extract_values(varcorr_full_result)

# STEP 6: COMPARE WITH glmmTMB
glmm_varcorr <- VarCorr(fit_cs_glmm)
print(glmm_varcorr)

# STEP 7: OPTIMIZER COMPARISON
glmm_theta <- getME(fit_cs_glmm, "theta")
print(glmm_theta)
print(logLik(fit_cs_glmm))

lme4_theta_default <- getME(model, "theta")
print(lme4_theta_default)
print(logLik(model))

# STEP 8: CUSTOM OPTIMIZATION TEST
lmod <- lFormula(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)
original_start <- getME(model, "theta")
print(original_start)

devfun <- do.call(mkLmerDevfun, lmod)
dev_default <- devfun(original_start)
dev_glmm <- devfun(glmm_theta)
print(dev_glmm)

opt_custom <- optimizeLmer(devfun, start = glmm_theta)
print(opt_custom$par)
print(opt_custom$fval)
 
fit_custom <- mkMerMod(environment(devfun), opt_custom, lmod$reTrms, fr = lmod$fr, mc = match.call())
print(logLik(fit_custom))
print(VarCorr(fit_custom))
 
# STEP 9: ALTERNATIVE OPTIMIZER
fit_alt <- lmer(Reaction ~ Days + cs(Days | Subject, hom = FALSE), 
                sleepstudy, REML = FALSE,
                start = glmm_theta,  
                control = lmerControl(optimizer = "bobyqa"))
print(fit_alt)




# Test mkVarCorr_for_structure Pipeline

# Extract key components
structure_info <- lme4:::extract_structure_info(fit_cs)
cs_structure <- structure_info$structures[[1]]
theta_values <- getME(fit_cs, "theta")
sigma_value <- sigma(fit_cs)
cnms_values <- getME(fit_cs, "cnms")


print(cs_structure)
cat("Theta values:", theta_values, "\n")
cat("Sigma value:", sigma_value, "\n")
cat("Term names:", cnms_values[[1]], "\n\n")

# test mkVarCorr_for_structure in isolation

# Get the theta block for the first term
theta_block <- theta_values  
term_cnms <- cnms_values[[1]]

cat("theta_block:", theta_block, "\n")
cat("sc:", sigma_value, "\n")
cat("term_cnms:", term_cnms, "\n\n")

# Call method directly
varcorr_result <- lme4:::mkVarCorr_for_structure(cs_structure, theta_block, sigma_value, term_cnms)
print(varcorr_result)
cat("  Standard deviations:", attr(varcorr_result, "stddev"), "\n")
print(attr(varcorr_result, "correlation"))

# Test full mkVarCorr pipeline


# Extract all components needed for mkVarCorr
nc <- lengths(cnms_values)
nms <- names(getME(fit_cs, "flist"))[attr(getME(fit_cs, "flist"), "assign")]

cat("sc:", sigma_value, "\n")
cat("nc:", nc, "\n") 
cat("theta:", theta_values, "\n")
cat("nms:", nms, "\n")

# Call mkVarCorr directly
mkvarcorr_result <- mkVarCorr(sc = sigma_value,
                              cnms = cnms_values,
                              nc = nc,
                              theta = theta_values,
                              nms = nms,
                              structure_info = structure_info)


print(mkvarcorr_result)

varcorr_full_result <- VarCorr(fit_cs)
print(varcorr_full_result)

# Extract std devs and correlations from each method
extract_values <- function(vc_result, name) {
  if (is.list(vc_result) && length(vc_result) > 0) {
    first_term <- vc_result[[1]]
    stddev <- attr(first_term, "stddev")
    corr <- attr(first_term, "correlation")
    if (is.matrix(corr) && nrow(corr) >= 2) {
      corr_val <- corr[2,1]
    } else {
      corr_val <- NA
    }
  } else {
    stddev <- attr(vc_result, "stddev")
    corr <- attr(vc_result, "correlation")
    if (is.matrix(corr) && nrow(corr) >= 2) {
      corr_val <- corr[2,1]
    } else {
      corr_val <- NA
    }
  }
  
  cat(name, ":\n")
  cat("  Std devs:", if(!is.null(stddev)) stddev else "NULL", "\n")
  cat("  Correlation:", if(!is.na(corr_val)) corr_val else "NA", "\n\n")
  
  return(list(stddev = stddev, correlation = corr_val))
}


result1 <- extract_values(varcorr_result, "mkVarCorr_for_structure")
result2 <- extract_values(mkvarcorr_result, "mkVarCorr")
result3 <- extract_values(varcorr_full_result, "VarCorr.merMod")



# Compare with glmmTMB

fit_cs_glmm <- glmmTMB(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)
glmm_varcorr <- VarCorr(fit_cs_glmm)


print(glmm_varcorr)

####################

# Manual lme4 Optimizer Test with glmmTMB Starting Values
fit_cs_glmm <- glmmTMB(Reaction ~ Days + cs(Days | Subject), sleepstudy, REML = FALSE)
glmm_theta <- getME(fit_cs_glmm, "theta")
print("glmmTMB theta:")
print(glmm_theta)
print("glmmTMB logLik:")
print(logLik(fit_cs_glmm))

# Run lme4 with default starting values 
fit_cs_default <- lmer(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)
lme4_theta_default <- getME(fit_cs_default, "theta")
print("lme4 default theta:")
print(lme4_theta_default)
print("lme4 default logLik:")
print(logLik(fit_cs_default))

# parse formula and create model structure
lmod <- lFormula(Reaction ~ Days + cs(Days | Subject, hom = FALSE), sleepstudy, REML = FALSE)
original_start <- getME(fit_cs_default, "theta")
print(original_start)
devfun <- do.call(mkLmerDevfun, lmod)
env_info <- ls(environment(devfun))
dev_default <- devfun(original_start)

# Try glmmTMB starting values 
dev_glmm <- devfun(glmm_theta)
print("Deviance at glmmTMB start:")
print(dev_glmm)

opt_custom <- optimizeLmer(devfun, start = glmm_theta)
print("Final theta:")
print(opt_custom$par)
print("Final deviance:")
print(opt_custom$fval)
 
fit_custom <- mkMerMod(environment(devfun), opt_custom, lmod$reTrms, fr = lmod$fr, mc = match.call())
 
print("Custom fit logLik:")
print(logLik(fit_custom))
print("Custom fit VarCorr:")
print(VarCorr(fit_custom))
 
# try alternative optimizer
fit_alt <- lmer(Reaction ~ Days + cs(Days | Subject, hom = FALSE), 
                 sleepstudy, REML = FALSE,
                 start = glmm_theta,  
                 control = lmerControl(optimizer = "bobyqa"))
print(fit_alt)


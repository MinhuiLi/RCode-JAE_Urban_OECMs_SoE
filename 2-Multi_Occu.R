# ==============================================================================
# Multi‑variable Occupancy Model for Scale of Effect Analysis
# Script: 2_MultiVar_Occu.R
#
# Description:
#   Fit multi‑covariate occupancy models using pre‑identified optimal‑scale predictors
#
# Workflow dependency:
#   RUN AFTER Script 1 (1_Identify_SoE.R).
#   var_names are optimal‑scale environmental variables extracted from <spe>_EOS_50_2000.xlsx
#
# Usage:
#   !! RUN ONE SPECIES AT A TIME !!
#   Set all parameters in CONFIG section before execution.
#
# Input:
#   - Species‑specific excel table containing optimal‑scale env covariates
#   - External custom R functions (WindSel.R, SiteCovSel.R etc.)
#
# Output:
#  1. <spe>_multi_variable_result.txt     # Model selection report
#  2. <spe>_output.RData                  # SiteCovSel output object
#  3. <spe>_aic_avg_model.RData           # Model‑averaged occupancy model
#  4. <spe>_aic_avg_ci.txt                # 95% CI for model‑averaged coefficients
#  5. <spe>_psi_Env_aicavg.csv            # Predicted occupancy psi
#  6. <spe>_p_aicavg.csv                  # Predicted detection probability p
# ==============================================================================

# ------------------------------------------------------------------------------
# CONFIGURATION: ALL USER‑MODIFIABLE PARAMETERS HERE
# !! CHANGE THESE BEFORE RUNNING FOR A NEW SPECIES !!
# ------------------------------------------------------------------------------
## File paths (consistent with script 1, use / on all OS)
raw_data_path  <- "."    # e.g., "./data/raw"
result_data_path <- "."  # e.g., "./results/occupancy"
fun_script_path <- "." # Folder storing custom WindSel, SiteCovSel etc.

## Species & input file
spe <- "CE"
input_excel_file <- "CE_Pallas_s_squirrel_env_detection_data.xlsx"

## Analysis parameters
n_observations <- 20
collinear_threshold <- 0.7

# Optimal‑scale covariates: extracted from Script1 EOS output <spe>_EOS_50_2000.xlsx
var_names <- c('CH300',
               'BD50',
               'MPS250',
               'PD600',
               'RD100',
               'WSR300',
               'NP350',
               'ED50',
               'ISR250',
               'WD1250',
               'NDVI300',
               'CO700')

# ------------------------------------------------------------------------------
## Step 1: Environment Setup & Data Import (Shared for Both Models)
## Core Goal: Load packages/scripts, set paths, import raw data
# ------------------------------------------------------------------------------
### Check & load required packages
required_packages <- c("unmarked","AICcmodavg","MuMIn","lattice","ggcorrplot",
                       "caret","dplyr","corrplot","mecofun","pROC","readxl")
for(pkg in required_packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    stop(sprintf("Package '%s' missing. Install with install.packages('%s')", pkg, pkg))
  }
  library(pkg, character.only = TRUE)
}

### Source custom functions
source(file.path(fun_script_path, "WindSel.R"))
source(file.path(fun_script_path, "SiteCovSel.R"))

# Ensure output directory exists
if (!dir.exists(result_data_path)) dir.create(result_data_path, recursive = TRUE)

### Import raw data
dataraw <- read_xlsx(file.path(raw_data_path, input_excel_file))

# ------------------------------------------------------------------------------
## Step 2: Data Preprocessing (Shared for Both Models)
## Core Goal: Extract response/observation covariates, define target variables
# ------------------------------------------------------------------------------
### Extract response variables (20 repeated occupancy observations)
y <- dataraw[,paste0("ob", 1:n_observations)]

### Extract observation‑level covariates (4 wind metrics, 20 observations each)
obsCovs <- list(
  wind1 = dataraw[,paste0("wind1.", 1:20)],
  wind2 = dataraw[,paste0("wind2.", 1:20)],
  wind3 = dataraw[,paste0("wind3.", 1:20)],
  wind4 = dataraw[,paste0("wind4.", 1:20)]
)

### Define full site covariate list (including identifiers + env variables)
var_names_SC <- c('X', 'camera', 'sitecode', var_names)

### Subset data for variable selection
xSel <- dataraw[,var_names]          # Env variables only
xSel_SC <- dataraw[,var_names_SC]    # Env + site identifiers

### Filter collinear site covariates (r < collinear_threshold)
predTest <- unmarkedFrameOccu(y=y, siteCovs=xSel, obsCovs=obsCovs)
var_sel <- select07_unm(predTest, names(predTest@siteCovs), threshold = collinear_threshold)
pred_sel <- var_sel$pred_sel  # Weakly correlated predictors
xraw <- dataraw[,pred_sel]

### variable mapping table: original variable name ↔ short name s1‑sN
var_num <- paste0("s", 1:ncol(xraw))
var_mapping <- data.frame(original_var = pred_sel, short_name = var_num)

# ------------------------------------------------------------------------------
## Module A: Standardized Model (0‑1 Min‑Max Scaling)
## Core Goal: Fit model with covariates scaled to 0‑1
# ------------------------------------------------------------------------------
### 1. Rename & standardize predictors (min‑max scaling to 0‑1)
x_sd <- xraw
names(x_sd) <- paste0("s", 1:ncol(x_sd))

# Min‑max standardization, handle zero‑range constant variable
for (col in colnames(x_sd)) {
  col_vec <- x_sd[[col]]
  rng <- max(col_vec, na.rm = TRUE) - min(col_vec, na.rm = TRUE)
  if(rng > 0){
    x_sd[[col]] <- (col_vec - min(col_vec, na.rm = TRUE)) / rng
  }else{
    warning(sprintf("Variable %s has zero range; kept raw values", col))
    x_sd[[col]] <- col_vec
  }
}

### 2. Select best wind covariates (custom function)
WindSel(x_sd)
WindBest_sd <- '~  wind1 + wind2 + wind4  ~' # Read from Line 140 output 

### 3. Create unmarked frame (standardized)
predTS_sd  <- unmarkedFrameOccu(y=y, siteCovs=x_sd, obsCovs=obsCovs)

### 4. Select best site covariates & export selection results
output_sd <- SiteCovSel(x_sd, WindBest_sd)
result_str_sd <- c(
  "=== Standardized Model Results (0‑1 Scaling) ===",
  "Detection Covariate Selection:", WindBest_sd, "",
  "Site Covariate (r<0.7):", var_sel$pred_sel, "",
  paste(spe, "STANDARDIZED RESULT ##########"), "",
  capture.output(print(output_sd)), ""
)
writeLines(result_str_sd, file.path(result_data_path, paste0(spe, "_multi_variable_result.txt")))
save(output_sd, file = file.path(result_data_path, paste0(spe, "_output.RData")))

### 5. Fit top models
top_models_sd <- output_sd$top_models
preds_sd <- list()
for (top_model in top_models_sd) {
  formula_str <- paste(unlist(strsplit(top_model, ",")), collapse = " + ")
  full_formula_str <- paste(WindBest_sd, formula_str)
  pred <- occu(as.formula(full_formula_str), data = predTS_sd)
  preds_sd <- append(preds_sd, list(pred))
}

### 6. Weighted average model
aic_avg_sd <- model.avg(preds_sd)
save(aic_avg_sd, file = file.path(result_data_path, paste0(spe, "_aic_avg_model.RData")))

# Save 95% confidence interval to file
sink(file.path(result_data_path, paste0(spe,"_aic_avg_ci.txt")))
print(confint(aic_avg_sd, method = 'normal', level = 0.95))
sink()

### 7. Predict Psi (occupancy) & P (detection)
psi_aicavg_sd <- predict(aic_avg_sd, type="state")  # Occupancy (Psi)
psi_aicavg_sd <- data.frame(psi_fit = psi_aicavg_sd$fit, psi_se.fit = psi_aicavg_sd$se.fit)
p_aicavg_sd <- predict(aic_avg_sd, type="det")  # Detection (P)
p_aicavg_sd <- data.frame(p_fit = p_aicavg_sd$fit, p_se.fit = p_aicavg_sd$se.fit)

### 8. Export standardized results
psi_Env_sd <- cbind.data.frame(xSel_SC , psi_aicavg_sd)
write.csv(psi_Env_sd, file.path(result_data_path, paste0(spe,'_psi_Env_aicavg.csv')), row.names=F)
write.csv(p_aicavg_sd, file.path(result_data_path, paste0(spe,'_p_aicavg.csv')), row.names=F)

### 9. Stratified K‑fold cross‑validation (AUC performance)
## Derive true observed occupancy status per site (any detection = occupied)
k_fold <- 4

true_psi_all <- apply(y, 1, function(x) {
  valid_values <- na.omit(x)
  if (length(valid_values) == 0) return(NA)
  if (any(valid_values == 1)) {
    return(1)
  } else {
    return(0)
  }
})
true_psi_all <- factor(true_psi_all, levels = c(0, 1))

pos_indices <- which(true_psi_all == 1)
neg_indices <- which(true_psi_all == 0)

folds <- vector("list", k_fold)
pos_folds <- createFolds(pos_indices, k = k_fold, list = TRUE)
neg_folds <- createFolds(neg_indices, k = k_fold, list = TRUE)

for (i in 1:k_fold) {
  folds[[i]] <- c(pos_folds[[i]], neg_folds[[i]])
}

auc_scores <- numeric(k_fold)
for (i in 1:k_fold) {
  train_indices <- unlist(folds[-i])
  test_indices  <- unlist(folds[i])
  
  train_umf <- predTS_sd[train_indices, ]
  test_umf  <- predTS_sd[test_indices, ]
  
  pred_psi <- predict(aic_avg_sd, type = "state", newdata = test_umf)$fit
  true_psi <- as.numeric(apply(test_umf@y, 1, function(x) any(x == 1)))
  auc_scores[i] <- auc(response = true_psi, predictor = pred_psi)
}

# Print CV‑AUC metrics in console
cat("\n==== K‑fold Cross‑Validation AUC Results ====\n")
cat(paste0("Individual fold AUC: ", paste(round(auc_scores,4), collapse = ", "), "\n"))
cat(paste0("Mean AUC: ", round(mean(auc_scores),4), "\n"))
cat(paste0("Variance of AUC: ", round(var(auc_scores),4), "\n"))

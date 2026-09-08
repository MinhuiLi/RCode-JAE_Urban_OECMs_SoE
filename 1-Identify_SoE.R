# ==============================================================================
# Occupancy Model for Scale of Effect Analysis
# Script: 1_Identify_SoE.R
#
# Description:
#   Fit single‑species occupancy models across multiple spatial scales;
#   export model outputs to txt; parse txt to AIC‑weight excel;
#   extract Equivalent Optimal Models (EOM) and compute Scale‑of‑Effect (SOE).
#
# Usage:
#   !! RUN ONE SPECIES AT A TIME !!
#   Manually set parameters in CONFIG section below for target species.
#   To batch multiple species: wrap full script inside a higher‑level loop
#   (see commented example in CONFIG).
#
# Input:  CE_Pallas_s_squirrel_env_detection_data.xlsx (species‑specific table)
# Output:
#   1. <spe>_<env>_AIC_estimates50_2000.txt      raw model output text
#   2. <spe>_AIC_weight_50_2000.xlsx             parsed AICc & weight table
#   3. <spe>_EOS_50_2000.xlsx                    EOM & SOE summary
#
# Notes for public release (GitHub/Zenodo):
#   - User MUST update raw_data_path, result_data_path to local file paths
#   - Input excel filename and target species code are defined in CONFIG
#   - This script is step‑1 of sequential analysis workflow.
# ==============================================================================

# ------------------------------------------------------------------------------
# CONFIGURATION: ALL USER‑MODIFIABLE PARAMETERS HERE
# !! CHANGE THESE BEFORE RUNNING FOR A NEW SPECIES !!
# ------------------------------------------------------------------------------
## File paths (set to your local directories)
raw_data_path  <- "."   # e.g., "./data/raw"
result_data_path <- "." # e.g., "./results/occupancy"

## Species & input file
spe <- "CE"
input_excel_file <- "CE_Pallas_s_squirrel_env_detection_data.xlsx"

## Analysis parameters
env_list        <- c('CH', 'NDVI', 'WSR','WD', 'RD', 'BD', 'ISR',
                     'NP', 'PD', 'MPS', 'ED', 'CO')   # env vars for model fitting
env_list_full   <- c('CH', 'NDVI', 'WSR','WD', 'RD', 'BD', 'ISR',
                     'NP', 'PD', 'MPS', 'ED', 'CO')   # full set for EOM parsing
scales          <- seq(50, 2000, by = 50)             # spatial scale sequence
n_scales        <- length(scales)
n_observations  <- 20                                 # repeated detection occasions
EOM_delta_thresh <- 2                                 # delta AICc threshold for EOM

## Detection formula: occu(~ detection_formula ~ occupancy_formula)
WindBest <- "~ 1 ~"

# ------------------------------------------------------------------------------
# Module 1: Load Packages, Import Data & Set Parameters
# Core Goal: Load dependencies, import raw data, define analysis parameters
# ------------------------------------------------------------------------------
## Load core packages (grouped by functionality)
### Occupancy modeling
library(unmarked)       # Unmarked occupancy model fitting
library(AICcmodavg)     # AICc calculation & model selection
library(MuMIn)          # Multi‑model inference & averaging
### Data processing
library(dplyr)          # Data cleaning/manipulation
library(readxl)         # Excel import
library(writexl)        # Excel export
library(tidyverse)      # Integrated data tools


# Ensure output directory exists
if (!dir.exists(result_data_path)) dir.create(result_data_path, recursive = TRUE)

## Import raw data
dataraw <- read_excel(file.path(raw_data_path, input_excel_file))

## Extract key variables for occupancy modeling
### Response variables (repeated observations)
y <- dataraw[, paste0("ob", 1:n_observations)]

### Observation‑level covariates (wind variables)
obsCovs <- list(
  wind1 = dataraw[, paste0("wind1.", 1:n_observations)],
  wind2 = dataraw[, paste0("wind2.", 1:n_observations)],
  wind3 = dataraw[, paste0("wind3.", 1:n_observations)],
  wind4 = dataraw[, paste0("wind4.", 1:n_observations)]
)

# ------------------------------------------------------------------------------
# Module 2: Fit Models & Identify Scale of Effect (Export TXT)
# Core Goal: Batch fit occupancy models across scales, export results to TXT
# ------------------------------------------------------------------------------
## Define batch processing function
process_env <- function(env) {
  ### Generate dynamic env variable names (e.g., WSR50, WSR100)
  var_names <- paste0(env, scales)
  xraw <- dataraw[, var_names]
  s <- `colnames<-`(xraw, paste0("s", seq_along(scales)))
  
  ### Create unmarked frame (standard format for occupancy models)
  predUMF <- unmarkedFrameOccu(y = y, siteCovs = s, obsCovs = obsCovs)
  
  ### Fit occupancy models for each scale
  model_list <- lapply(seq_along(scales), function(i) {
    formula_str <- paste0(WindBest, paste0("s", i))
    occu(as.formula(formula_str), data = predUMF)
  })
  
  ### Perform model averaging & generate report
  avg_model <- model.avg(model_list)
  generate_report <- function() {
    c(
      paste(spe, "_", env, "WEIGHT ##########"), "",
      capture.output(summary(avg_model)), "",
      paste(spe, "_", env, "ESTIMATE ##########"), "",
      sapply(seq_along(model_list), function(i) {
        c(paste("### pred", scales[i], "estimate:"),
          capture.output(show(model_list[[i]]@estimates)), "")
      })
    )
  }
  
  ### Export report to TXT file
  out_txt <- file.path(result_data_path, paste0(spe, "_", env, "_AIC_estimates50_2000.txt"))
  writeLines(generate_report(), out_txt)
}

## Execute batch model fitting
lapply(env_list, process_env)

# ------------------------------------------------------------------------------
# Module 3: Parse Model Output (TXT → Clean Excel)
# Core Goal: Extract AICc/delta/weight from TXT, organize into structured Excel
#
# WARNING: text‑parsing is fragile; depends on printed output format of MuMIn.
# Do not modify printed console output format when running.
# ------------------------------------------------------------------------------
## Switch to result path
setwd(result_data_path)

## Extract metrics from TXT files
### Get all TXT result files
file_paths <- list.files(pattern = "\\.txt$", full.names = TRUE)

psi_s_sorted <- sort(paste0('s', 1:n_scales))
result_df <- data.frame()

### Parse each TXT file
for (file_path in file_paths) {
  lines <- readLines(file_path)
  # Extract species‑env variable name
  aic_line <- grep(paste0(spe, ".* WEIGHT"), lines)[1]
  variable_name <- sub(paste0(spe, "(.*) WEIGHT.*"), "\\1", lines[aic_line])
  row_name <- paste0(spe, variable_name)
  
  # Extract Component Models section
  comp_models_start <- grep("Component models:", lines)[1] + 2
  comp_models_lines <- lines[comp_models_start:(comp_models_start + 99)]
  
  # Parse model metrics (df/logLik/AICc/delta/weight)
  for (line in comp_models_lines) {
    parts <- strsplit(line, "\\s+")[[1]]
    if (length(parts) >= 6) {
      temp_df <- data.frame(
        Species_variable = row_name,
        Component_Model = psi_s_sorted[as.numeric(parts[1])],
        df = parts[2],
        logLik = parts[3],
        AICc = parts[4],
        delta = parts[5],
        weight = parts[6],
        stringsAsFactors = FALSE
      )
      result_df <- rbind(result_df, temp_df)
    }
  }
}


scale_values_eom <- seq(50, 2000, by = 50)
df <- result_df %>%
  mutate(scale = scale_values_eom[match(Component_Model, paste0("s", 1:length(scale_values_eom)))])
dfc <- df[!is.na(df$scale), ]

## Export cleaned data to Excel
write_xlsx(dfc, paste0(spe, "_AIC_weight_50_2000.xlsx"))

# ------------------------------------------------------------------------------
# Module 4: Extract Equivalent Optimal Models (EOMs) & Calculate SOE
# Core Goal: Filter EOMs (delta ≤2), compute Scale of Effect (SOE)
# ------------------------------------------------------------------------------
## Import cleaned Excel data (output from Module 3)
data <- read_excel(file.path(result_data_path, paste0(spe, "_AIC_weight_50_2000.xlsx")))

## Initialize dataframe for EOM results
eom_result_df <- data.frame()

## Iterate over each env variable to identify EOMs
for (env in env_list_full) {
  spe_env <- paste0(spe, " _ ", env)
  dataEnv <- data[data$Species_variable == spe_env, ] %>% arrange(scale)
  
  # Filter EOMs (delta AICc ≤ threshold defined in CONFIG)
  EOM <- dataEnv[as.numeric(dataEnv$delta) <= EOM_delta_thresh, ]
  
  # Calculate SOE & EOM range (if EOMs exist)
  if (nrow(EOM) > 0) {
    temp_df <- data.frame(
      Species = spe,
      Env = env,
      Species_variable = spe_env,
      SOE = EOM$scale[which.min(EOM$AICc)],  # Optimal scale (min AICc)
      EOSmin = min(EOM$scale),              # Min EOM scale
      EOSmax = max(EOM$scale),              # Max EOM scale
      EOSs = paste(EOM$scale, collapse = ", "),  # All EOM scales
      NO.EOSs = nrow(EOM),                  # Number of EOMs
      stringsAsFactors = FALSE
    )
    eom_result_df <- rbind(eom_result_df, temp_df)
  }
}

## Print & export EOM results
print(eom_result_df)
write_xlsx(eom_result_df, paste0(spe, "_EOS_50_2000.xlsx"))

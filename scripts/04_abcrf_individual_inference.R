# ==============================================================================
# Individual-Level Model Inference Using ABC Random Forest
# ==============================================================================
# AUTHOR: Yitong
# 
# PURPOSE:
#   - Train ABC-RF classifier on simulated data
#   - Infer individual posterior model probabilities
#   - Perform hierarchical Bayesian Model Selection (RFX-BMS)
#   - Generate figures and tables for manuscript
#
# WORKFLOW:
#   Part I: Primary dataset analysis (Sections 0-8)
#   Part II: Stengard replication dataset (Sections 9-13)
# ==============================================================================

###############################################################
# SECTION 0: Load Packages
###############################################################

pacman::p_load(
  # Data wrangling
  dplyr, tidyr, stringr, magrittr, data.table, tibble,knitr,kableExtra,
  # Visualisation
  ggplot2, ggpubr, gridExtra, cowplot, forcats, scales,
  # ABC-RF
  ranger, abcrf,
  # Parallel
  foreach, doParallel, parallel,
  # File I/O
  readxl,readr,
  # For rowSums2 etc. (though often part of base or loaded by others)
  matrixStats,devtools,patchwork
)

###############################################################
# SECTION 1: Setup Parallel Computing
###############################################################
# Reserve 4 cores for system, use remaining for ABC-RF computation

num_tasks <- parallel::detectCores(logical = TRUE) - 4L
if (num_tasks < 1) num_tasks <- 1L
cl <- parallel::makeCluster(num_tasks)
registerDoParallel(cl)

###############################################################
# SECTION 2: Source Custom Functions
###############################################################
# VB_bms.R: Random-effects Bayesian Model Selection
# abcrf_helpers.R: ABC-RF utility functions

source("R/VB_bms.R")
source("R/abcrf_helpers.R")
source("R/plot_helpers.R")
###############################################################
# SECTION 3: Load Data (Primary Dataset)
###############################################################
# simulations: Simulated summary statistics with known model labels
# observed: Human data summary statistics (model labels to be inferred)

simulations <- as.data.table(readRDS("data/Simulate_Summary_dt.rds"))
observed    <-  as.data.table(readRDS("data/Human_Summary_dt.rds"))
plot_dir_tiff    <- "fig/tiff"
plot_dir_png    <- "fig/png"
default_dpi <- 300
plot_width  <- 8
plot_height <- 5
###############################################################
# SECTION 4: Load Trial-Level Data
###############################################################
# Import experimental data for format-specific analyses

raw_path <- "data/df.csv"
S_path   <- "data/s.csv"

df <- load_clean_data(raw_path, type = "experiment")
s  <- load_clean_data(S_path,   type = "stengard")

simulations$model <- factor(simulations$model)

###############################################################
# SECTION 5: Train ABC-RF Classifier
###############################################################
# Formula: 40 summary statistics predict model labels
#   - Regression slopes/intercepts (18): BR/HR/FAR sensitivity
#   - Empirical divergence (7): EMP_B, EMP_BO, EMP_HO, etc.
#   - Probability distortion (7): PD_*
#   - Additivity indices (7): Ad_*
#   - Response variance (1): mean_variance
# Hyperparameters: 500 trees, parallel processing enabled

model_formula <- factor(model) ~
  slope_BR0.01 + slope_BR0.4 + slope_BR0.97 +
  slope_HR0.53 + slope_HR0.72 + slope_HR0.9 +
  slope_FAR0.11 + slope_FAR0.31 + slope_FAR0.42 +
  intercept_BR0.01 + intercept_BR0.4 + intercept_BR0.97 +
  intercept_HR0.53 + intercept_HR0.72 + intercept_HR0.9 +
  intercept_FAR0.11 + intercept_FAR0.31 + intercept_FAR0.42 +
  EMP_B + EMP_BO + EMP_HO + EMP_FO + EMP_JO + EMP_LS + EMP_H +
  PD_B + PD_BO + PD_HO + PD_FO + PD_JO + PD_LS + PD_H +
  Ad_B + Ad_BO + Ad_HO + Ad_FO + Ad_JO + Ad_LS + Ad_H +
  mean_variance


abc_model <- abcrf::abcrf(
  formula = model_formula,
  data    = simulations,
  ntree   = 500,
  ncores  = num_tasks,
  paral   = FALSE
)

###############################################################
# SECTION 5.1: OOB Validation
###############################################################
# Evaluate classifier performance using Out-of-Bag predictions

pred_oob <- predict(abc_model$model.rf, data = simulations, type = "response")$predictions
true_mod <- simulations$model 

# Standardize model names and assign to families
mods_in_data <- sort(unique(c(as.character(true_mod), as.character(pred_oob))))
std_names <- clean_model_names(mods_in_data)   
names(std_names) <- mods_in_data               
fam_std <- assign_family(std_names)     
fam_map_full <- fam_std
names(fam_map_full) <- names(std_names) 

# Compute family-level OOB metrics
res_fam <- family_oob_metrics(pred_oob, true_mod, fam_map_full)
res_fam$confusion
res_fam$err_micro
res_fam$err_macro
res_fam$prior_chance

###############################################################
# SECTION 6: Compute Posterior Probabilities
###############################################################
# Apply trained ABC-RF to observed data

posterior_matrix <- compute_posterior(observed, simulations, abc_model)
model_names_from_posterior <- clean_model_names(colnames(posterior_matrix))
posterior_long <- to_long_posterior(posterior_matrix, observed)

# Merge with format information
df_unique <- df %>%
  distinct(subject, format) %>%
  mutate(subject = as.character(subject))

posterior_long <- posterior_long %>%
  mutate(id = as.character(id)) %>%
  left_join(df_unique, by = c("id" = "subject")) 

saveRDS(posterior_long, "data/posterior_long.rds")

###############################################################
# SECTION 7: Hierarchical Bayesian Model Selection (RFX-BMS)
###############################################################
# Transform posteriors to log-Bayes factors
# Run two-stage BMS separately for probability and frequency formats

post_p <- build_matrix(posterior_long, "probability")
post_f <- build_matrix(posterior_long, "frequency")

K_models <- ncol(posterior_matrix)
eps      <- 1e-10

m_p <- log((post_p /  (1 / K_models)) + eps)
m_f <- log((post_f / (1 / K_models)) + eps)

# Execute RFX-BMS
out_p <- two_stage_bms(
  m_p          = m_p,
  model_names  = model_names_from_posterior,
  fam_map      = fam_map,
  n_samples    = 1e6
)

out_f <- two_stage_bms(
  m_p          = m_f,
  model_names  = model_names_from_posterior,
  fam_map      = fam_map,
  n_samples    = 1e6
)

###############################################################
# SECTION 8: Generate Figures (Experiment)
###############################################################

# Display results
print_family_table(out_p, title = "Prob-based Family BMS")
print_model_table(out_p,  title = "Prob-based Model BMS")

print_family_table(out_f, title = "Freq-based Family BMS")
print_model_table(out_f,  title = "Freq-based Model BMS")

# Create plots
pa <- plot_two_stage_bms(out_p, title = "Probability Format")
pb <- plot_two_stage_bms(out_f, title = "Frequency Format")
final_plot <- (pa + pb) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Figure_13.tiff"), final_plot,
       width = 8, height = 5, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_13.png"), final_plot,
       width = 8, height = 5, dpi = 300)

###############################################################
# ==============================================================================
# PART II: REPLICATION ANALYSIS (STENGARD DATASET)
# ==============================================================================
# Repeat entire pipeline for independent replication dataset
###############################################################

###############################################################
# SECTION 9: Load Data (Stengard Dataset)
###############################################################

simulations_s <- as.data.table(readRDS("data/Simulate_Summary_dts.rds"))
observed_s    <-  as.data.table(readRDS("data/Human_Summary_dts.rds"))

df_s <- read_csv("data/s.csv") |>
  transmute(
    subject, format,trial,
    true_posterior,response,
    BR = br,
    HR = hr,
    FAR = far,
    subject_s = paste0(subject, "_", ifelse(format == "frequency", 1, 2))
  )

simulations_s$model <- factor(simulations_s$model)

###############################################################
# SECTION 10: Train ABC-RF (Stengard Dataset)
###############################################################
# Note: Different BR/HR/FAR levels than primary dataset
# Formula uses 43 predictors (no mean_variance term)

model_formula_s <- factor(model) ~ 
  slope_BR0.1+slope_BR0.3+slope_BR0.5+slope_BR0.7+slope_BR0.9+
  intercept_BR0.1+intercept_BR0.3+intercept_BR0.5+
  intercept_BR0.7+intercept_BR0.9+slope_HR0.5+
  slope_HR0.7+slope_HR0.9+intercept_HR0.5+intercept_HR0.7+
  intercept_HR0.9+slope_FAR0.1+slope_FAR0.3+slope_FAR0.5+
  intercept_FAR0.1+intercept_FAR0.3+intercept_FAR0.5+ 
  EMP_B+EMP_BO + EMP_HO + EMP_FO + EMP_JO + EMP_LS +EMP_H+
  PD_B +PD_BO + PD_HO + PD_FO + PD_JO +PD_LS+ PD_H+
  Ad_B + Ad_BO + Ad_HO + Ad_FO + Ad_JO + Ad_LS + Ad_H 

abc_model_s <- abcrf::abcrf(
  formula = model_formula_s,
  data    = simulations_s,
  ntree   = 500,
  ncores  = num_tasks,
  paral   = FALSE
)

###############################################################
# SECTION 10.1: OOB Validation (Stengard)
###############################################################

pred_oob_s <- predict(abc_model_s$model.rf, data = simulations_s, type = "response")$predictions
true_mod_s <- simulations_s$model 

mods_in_data_s <- sort(unique(c(as.character(true_mod_s), as.character(pred_oob_s))))
std_names_s <- clean_model_names(mods_in_data_s)   
names(std_names_s) <- mods_in_data_s               
fam_std_s <- assign_family(std_names_s)     
fam_map_full_s <- fam_std_s
names(fam_map_full_s) <- names(std_names_s) 

res_fam_s <- family_oob_metrics(pred_oob_s, true_mod_s, fam_map_full_s)
res_fam_s$confusion
res_fam_s$err_micro
res_fam_s$err_macro
res_fam_s$prior_chance

###############################################################
# SECTION 11: Compute Posteriors (Stengard)
###############################################################

posterior_matrix_s <- compute_posterior(observed_s, simulations_s, abc_model_s)
model_names_from_posterior_s <- clean_model_names(colnames(posterior_matrix_s))
posterior_long_s <- to_long_posterior(posterior_matrix_s, observed_s)

df_unique_s <- df_s %>%
  distinct(subject_s, format) %>%
  mutate(subject = as.character(subject_s))

posterior_long_s <- posterior_long_s %>%
  mutate(id = as.character(id)) %>%
  left_join(df_unique_s, by = c("id" = "subject_s")) 

saveRDS(posterior_long_s, "data/posterior_long_s.rds")

###############################################################
# SECTION 12: RFX-BMS (Stengard)
###############################################################

post_s_p <- build_matrix(posterior_long_s, "probability")
post_s_f <- build_matrix(posterior_long_s, "frequency")

K_models_s <- ncol(posterior_matrix_s)
eps      <- 1e-10

m_s_p <- log((post_s_p /  (1 / K_models_s)) + eps)
m_s_f <- log((post_s_f / (1 / K_models_s)) + eps)

out_s_p <- two_stage_bms(
  m_p          = m_s_p,
  model_names  = model_names_from_posterior_s,
  fam_map      = fam_map,
  n_samples    = 1e6
)
out_s_f <- two_stage_bms(
  m_p          = m_s_f,
  model_names  = model_names_from_posterior_s,
  fam_map      = fam_map,
  n_samples    = 1e6
)

###############################################################
# SECTION 13: Generate Figures (Stengard)
###############################################################
# Output: Figure_6.png + LaTeX tables

# Display results
print_family_table(out_s_p, title = "Prob-based Family BMS")
print_model_table(out_s_p,  title = "Prob-based Model BMS")

print_family_table(out_s_f, title = "Freq-based Family BMS")
print_model_table(out_s_f,  title = "Freq-based Model BMS")

# Create plots
c <- plot_two_stage_bms(out_s_p, title = "Probability Format")
d <- plot_two_stage_bms(out_s_f, title = "Frequency Format")
final_plot_s <- (c + d) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Figure_6.tiff"), final_plot_s,
       width = 8, height = 5, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_6.png"), final_plot_s,
       width = 8, height = 5, dpi = 300)

###############################################################
# END OF SCRIPT
###############################################################
# Output files:
#   - posterior_long.rds, posterior_long_s.rds
#   - Figure_11.png (primary), Figure_6.png (Stengard)
#   - family_bms_*.tex, model_bms_*.tex (8 tables total)
# Note: Close parallel cluster with parallel::stopCluster(cl) when done
###############################################################
###############################################################
# 0 · PACKAGE SET-UP  (explicit namespaces → fewer clashes)
###############################################################
pacman::p_load(
  # Data wrangling
  dplyr, tidyr, stringr, magrittr, data.table, tibble,
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
# 1. PARALLEL BACKEND INITIALIZATION
###############################################################
num_tasks <- parallel::detectCores(logical = TRUE) - 4L   # keep 4 threads free
if (num_tasks < 1) num_tasks <- 1L # Ensure at least one task
cl <- parallel::makeCluster(num_tasks)
registerDoParallel(cl)
###############################################################
# 2. SOURCE CUSTOM FUNCTIONS
# Import project-specific utilities for model fitting and ABC-RF
###############################################################
source("R/VB_bms.R")
source("R/abcrf_helpers.R")
###############################################################
# 3. DATA LOADING
# Load both simulated and human summary data
###############################################################
simulations <- as.data.table(readRDS("data/Simulate_Summary_dt.rds"))
observed    <-  as.data.table(readRDS("data/Human_Summary_dt.rds"))

###############################################################
# 4. EXPERIMENTAL DATA PREPARATION
# Read trial-level human responses and merge key variables
###############################################################
df <- readxl::read_excel("data/df.xlsx") |>
  transmute(
    subject, format, n_trial, rt,
    BR = br,
    HR = hr,
    FAR = far,
    response = round(inputvalue / 100, 2),
    true_posterior = round(correctAnswer / 100, 2)
  )

simulations$model <- factor(simulations$model)
###############################################################
# 5. BUILD ABC-RF CLASSIFIER
# Train random forest on simulation summaries to predict model labels

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


abc_model <- abcrf::abcrf( # Use abcrf:: for clarity
  formula = model_formula,
  data    = simulations,
  ntree   = 500,
  ncores  = num_tasks, # ranger parallelism
  paral   = FALSE      # abcrf's own parallelism, ranger handles it
)
###############################################################

###############################################################
# 6. COMPUTE POSTERIOR PROBABILITIES
# Apply trained ABC-RF to observed data to get posterior matrix

posterior_matrix <- compute_posterior(observed, simulations, abc_model)
model_names_from_posterior <- clean_model_names(colnames(posterior_matrix))

posterior_long <- to_long_posterior(posterior_matrix, observed)
df_unique <- df %>%
  distinct(subject, format) %>%
  mutate(subject = as.character(subject))

posterior_long <- posterior_long %>%
  mutate(id = as.character(id)) %>%
  left_join(df_unique, by = c("id" = "subject")) 
###############################################################

###############################################################
# 7. MODEL-LEVEL RANDOM EFFECTS BMS
# 7.1 Build log-evidence (log Bayes factor) matrices for each format

post_p <- build_matrix(posterior_long, "probability")
post_f <- build_matrix(posterior_long, "frequency")

# 7.1 · Build pseudo log-evidence matrix (log Bayes factors) m_models
K_models <- ncol(posterior_matrix)
eps      <- 1e-10 # Small constant for numerical stability (avoid log(0))
# Log Bayes Factor for each model vs. a uniform model prior (1/K_models)
m_p <- log((post_p /  (1 / K_models)) + eps)
m_f <- log((post_f / (1 / K_models)) + eps)

# 7.2 · Run model-level RFX-BMS
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

###############################################################
# 8. VISUALIZE MODEL FIT
# Generate side-by-side BMS plots for probability and frequency formats
print_family_table(out_p, title = "Prob-based Family BMS")
print_model_table(out_p,  title = "Prob-based Model BMS")

print_family_table(out_f, title = "Freq-based Family BMS")
print_model_table(out_f,  title = "Freq-based Model BMS")

pa <- plot_two_stage_bms(out_p, title = "Probability Format")
pb <- plot_two_stage_bms(out_f, title = "Frequency Format")


final_plot <- (pa + pb) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

# ========== Figure 6: Model Comparison (Model Fitting) for Experiment ==========

ggsave("fig/Figure_10.png", final_plot,width = 10, height = 6, dpi = 300)

###############################################################

###############################################################
# 9. stangard DATA ANALYSIS
# Repeat ABC-RF and BMS pipeline for secondary (Stengard) dataset
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

abc_model_s <- abcrf::abcrf( # Use abcrf:: for clarity
  formula = model_formula_s,
  data    = simulations_s,
  ntree   = 500,
  ncores  = num_tasks, # ranger parallelism
  paral   = FALSE      # abcrf's own parallelism, ranger handles it
)

posterior_matrix_s <- compute_posterior(observed_s, simulations_s, abc_model_s)
model_names_from_posterior_s <- clean_model_names(colnames(posterior_matrix_s))
posterior_long_s <- to_long_posterior(posterior_matrix_s, observed_s)

df_unique_s <- df_s %>%
  distinct(subject_s, format) %>%
  mutate(subject = as.character(subject_s))

posterior_long_s <- posterior_long_s %>%
  mutate(id = as.character(id)) %>%
  left_join(df_unique_s, by = c("id" = "subject_s")) 

post_s_p <- build_matrix(posterior_long_s, "probability")
post_s_f <- build_matrix(posterior_long_s, "frequency")

K_models_s <- ncol(posterior_matrix_s)
eps      <- 1e-10 # Small constant for numerical stability (avoid log(0))

# Log Bayes Factor for each model vs. a uniform model prior (1/K_models)
m_s_p <- log((post_s_p /  (1 / K_models_s)) + eps)
m_s_f <- log((post_s_f / (1 / K_models_s)) + eps)

# 7.3 · Assemble and display Family and model-level BMS results
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

print_family_table(out_s_p, title = "Prob-based Family BMS")
print_model_table(out_s_p,  title = "Prob-based Model BMS")

print_family_table(out_s_f, title = "Freq-based Family BMS")
print_model_table(out_s_f,  title = "Freq-based Model BMS")
# ========== Figure 8: Model Comparison (Model Fitting) for Experiment ==========
c <- plot_two_stage_bms(out_s_p, title = "Probability Format")
d <- plot_two_stage_bms(out_s_f, title = "Frequency Format")

final_plot_s <- (c + d) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

# ========== Figure 6: Model Comparison (Model Fitting) for Experiment ==========

ggsave("fig/Figure_6.png", final_plot_s,width = 10, height = 6, dpi = 300)

###############################################################
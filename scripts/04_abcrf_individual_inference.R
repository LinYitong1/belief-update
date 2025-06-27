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
  matrixStats
)

###############################################################
# 1 · PARALLEL BACKEND (Windows-safe PSOCK)
###############################################################
num_tasks <- parallel::detectCores(logical = TRUE) - 4L   # keep 4 threads free
if (num_tasks < 1) num_tasks <- 1L # Ensure at least one task
cl <- parallel::makeCluster(num_tasks)
registerDoParallel(cl)

source("R/VB_bms.R")
source("R/abcrf_helpers.R")
################################################################
# 2 · LOAD DATA
################################################################
simulations <- as.data.table(readRDS("data/Simulate_Summary_dt.rds"))
observed    <-  as.data.table(readRDS("data/Human_Summary_dt.rds"))

###############################################################
# 3 · Experiment Data
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
################################################################
# 3 · FIT ABC-RF CLASSIFIER
################################################################
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

posterior_matrix <- compute_posterior(observed, simulations, abc_model)
posterior_long <- to_long_posterior(posterior_matrix, observed)

df_unique <- df %>%
  distinct(subject, format) %>%
  mutate(subject = as.character(subject))

posterior_long <- posterior_long %>%
  mutate(id = as.character(id)) %>%
  left_join(df_unique, by = c("id" = "subject")) 


post_p <- build_matrix(posterior_long, "probability")
post_f <- build_matrix(posterior_long, "frequency")

# 7.1 · Build pseudo log-evidence matrix (log Bayes factors) m_models
K_models <- ncol(posterior_matrix)
eps      <- 1e-10 # Small constant for numerical stability (avoid log(0))
# Log Bayes Factor for each model vs. a uniform model prior (1/K_models)
m_p <- log((post_p /  (1 / K_models)) + eps)
m_f <- log((post_f / (1 / K_models)) + eps)

# 7.2 · Run model-level RFX-BMS
out_models_p <- VB_bms_1(m_p, n_samples = 1e5) # n_samples for calculating exceedance probabilities
out_models_f <- VB_bms_1(m_f, n_samples = 1e5) 

# 7.3 · Assemble and display model-level BMS results
#RENAME Model Name for understand
model_names_from_posterior <- colnames(posterior_matrix) %>%
  str_remove("_NH$") %>%
  str_replace("^BS_R$", "A_BS") %>%
  str_replace("B$", "Bayes") %>%
  str_replace("M_L$", "Two Stage") %>%
  str_replace("^MIN", "Hybrid")

df_models_bms_p <- tibble(
  model = model_names_from_posterior,
  alpha = out_models_p$alpha, # Parameters of the Dirichlet distribution
  r     = out_models_p$r,     # Estimated model frequencies in the population
  xp    = out_models_p$xp,    # Exceedance probabilities
  F1    = out_models_p$F1,    # Evidence of alternative
  F0    = out_models_p$F0,    # Evidence of null (equal model freqs)
  pxp   = out_models_p$pxp    # Protected exceedance probabilities (if available, else same as xp)
) %>%
  arrange(desc(r)) # Sort by estimated model frequency

print("Model-level RFX-BMS for Probability Format Results:")
print(df_models_bms_p)

df_models_bms_p<-df_models_bms_p%>%
  mutate(group = case_when(
    model %in% "Bayes" ~ "Bayes",
    model %in% c("HO","BO","FO","JO","LS","H") ~ "Single Heuristic",
    model %in% c("MH","AH","LE") ~ "Multiple Heuristic",
    model == "LA" ~ "Linear Averaging",
    model %in% c("BS","A_BS") ~ "Bayesian Sampler",
    model == "R" ~ "Random",
    model %in% c("Hybrid_BS_BO","Hybrid_BS_HO","Hybrid_BS_FO",
                 "Hybrid_BS_JO","Hybrid_BS_LS","Hybrid_BS_H","Two Stage") ~ "Hybrid Models"
  ))

df_models_bms_f <- tibble(
  model = model_names_from_posterior,
  alpha = out_models_f$alpha, # Parameters of the Dirichlet distribution
  r     = out_models_f$r,     # Estimated model frequencies in the population
  xp    = out_models_f$xp,    # Exceedance probabilities
  F1    = out_models_f$F1,    # Evidence of alternative
  F0    = out_models_f$F0,    # Evidence of null (equal model freqs)
  pxp   = out_models_f$pxp    # Protected exceedance probabilities (if available, else same as xp)
) %>%
  arrange(desc(r)) # Sort by estimated model frequency

print("Model-level RFX-BMS for Frequency Format Results:")
print(df_models_bms_f)

df_models_bms_f<-df_models_bms_f%>%
  mutate(group = case_when(
    model %in% "Bayes" ~ "Bayes",
    model %in% c("HO","BO","FO","JO","LS","H") ~ "Single Heuristic",
    model %in% c("MH","AH","LE") ~ "Multiple Heuristic",
    model == "LA" ~ "Linear Averaging",
    model %in% c("BS","A_BS") ~ "Bayesian Sampler",
    model == "R" ~ "Random",
    model %in% c("Hybrid_BS_BO","Hybrid_BS_HO","Hybrid_BS_FO",
                 "Hybrid_BS_JO","Hybrid_BS_LS","Hybrid_BS_H","Two Stage") ~ "Hybrid Models"
  ))


# -------- Generate the two plots --------
pa <- plot_bms(df_models_bms_p, "C. Probability Format")
pb <- plot_bms(df_models_bms_f, "D. Frequency Format")

# -------- Combine them --------
final_plot <- pa + pb + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

# -------- Save high-resolution version --------

ggsave("figs/Figure_model_fit_exp.png", final_plot,width = 10, height = 6, dpi = 300)
###############################################################
# 4 · Stengard Data
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

# 7.2 · Run model-level RFX-BMS
out_models_s_p <- VB_bms_1(m_s_p, n_samples = 1e5) # n_samples for calculating exceedance probabilities
out_models_s_f <- VB_bms_1(m_s_f, n_samples = 1e5) 

# 7.3 · Assemble and display model-level BMS results
#RENAME Model Name for understand
model_names_from_posterior_s <- colnames(posterior_matrix_s) %>%
  str_remove("_NH$") %>%
  str_replace("^BS_R$", "A_BS") %>%
  str_replace("B$", "Bayes") %>%
  str_replace("M_L$", "Two Stage") %>%
  str_replace("^MIN", "Hybrid")

df_models_bms_s_p <- tibble(
  model = model_names_from_posterior_s,
  alpha = out_models_s_p$alpha, # Parameters of the Dirichlet distribution
  r     = out_models_s_p$r,     # Estimated model frequencies in the population
  xp    = out_models_s_p$xp,    # Exceedance probabilities
  F1    = out_models_s_p$F1,    # Evidence of alternative
  F0    = out_models_s_p$F0,    # Evidence of null (equal model freqs)
  pxp   = out_models_s_p$pxp    # Protected exceedance probabilities (if available, else same as xp)
) %>%
  arrange(desc(r)) # Sort by estimated model frequency

print("Model-level RFX-BMS for Probability Format Results:")
print(df_models_bms_s_p)

df_models_bms_s_p<-df_models_bms_s_p%>%
  mutate(group = case_when(
    model %in% "Bayes" ~ "Bayes",
    model %in% c("HO","BO","FO","JO","LS","H") ~ "Single Heuristic",
    model %in% c("MH","AH","LE") ~ "Multiple Heuristic",
    model == "LA" ~ "Linear Averaging",
    model %in% c("BS","A_BS") ~ "Bayesian Sampler",
    model == "R" ~ "Random",
    model %in% c("Hybrid_BS_BO","Hybrid_BS_HO","Hybrid_BS_FO",
                 "Hybrid_BS_JO","Hybrid_BS_LS","Hybrid_BS_H","Two Stage") ~ "Hybrid Models"
  ))

df_models_bms_s_f <- tibble(
  model = model_names_from_posterior_s,
  alpha = out_models_s_f$alpha, # Parameters of the Dirichlet distribution
  r     = out_models_s_f$r,     # Estimated model frequencies in the population
  xp    = out_models_s_f$xp,    # Exceedance probabilities
  F1    = out_models_s_f$F1,    # Evidence of alternative
  F0    = out_models_s_f$F0,    # Evidence of null (equal model freqs)
  pxp   = out_models_s_f$pxp    # Protected exceedance probabilities (if available, else same as xp)
) %>%
  arrange(desc(r)) # Sort by estimated model frequency

print("Model-level RFX-BMS for Frequency Format Results:")
print(df_models_bms_s_f)

df_models_bms_s_f<-df_models_bms_s_f%>%
  mutate(group = case_when(
    model %in% "Bayes" ~ "Bayes",
    model %in% c("HO","BO","FO","JO","LS","H") ~ "Single Heuristic",
    model %in% c("MH","AH","LE") ~ "Multiple Heuristic",
    model == "LA" ~ "Linear Averaging",
    model %in% c("BS","A_BS") ~ "Bayesian Sampler",
    model == "R" ~ "Random",
    model %in% c("Hybrid_BS_BO","Hybrid_BS_HO","Hybrid_BS_FO",
                 "Hybrid_BS_JO","Hybrid_BS_LS","Hybrid_BS_H","Two Stage") ~ "Hybrid Models"
  ))


# -------- Generate the two plots --------
pc <- plot_bms(df_models_bms_s_p, "C. Probability Format")
pd <- plot_bms(df_models_bms_s_f, "D. Frequency Format")

# -------- Combine them --------
final_plot_s <- pc + pd + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

# -------- Save high-resolution version --------

ggsave("figs/Figure_model_fit_s.png", final_plot_s,width = 10, height = 6, dpi = 300)

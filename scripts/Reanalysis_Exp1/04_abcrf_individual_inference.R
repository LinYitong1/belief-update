# Dependencies
pacman::p_load(
  dplyr, tidyr, stringr, magrittr, data.table, tibble,knitr,kableExtra,
  ggplot2, ggpubr, gridExtra, cowplot, forcats, scales,
  ranger, abcrf,
  foreach, doParallel, parallel,
  readxl,readr,
  matrixStats,devtools,patchwork
)

# Parallel backend
num_tasks <- parallel::detectCores(logical = TRUE) - 4L
if (num_tasks < 1) num_tasks <- 1L
cl <- parallel::makeCluster(num_tasks)
registerDoParallel(cl)

# Analysis helpers
source("R/VB_bms.R")
source("R/abcrf_helpers.R")
source("R/plot_helpers.R")

simulations <- as.data.table(readRDS("data/Exp1/Simulate_Summary_dt_r05.rds"))
observed    <-  as.data.table(readRDS("data/Exp1/Human_Summary_dt.rds"))

plot_dir_tiff    <- "fig/tiff"
plot_dir_png    <- "fig/png"
default_dpi <- 300
plot_width  <- 8
plot_height <- 5

raw_path <- "data/Exp1/df.csv"
S_path   <- "data/Stengard/s.csv"

df <- load_clean_data(raw_path, type = "experiment")
s  <- load_clean_data(S_path,   type = "stengard")

simulations$model <- factor(simulations$model)

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

pred_oob <- predict(abc_model$model.rf, data = simulations, type = "response")$predictions
true_mod <- simulations$model

mods_in_data <- sort(unique(c(as.character(true_mod), as.character(pred_oob))))
std_names <- clean_model_names(mods_in_data)
names(std_names) <- mods_in_data
fam_std <- assign_family(std_names)
fam_map_full <- fam_std
names(fam_map_full) <- names(std_names)

# Compute family-level errors from those training-data predictions.
res_fam <- family_oob_metrics(pred_oob, true_mod, fam_map_full)
res_fam$confusion
res_fam$err_micro
res_fam$err_macro
res_fam$prior_chance

# Primary experiment model probabilities
posterior_matrix <- compute_posterior(observed, simulations, abc_model)
model_names_from_posterior <- clean_model_names(colnames(posterior_matrix))
posterior_long <- to_long_posterior(posterior_matrix, observed)

df_unique <- df %>%
  distinct(subject, format) %>%
  mutate(subject = as.character(subject))

posterior_long <- posterior_long %>%
  mutate(id = as.character(id)) %>%
  left_join(df_unique, by = c("id" = "subject"))

saveRDS(posterior_long, "data/Exp1/posterior_long.rds")

# Primary experiment family and model BMS
# Divide model probabilities by the uniform model prior, add eps, and take logs.
# Run two-stage BMS separately for probability and frequency formats

post_p <- build_matrix(posterior_long, "probability")
post_f <- build_matrix(posterior_long, "frequency")
K_models <- ncol(posterior_matrix)
eps      <- 1e-10
m_p <- log((post_p /  (1 / K_models)) + eps)
m_f <- log((post_f / (1 / K_models)) + eps)
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

# experiment 1 BMS figures and tables
print_family_table(out_p, title = "Prob-based Family BMS")
print_model_table(out_p,  title = "Prob-based Model BMS")

print_family_table(out_f, title = "Freq-based Family BMS")
print_model_table(out_f,  title = "Freq-based Model BMS")

pa <- plot_two_stage_bms(out_p, title = "Probability Format")
pb <- plot_two_stage_bms(out_f, title = "Frequency Format")
final_plot <- (pa + pb) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Figure_13.tiff"), final_plot,
       width = 8, height = 5, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_13.png"), final_plot,
       width = 8, height = 5, dpi = 300)

# Stengard analysis
simulations_s <- as.data.table(readRDS("data/Stengard/Simulate_Summary_dts_r05.rds"))
observed_s    <-  as.data.table(readRDS("data/Stengard/Human_Summary_dts.rds"))

df_s <- read_csv("data/Stengard/s.csv") |>
  transmute(
    subject, format,trial,
    true_posterior,response,
    BR = br,
    HR = hr,
    FAR = far,
    subject_s = paste0(subject, "_", ifelse(format == "frequency", 1, 2))
  )

simulations_s$model <- factor(simulations_s$model)

# Stengard ABC-RF classifier
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

# Stengard classifier predictions on the training data
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

# Stengard model probabilities
posterior_matrix_s <- compute_posterior(observed_s, simulations_s, abc_model_s)
model_names_from_posterior_s <- clean_model_names(colnames(posterior_matrix_s))
posterior_long_s <- to_long_posterior(posterior_matrix_s, observed_s)

df_unique_s <- df_s %>%
  distinct(subject_s, format) %>%
  mutate(subject = as.character(subject_s))

posterior_long_s <- posterior_long_s %>%
  mutate(id = as.character(id)) %>%
  left_join(df_unique_s, by = c("id" = "subject_s"))

saveRDS(posterior_long_s, "data/Stengard/posterior_long_s.rds")

# family and model BM
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

# Stengard BMS figures and tables
# Save Figure 6 in PNG and TIFF formats; print BMS tables.

print_family_table(out_s_p, title = "Prob-based Family BMS")
print_model_table(out_s_p,  title = "Prob-based Model BMS")
print_family_table(out_s_f, title = "Freq-based Family BMS")
print_model_table(out_s_f,  title = "Freq-based Model BMS")

c <- plot_two_stage_bms(out_s_p, title = "Probability Format")
d <- plot_two_stage_bms(out_s_f, title = "Frequency Format")
final_plot_s <- (c + d) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Figure_6.tiff"), final_plot_s,
       width = 8, height = 5, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_6.png"), final_plot_s,
       width = 8, height = 5, dpi = 300)

# Saved outputs:
# posterior_long.rds and posterior_long_s.rds
# Figures 13 (experiment 1) and 6 (Stengard), in PNG and TIFF formats
# BMS tables are printed, not saved as LaTeX files by this script.
# Note: Close parallel cluster with parallel::stopCluster(cl) when done

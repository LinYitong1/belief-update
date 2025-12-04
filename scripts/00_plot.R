# ================================================================
# Belief Updating Analysis Script (Figure)
# Author: Yitong Lin
# ================================================================

# 0. Setup -------------------------------------------------------------------
# Load necessary packages
pacman::p_load(
  readxl, readr, dplyr, tidyr, ggplot2, scales, viridisLite, patchwork,
  lme4, broom, gt, data.table, PairedData, gridExtra, RColorBrewer,ggpubr
)
options(stringsAsFactors = FALSE)

# Define output directory and plotting defaults
plot_dir_tiff    <- "fig/tiff"
plot_dir_png    <- "fig/png"
default_dpi <- 300
plot_width  <- 6
plot_height <- 4.5

# Source model and utility functions
source("R/models.R")             
source("R/plot_helpers.R") 
source("R/metrics.R")

# 1. Data Import ------------------------------------------------------------
raw_path <- "data/df.csv"
S_path   <- "data/s.csv"
Si_path  <- "data/sirota.csv"

df <- load_clean_data(raw_path, type = "experiment")
s  <- load_clean_data(S_path,   type = "stengard")
si <- load_clean_data(Si_path,  type = "sirota")

dt <- as.data.table(readRDS("data/prediction_dt.rds"))
dt_s<-as.data.table(readRDS("data/prediction_dt_s.rds"))

# ========== Figure 1: Distribution of participants’ estimated posterior probabilities (Sirota et al. 2024). ==========
# Compute heuristic predictions and reference lines
si_m <- compute_heuristics(si)
vlines <- extract_reference_lines(si_m) 

# Plot human response distribution and overlay heuristic anchors
p1 <- ggplot(si_m, aes(x = response_pct)) +
  geom_histogram(binwidth = 1.5, fill = "#5b5e6e") +
  geom_vline(data = vlines, aes(xintercept = value, colour = heuristic),
             linetype = "dashed", linewidth = 1.5) +
  facet_wrap(~ format, ncol = 1, scales = "free_y") +
  scale_x_continuous("Estimated Probability (%)", limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_y_continuous("N", limits = c(0, 60), breaks = seq(0, 60, 20)) +
  scale_colour_manual(values = okabe_ito, name = NULL,
                      guide = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold"),
    legend.position    = "bottom",
    legend.key.height  = unit(4, "mm"),
    legend.text        = element_text(size = 12)
  )


ggsave(file.path(plot_dir_tiff, "Figure_1.tiff"), p1,
       width = plot_width, height = plot_height, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_1.png"), p1,
       width = plot_width, height = plot_height, dpi = 300)

# ========== Figure 2: Response distributions and model predictions (Stengard)  ==========
# Define example tasks for simulation
tasks <- list(
  list(BR = 0.1, HR = 0.5, FAR = 0.3),
  list(BR = 0.7, HR = 0.9, FAR = 0.5),
  list(BR = 0.9, HR = 0.9, FAR = 0.1)
)

# Compute heuristics and task IDs
S_clean <- compute_heuristics(s) %>%
  mutate(task_id = create_task_id(BR, HR, FAR))

# Simulate model predictions for each task
param_all <- do.call(rbind, lapply(tasks, function(tk) {
  tt <- s %>% filter(BR == tk$BR, HR == tk$HR, FAR == tk$FAR) %>% slice(1)
  simulate_bs_models(tt)
}))

# Compute and reshape Bayesian Sampler predictions
pv_all <- predict_bs_models(param_all)

# Generate comparison panel plots for each task
p_panels <- lapply(seq_along(tasks), function(i) {
  plot_task_panel(paste0(
    "BR=", tasks[[i]]$BR,
    ", HR=", tasks[[i]]$HR,
    ", FAR=", tasks[[i]]$FAR
  ), include_legend = (i==1))
})

# Combine and save Figure 2
p_combined <- wrap_plots(p_panels) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.title      = element_text(size = 5, hjust = 0.5)
  )

ggsave(file.path(plot_dir_tiff, "Figure_2.tiff"), p_combined,
       width = plot_width, height = plot_height, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_2.png"), p_combined,
       width = plot_width, height = plot_height, dpi = 300)
# ========== Figure 3: Mean Response vs True Posterior (Stengård) ==========
plot_mean_vs_tp(
  df = s,
  group_vars = c("format", "BR", "HR", "FAR"),
  color_var  = "format",
  legend_title = "Format",
  file_stub = "Figure_3"
)

# ========== Figure 4: Mean Response vs True Posterior for non match (Stengård) ==========
s_nm <- filter_nonmatches(s)

# Plot non-matching trials
plot_mean_vs_tp(
  df = s_nm,
  group_vars = c("format", "BR", "HR", "FAR"),
  color_var  = "format",
  legend_title = "Format",
  file_stub = "Figure_4"
)

# ========== Figure 5: Model Comparison (accuracy) for Stengard  ==========
# ----- Compute Accuracy -----
# Accuracy = % of predictions within 3% of the true posterior
A_Simulate_s <- compute_CLC_summary(dt_s)
plot_model_accuracy(A_Simulate_s, file_stub = "Figure_5")

# ========== Figure 6: Group-level model selection results for Stengard ==========
# NOTE: This figure is generated by the script '04_abcrf_individual_inference.R 
# SECTION 13: Generate Figures (Stengard)

# ========= Figure 7: MH Posterior predictive distributions by format and condition(Stengard)=========
# NOTE: This figure is generated by the script '05_posterior_MH.R 
# SECTION 3: POSTERIOR INFERENCE - STENGARD REPLICATION DATASET

# ========= Figure 8: HABS Posterior predictive distributions by format and condition(Stengard)=========
# NOTE: This figure is generated by the script '05_posterior.R 
# SECTION 3: POSTERIOR INFERENCE - STENGARD REPLICATION DATASET

# ========= Figure 9: Mean Response vs True Posterior for Experiment Data ==========
plot_mean_vs_tp(
  df = df,
  group_vars = c("format", "BR", "HR", "FAR"),
  color_var  = "format",
  legend_title = "Format",
  file_stub = "Figure_9"
)

# ========== Figure 10: Mean Response vs True Posterior for non match (Experiment) ==========
df_nm <- filter_nonmatches(df)

# Plot non-matching trials
plot_mean_vs_tp(
  df = df_nm,
  group_vars = c("format", "BR", "HR", "FAR"),
  color_var  = "format",
  legend_title = "Format",
  file_stub = "Figure_10"
)

# ========== Figure 11: Variance of response within each subject's trial set (Experiment) ==========

human_dt <- as.data.table(df)
V_Human <- compute_variance_summary(
  dt = human_dt,
  column = "response",
  group_vars = c("subject", "BR", "HR", "FAR","format"),
  summary_vars = c("subject","format")
)
p_V_Human<-ggplot(V_Human, aes(x = mean_variance)) +
  geom_histogram(binwidth = 0.002, fill = "#4575b4", color = "black") +
  labs(
    x = "Mean Variance Across Repeated Items",
    y = "Number of Participants") +
  theme_minimal()

ggsave(file.path(plot_dir_tiff, "Figure_11.tiff"), p_V_Human,
       width = plot_width, height = plot_height, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_11.png"), p_V_Human,
       width = plot_width, height = plot_height, dpi = 300)

# ========== Figure 12: Model Comparison (accuracy) for Experiment  ==========
# ----- Compute Accuracy -----
# Accuracy = % of predictions within 3% of the true posterior
A_Simulate <- compute_CLC_summary(dt)
plot_model_accuracy(A_Simulate, file_stub = "Figure_12")

# ========== Figure 13: Group-level model selection results for Experiment ==========
# NOTE: This figure is generated by the script '04_abcrf_individual_inference.R 
# SECTION 8: Generate Figures (Experiment)

# ========= Figure 14: HABS Posterior predictive distributions by format and condition(Experiment)=========
# NOTE: This figure is generated by the script '05_posterior.R 
# SECTION 8: POSTERIOR INFERENCE - STENGARD REPLICATION DATASET


# ========= Appendix D:Mean and Accuracy=========
df <- df %>%
  mutate(
    se = (response - true_posterior)^2,
    e = abs(response - true_posterior)
  )

# Panel A: Mean Judgments
pd_mean <- compute_paired_profiles(df, "response")
p1 <- plot_paired_profile(pd_mean, "A. Mean Judgments by Format", "Mean Judgment")

# Panel B: Absolute Deviations
pd_dev <- compute_paired_profiles(df, "e")
p2 <- plot_paired_profile(pd_dev, "B. Absolute Deviations by Format", "Mean Absolute Deviation")

combined <- ggarrange(p1, p2, ncol = 2, nrow = 1, align = "hv")
ggsave(file.path(plot_dir_tiff, "Apendix_D.tiff"),combined,
       width = plot_width, height = plot_height, dpi = 600)
ggsave(file.path(plot_dir_png, "Apendix_D.png"), combined,
       width = plot_width, height = plot_height, dpi = 300)

# ========= Appendix E1: Posterior distributions of Mixture Heuristic model parameters (Stengard).=========
# ========= Appendix E2: Posterior distributions of HM_Mixture model parameters (Stengard).=========
# ========= Appendix E3: Posterior distributions of HM_Mixture model parameters (Expeirment).=========

# ========= Appendix F: Posterior predictive p-value distributions=========
install.packages("magick")
library(cowplot)
library(magick)
p_s_mh <- ggdraw() + draw_image("fig/png/PPP_hist_Stengard_indiv_MH.png")
p_s_hybrid <- ggdraw() + draw_image("fig/png/PPP_hist_Stengard_indiv.png")
p_e_mh <- ggdraw() + draw_image("fig/png/PPP_hist_exp_indiv_MH.png")
p_e_hybrid<- ggdraw() + draw_image("fig/png/PPP_hist_Experimental_indiv.png")

ppp_combined <- plot_grid(
  p_s_mh,p_s_hybrid,
  p_e_mh, p_e_hybrid,
  ncol = 2,
  labels = c("a", "b", "c", "d")  
)

ggsave(file.path(plot_dir_tiff, "Apendix_F.tiff"),ppp_combined,
       width = plot_width, height = plot_height, dpi = 600)
ggsave(file.path(plot_dir_png, "Apendix_F.png"), ppp_combined,
       width = plot_width, height = plot_height, dpi = 300)

# ======================== End of Script ============================

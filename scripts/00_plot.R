# ================================================================
# Belief Updating Analysis Script (Publication Ready)
# Author: Yitong Lin
# Date: 2025-06-09
# ================================================================

# 0. Setup -------------------------------------------------------------------
# Load necessary packages
pacman::p_load(
  readxl, readr, dplyr, tidyr, ggplot2, scales, viridisLite, patchwork,
  lme4, broom, gt, data.table, PairedData, gridExtra, RColorBrewer,ggpubr
)
options(stringsAsFactors = FALSE)

# Define output directory and plotting defaults
plot_dir    <- "fig"
default_dpi <- 300
plot_width  <- 6
plot_height <- 4.5

# Source model and utility functions
source("R/models.R")              # Contains BS, BS_R, simulate_relative, etc.
source("R/plot_helpers.R") # Contains helpers like compute_heuristics, plot_task_panel

# 1. Data Import ------------------------------------------------------------
raw_path <- "data/df.xlsx"
S_path   <- "data/s.csv"
Si_path  <- "data/sirota.csv"

df <- load_clean_data(raw_path, type = "experiment")
s  <- load_clean_data(S_path,   type = "stengard")
si <- load_clean_data(Si_path,  type = "sirota")

# ========== Figure 1: Overall Estimate Distribution with Heuristic Anchors ==========
# Compute heuristic predictions and reference lines
si_m <- compute_heuristics(si)
vlines <- extract_reference_lines(si_m) 

# Plot human response distribution and overlay heuristic anchors
p1 <- ggplot(si_m, aes(x = response_pct)) +
  geom_histogram(binwidth = 1.5, fill = "#5b5e6e") +
  geom_vline(data = vlines, aes(xintercept = value, colour = heuristic),
             linetype = "dashed", linewidth = .9) +
  facet_wrap(~ format, ncol = 1, scales = "free_y") +
  scale_x_continuous("Estimated Probability (%)", limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_y_continuous("N", limits = c(0, 60), breaks = seq(0, 60, 20)) +
  scale_colour_manual(values = okabe_ito, name = NULL,
                      guide = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_line(linewidth = .3, colour = "grey85"),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold"),
    legend.position    = "bottom",
    legend.key.height  = unit(4, "mm"),
    legend.text        = element_text(size = 8)
  )

# Save Figure 1
ggsave(file.path(plot_dir, "Figure_1.png"), p1,
       width = plot_width, height = plot_height, dpi = default_dpi)

# ========== Figure 2: Human vs Bayesian Sampler Model, Three Example Tasks ==========
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
    plot.title      = element_text(size = 6.5, hjust = 0.5)
  )

ggsave(file.path(plot_dir, "Figure_2.png"),
       p_combined, width = 12, height = 4.5, dpi = default_dpi)

# ========== Figure 3: Mean Response vs True Posterior (Stengård) ==========
plot_mean_vs_tp(
  df = s,
  group_vars = c("format", "BR", "HR", "FAR"),
  color_var  = "format",
  legend_title = "Format",
  file_name = "Figure_3.png"
)

# ========== Figure 4: Mean Response vs True Posterior (Non-matching Trials Only) ==========
# Filter responses that do not match any heuristic within tolerance
s_nm <- filter_nonmatches(s)

# Plot non-matching trials
plot_mean_vs_tp(
  df = s_nm,
  group_vars = c("format", "BR", "HR", "FAR"),
  color_var  = "format",
  legend_title = "Format",
  file_name = "Figure_4.png"
)

# ========== Figure 5: Model Accuracy Comparison for Stengard Data in 03  ==========
# ========== Figure 6: Model Fitting Comparison for Stengard Data in 04  ==========
# ========= Figure 7: Mean Response and Accuracy for Experiment Data ==========

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

ggsave(file.path(plot_dir, "Figure_7.png"), combined,
       width = plot_width, height = plot_height, dpi = default_dpi)
# ========== Figure 8: Mean Response vs True Posterior (Experiment)  ==========
plot_mean_vs_tp(
  df = df,
  group_vars = c("format", "BR", "HR", "FAR"),
  color_var  = "format",
  legend_title = "Format",
  file_name = "Figure_8.png"
)
# ========== Figure 9: Mean Response vs True Posterior (Non-matching Trials Only)  ==========
# Filter responses that do not match any heuristic within tolerance
df_nm <- filter_nonmatches(df)

# Plot non-matching trials
plot_mean_vs_tp(
  df = df_nm,
  group_vars = c("format", "BR", "HR", "FAR"),
  color_var  = "format",
  legend_title = "Format",
  file_name = "Figure_9.png"
)
# ========== Figure 10: Model Accuracy Comparison for Experiment Data in 03  ==========
# ========== Figure 11: Model Fitting Comparison for Experiment Data in 04  ==========

# ======================== End of Script ============================

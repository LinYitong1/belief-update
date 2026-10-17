# ABC parameter inference and posterior predictive validation

# 1. Setup and configuration
set.seed(20261011)
# 1.1 Load required packages
load_packages <- function() {
  required <- c(
    "data.table", "tidyverse", "abc", "future", "furrr", "progressr",
    "ggplot2", "patchwork", "readxl", "viridisLite", "ggforce"
  )
  purrr::walk(required, ~{
    if (!require(.x, character.only = TRUE)) install.packages(.x)
    library(.x, character.only = TRUE)
  })
}

load_packages()

# 1.2 Source analysis helper scripts
source("R/metrics.R")
source("R/models.R")
source("R/plot_helpers.R")
source("R/table_helper.R")
# Posterior helper functions support subject- and format-level inference/PPC.
source("R/posterior_helper.R")
`.` <- function(...) list(...)
# 1.3 Figure output settings
plot_dir_tiff <- "fig/tiff"
plot_dir_png  <- "fig/png"
default_dpi   <- 300
plot_width    <- 6
plot_height   <- 4.5

if (!dir.exists(plot_dir_tiff)) dir.create(plot_dir_tiff, recursive = TRUE)
if (!dir.exists(plot_dir_png))  dir.create(plot_dir_png,  recursive = TRUE)

# 2. Load data and define parameter labels
# 2.1 Observed data: Experiment 1 and Stengård
human_dt   <- as.data.table(load_clean_data("data/Exp1/df.csv", type = "experiment"))
human_summary<-as.data.table(readRDS("data/Stengard/Human_Summary_dts.rds"))
human_data_format <- human_dt %>%
  distinct(subject, format)%>%
  dplyr::mutate(
    subject = as.double(subject),
    format  = as.character(format)
  )

human_dt_s <- as.data.table(load_clean_data("data/Stengard/s.csv",  type = "stengard"))
human_data_s_format <- human_dt_s %>%
  distinct(subject, format)%>%
  dplyr::mutate(
    subject = as.character(subject),
    format  = as.character(format)
  )
# 2.2 Simulated summaries and parameters: Experiment 1 (r05)
simulations    <- as.data.table(readRDS("data/Exp1/Simulate_Summary_dt_r05.rds")) %>%
  dplyr::filter(model == "M_L")
parameter_dt   <- readRDS("data/Exp1/parameter_dt_r05.rds")

# 2.3 Simulated summaries and parameters: Stengård (r05)
simulations_s  <- as.data.table(readRDS("data/Stengard/Simulate_Summary_dts_r05.rds")) %>%
  dplyr::filter(model == "M_L")
parameter_dt_s <- readRDS("data/Stengard/parameter_dt_s_r05.rds")

# 2.4 Saved model-comparison posterior data
posterior_long   <- readRDS("data/Exp1/posterior_long.rds")
posterior_long_s <- readRDS("data/Stengard/posterior_long_s.rds")

# 2.5 Model parameters and display labels
all_params <- c(paste0("p8_", 1:3), paste0("q_", 1:6), "N9", "v9")

old_names <- all_params

label_values <- c(
  "P(Integrated)", "P(Heuristic)", "P(Random Noise)",
  "P(BO)", "P(HO)", "P(FO)",
  "P(JO)", "P(LS)", "P(50%)",
  "Sample_Size", "Prior_Size"
)
label_map <- setNames(label_values, old_names)

# 3. Stengård: participant-level ABC and posterior predictive checks
# 3.1 Estimate participant-level posteriors and compute PPC summaries
stengard_subject_results <- run_subject_level_validation(
  dataset_name       = "Stengard (individual)",
  human_data         = human_dt_s,
  sim_params         = parameter_dt_s,
  simulations        = simulations_s,
  subject_col        = "subject",
  abc_tol            = 0.05,
  bins_overlap       = 30,
  ppp_alpha          = 0.05,
  calculate_variance = FALSE
)

abc_indiv_s       <- stengard_subject_results$abc_indiv
observed_indiv_s  <- stengard_subject_results$observed_indiv
post_pred_indiv_s <- stengard_subject_results$post_pred_indiv
posterior_indiv_s <- stengard_subject_results$posterior_indiv
stats_long_indiv_s <- stengard_subject_results$stats_long_indiv
ppp_by_subject_s   <- stengard_subject_results$ppp_by_subject
ppp_summary_s      <- stengard_subject_results$ppp_summary_overall
print(ppp_summary_s)

# 3.2 Plot participant-level posterior predictive p-values
alpha <- 0.05

ppp_hist_s <- ggplot(ppp_by_subject_s, aes(x = ppp)) +
  geom_histogram(
    bins     = 30,
    color    = "white",
    fill     = "#3C5488FF",
    alpha    = 0.9,
    position = "identity"
  ) +
  geom_vline(
    xintercept = alpha,
    linetype   = "dashed"
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "B",
    x     = "Posterior predictive p-value (PPP)",
    y     = "Count"
  ) +
  theme_classic(base_size = 10)

ggsave(file.path(plot_dir_tiff, "PPP_hist_Stengard_indiv.tiff"), ppp_hist_s,
       width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPP_hist_Stengard_indiv.png"),  ppp_hist_s,
       width = 6, height = 4, dpi = 300)

# 3.3 Compare observed and posterior predictive response distributions
ref_points_s <- tibble(
  BR  = c(0.1, 0.7, 0.9),
  HR  = c(0.5, 0.9, 0.9),
  FAR = c(0.3, 0.5, 0.1)
)

ref_lines_s <- get_reference_lines(human_dt_s, ref_points_s)

human_df_s_prob <- human_dt_s %>%
  filter(format == "probability") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  filter(condition %in% unique(ref_lines_s$condition))

human_df_s_freq <- human_dt_s %>%
  filter(format == "frequency") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  filter(condition %in% unique(ref_lines_s$condition))

post_pred_indiv_s <- post_pred_indiv_s %>%
  left_join(human_data_s_format, by = "subject")

model_df_prob_s <- post_pred_indiv_s %>%
  filter(format == "probability") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  filter(condition %in% unique(ref_lines_s$condition))

model_df_freq_s <- post_pred_indiv_s %>%
  filter(format == "frequency") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  filter(condition %in% unique(ref_lines_s$condition))

plot_prob_overlay_s <- plot_overlay(
  human_df  = human_df_s_prob,
  model_df  = model_df_prob_s,
  ref_lines = ref_lines_s,
  format_label = "Probability Format"
)

plot_freq_overlay_s <- plot_overlay(
  human_df  = human_df_s_freq,
  model_df  = model_df_freq_s,
  ref_lines = ref_lines_s,
  format_label = "Frequency Format"
)

final_plot_s <- plot_prob_overlay_s + plot_freq_overlay_s +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Figure_8.tiff"),
       final_plot_s, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "Figure_8.png"),
       final_plot_s, width = 6, height = 4, dpi = 300)

# 4. Experiment 1: participant-level ABC and posterior predictive checks

# 4.1 Estimate participant-level posteriors and compute PPC summaries
experimental_subject_results <- run_subject_level_validation(
  dataset_name       = "Experimental (individual)",
  human_data         = human_dt,
  sim_params         = parameter_dt,
  simulations        = simulations,
  subject_col        = "subject",
  abc_tol            = 0.05,
  bins_overlap       = 30,
  ppp_alpha          = 0.05,
  calculate_variance = TRUE
)

abc_indiv_exp       <- experimental_subject_results$abc_indiv
observed_indiv_exp  <- experimental_subject_results$observed_indiv
post_pred_indiv_exp <- experimental_subject_results$post_pred_indiv
posterior_indiv_exp <- experimental_subject_results$posterior_indiv
stats_long_indiv_exp <- experimental_subject_results$stats_long_indiv
ppp_by_subject_exp   <- experimental_subject_results$ppp_by_subject
ppp_summary_exp      <- experimental_subject_results$ppp_summary_overall
print(ppp_summary_exp)
# 4.2 Plot participant-level posterior predictive p-values
ppp_hist_exp <-ggplot(ppp_by_subject_exp, aes(x = ppp)) +
  geom_histogram(
    bins     = 30,
    color    = "white",
    fill = "#3C5488FF",
    alpha    = 0.9,
    position = "identity"
  ) +
  geom_vline(
    xintercept = alpha,
    linetype   = "dashed"
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "DS",
    x     = "Posterior predictive p-value (PPP)",
    y     = "Count",
    fill  = ""
  ) +
  theme_classic(base_size = 10)

ggsave(file.path(plot_dir_tiff, "PPP_hist_Experimental_indiv.tiff"),
       ppp_hist_exp, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPP_hist_Experimental_indiv.png"),
       ppp_hist_exp, width = 6, height = 4, dpi = 300)

# 4.3 Compare observed and posterior predictive response distributions
ref_points <- tibble(BR = c(0.01,  0.4,  0.97),
                     HR = c( 0.72 , 0.53 ,0.53 ),
                     FAR = c(0.42,0.31, 0.11))

ref_lines <- get_reference_lines(human_dt, ref_points)

human_df_prob <- human_dt %>%
  filter(format == "probability") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  filter(condition %in% unique(ref_lines$condition))   %>%
  dplyr::mutate(
    subject = as.double(subject),
    format  = as.character(format)
  )

human_df_freq <- human_dt %>%
  filter(format == "frequency") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  filter(condition %in% unique(ref_lines$condition))

post_pred_indiv_exp <- post_pred_indiv_exp %>%
  left_join(human_data_format, by = "subject")

model_df_prob_exp <- post_pred_indiv_exp %>%
  filter(format == "probability") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  filter(condition %in% unique(ref_lines$condition))%>%
  dplyr::mutate(
    subject = as.double(subject),
    format  = as.character(format)
  )

model_df_freq_exp <- post_pred_indiv_exp %>%
  filter(format == "frequency") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  filter(condition %in% unique(ref_lines$condition))

plot_prob_overlay <- plot_overlay(
  human_df  = human_df_prob,
  model_df  = model_df_prob_exp,
  ref_lines = ref_lines,
  format_label = "Probability Format"
)

plot_freq_overlay <- plot_overlay(
  human_df  = human_df_freq,
  model_df  = model_df_freq_exp,
  ref_lines = ref_lines,
  format_label = "Frequency Format"
)

final_plot <- plot_prob_overlay + plot_freq_overlay +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

final_plot

ggsave(file.path(plot_dir_tiff, "Exp1_ML_Predictive_Distributions.tiff"),
       final_plot, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "Exp1_ML_Predictive_Distributions.png"),
       final_plot, width = 6, height = 4, dpi = 300)

# 5. Format-level ABC parameter inference and prior/posterior plots
# 5.1 Experiment 1
# 5.1.1 Compute observed summary statistics by format
human_sum_fmt_exp <- compute_all_metrics(
  df         = human_dt,
  column     = "response",
  group_vars = c("format")
)

slope_int_fmt_exp <- compute_SI_by(
  dt         = human_dt,
  group_vars = c("format"),
  predictors = c("BR", "HR", "FAR"),
  column     = "response"
)

variance_int_fmt_exp <- compute_variance_summary(
  data.table::as.data.table(human_dt),
  column       = "response",
  group_vars   = c("format", "BR", "HR", "FAR"),
  summary_vars = c("format")
)

human_summary_fmt_exp <- human_sum_fmt_exp %>%
  dplyr::left_join(slope_int_fmt_exp,   by = "format") %>%
  dplyr::left_join(variance_int_fmt_exp, by = "format")

# 5.1.2 Infer parameter posteriors and plot prior/posterior distributions
abc_group <- compute_abc_by_format(
  human_summary = human_summary_fmt_exp,
  sim_params    = parameter_dt,
  simulations   = simulations,
  tol           = 0.05
)

posterior_means <- lapply(names(abc_group$post_draws), function(fmt) {
  abc_group$post_draws[[fmt]] %>%
    dplyr::summarise(dplyr::across(all_of(all_params), ~ mean(.x, na.rm = TRUE))) %>%
    dplyr::mutate(format = fmt)
}) %>%
  dplyr::bind_rows()

plot_list_group_exp <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_group(
    param_name = param,
    sim_params = parameter_dt,
    abc_fmt    =abc_group,
    x_label    = label_map[[param]]
  )
})

final_plot_group_exp <- patchwork::wrap_plots(plot_list_group_exp, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Appendix_E4.tiff"),
       final_plot_group_exp, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Appendix_E4.png"),
       final_plot_group_exp, width = 8, height = 8, dpi = 300)

# 5.2 Stengård: compute observed summary statistics by format
human_sum_fmt_s <- compute_all_metrics(
  df         = human_dt_s,
  column     = "response",
  group_vars = c("format")
)

slope_int_fmt_s <- compute_SI_by(
  dt         = human_dt_s,
  group_vars = c("format"),
  predictors = c("BR", "HR", "FAR"),
  column     = "response"
)

human_summary_fmt_s <- human_sum_fmt_s %>%
  dplyr::left_join(slope_int_fmt_s,   by = "format")

# 5.2.1 Infer parameter posteriors and plot prior/posterior distributions
abc_group_s <- compute_abc_by_format(
  human_summary = human_summary_fmt_s,
  sim_params    = parameter_dt_s,
  simulations   = simulations_s,
  tol           = 0.05
)

posterior_means_s <- lapply(names(abc_group_s$post_draws), function(fmt) {
  abc_group_s$post_draws[[fmt]] %>%
    dplyr::summarise(dplyr::across(all_of(all_params), ~ mean(.x, na.rm = TRUE))) %>%
    dplyr::mutate(format = fmt)
}) %>%
  dplyr::bind_rows()

plot_list_group_s <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_group (
    param_name = param,
    sim_params = parameter_dt_s,
    abc_fmt    =abc_group_s,
    x_label    = label_map[[param]]
  )
})
plot_list_group_s <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_group(
    param_name = param,
    sim_params = parameter_dt_s,
    abc_fmt    = abc_group_s,
    x_label    = label_map[[param]]
  )
})

final_plot_group_s <- patchwork::wrap_plots(plot_list_group_s, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Appendix_E2.tiff"),
       final_plot_group_s, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Appendix_E2.png"),
       final_plot_group_s, width = 8, height = 8, dpi = 300)

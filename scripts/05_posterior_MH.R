################################################################################
# 1. ENVIRONMENT SETUP
################################################################################

# ---- 1.1 Package management -------------------------------------------------
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

# ---- 1.2 Source custom helper functions ------------------------------------
source("R/metrics.R")
source("R/models.R")
source("R/plot_helpers.R")
source("R/table_helper.R")
# IMPORTANT: this is your new SUBJECT-LEVEL ONLY helper
source("R/posterior_helper.R")
`.` <- function(...) list(...)
# ---- 1.3 Output configuration ----------------------------------------------
plot_dir_tiff <- "fig/tiff"
plot_dir_png  <- "fig/png"
default_dpi   <- 300
plot_width    <- 6
plot_height   <- 4.5
out_path  <- "Table"
if (!dir.exists(plot_dir_tiff)) dir.create(plot_dir_tiff, recursive = TRUE)
if (!dir.exists(plot_dir_png))  dir.create(plot_dir_png,  recursive = TRUE)


################################################################################
# 2. DATA LOADING
################################################################################

# ---- 2.1 Empirical data ----------------------------------------------------
# Trial-level responses for both datasets
human_dt   <- as.data.table(load_clean_data("data/df.csv", type = "experiment"))
human_summary<-as.data.table(readRDS("data/Human_Summary_dts.rds"))
human_data_format <- human_dt %>%
  distinct(subject, format)%>%
  dplyr::mutate(
    subject = as.double(subject), 
    format  = as.character(format)  
  )


human_dt_s <- as.data.table(load_clean_data("data/s.csv",  type = "stengard"))


human_data_s_format <- human_dt_s %>%
  distinct(subject, format)%>%
  dplyr::mutate(
    subject = as.character(subject), 
    format  = as.character(format)  
  )
# ---- 2.2 Simulation data & priors: Experimental ----------------------------
simulations_MH    <- as.data.table(readRDS("data/Simulate_Summary_dt.rds")) %>%
  dplyr::filter(model == "MH")
parameter_dt   <- readRDS("data/parameter_dt.rds")

# ---- 2.3 Simulation data & priors: Stengard --------------------------------
simulations_s_MH   <- as.data.table(readRDS("data/Simulate_Summary_dts.rds")) %>%
  dplyr::filter(model == "MH")
parameter_dt_s <- readRDS("data/parameter_dt_s.rds")

# ---- 2.4 Model comparison results -------
posterior_long   <- readRDS("data/posterior_long.rds")
posterior_long_s <- readRDS("data/posterior_long_s.rds")

# ---- 2.5 Parameter specification -------------------------------------------
all_params <- c(paste0("p7_", 1:7))

old_names <- all_params

label_values <- c("P(BO)", "P(HO)", "P(FO)",
                  "P(JO)", "P(LS)", "P(50%)", "P(Random)")
label_map <- setNames(label_values, old_names)


################################################################################
# 3. Non Match Trial analysis for MH participants
################################################################################

best_fit_s <- posterior_long_s %>%
  group_by(id, format) %>%            
  slice_max(p, n = 1, with_ties = FALSE) %>%
  ungroup()

mh_subj_s <- best_fit_s %>%
  filter(model == "MH") %>%
  distinct(id) %>%                     
  rename(subject = id)   


df_s_mh <- human_dt_s %>%
  inner_join(mh_subj_s, by = "subject")


s_mh_nm <- filter_nonmatches(df_s_mh)


filtered_s_mh <- s_mh_nm %>%
  group_by(subject) %>%
  filter(n() >= 3, n_distinct(true_posterior) > 1) %>%
  ungroup()


slope_s_mh_non <- slope_t_summary(
  filtered_s_mh,
  mu       = 1, 
  out_path = file.path(out_path, "slope_s_mh_non.tex"),
  caption  = "Slope summary of non exact match for MH Participants in Stengard."
)

cat("\n=== Table: Slope summary of non exact match for MH Participants in Stengard ===\n")
print(slope_s_mh_non)

l_df_n <- lmer(
  response ~ true_posterior * format + (1 | subject),
  data = filtered_s_mh
)
summary(l_df_n)

################################################################################
# 4. SUBJECT-LEVEL Posterior Check (Stengard, MH)
################################################################################

stengard_subject_MH_results <- run_subject_level_MH_validation(
  dataset_name       = "Stengard (MH, individual)",
  human_data         = human_dt_s,
  sim_params         = parameter_dt_s,
  simulations        = simulations_s_MH,
  subject_col        = "subject",  
  abc_tol            = 0.05,
  ppp_alpha          = 0.05,
  calculate_variance = FALSE
)

# Unpack for convenience
abc_indiv_s_MH        <- stengard_subject_MH_results$abc_indiv
observed_indiv_s_MH   <- stengard_subject_MH_results$observed_indiv
post_pred_indiv_s_MH  <- stengard_subject_MH_results$post_pred_indiv
posterior_indiv_s_MH  <- stengard_subject_MH_results$posterior_indiv
stats_long_indiv_s_MH <- stengard_subject_MH_results$stats_long_indiv
ppp_by_subject_s_MH   <- stengard_subject_MH_results$ppp_by_subject
ppp_summary_s_MH      <- stengard_subject_MH_results$ppp_summary_overall

alpha <- 0.05

# PPP distribution (Stengard, MH)
ppp_hist_s_MH <- ggplot(ppp_by_subject_s_MH, aes(x = ppp)) +
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
    title = "PPP distribution (Stengard, individual level, MH)",
    x     = "Posterior predictive p-value (PPP)",
    y     = "Count"
  ) +
  theme_classic(base_size = 10)

ggsave(file.path(plot_dir_tiff, "PPP_hist_Stengard_indiv_MH.tiff"),
       ppp_hist_s_MH, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPP_hist_Stengard_indiv_MH.png"),
       ppp_hist_s_MH, width = 6, height = 4, dpi = 300)

# ---- 3.4 Pooled human vs posterior prediction (Stengard, MH) --------------

ref_points_s <- tibble(
  BR  = c(0.1, 0.7, 0.9), 
  HR  = c(0.5, 0.9, 0.9), 
  FAR = c(0.3, 0.5, 0.1)
)

ref_lines_s_MH <- get_reference_lines(human_dt_s, ref_points_s)

human_df_s_prob <- human_dt_s %>% 
  dplyr::filter(format == "probability") %>%
  dplyr::mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  dplyr::filter(condition %in% unique(ref_lines_s_MH$condition))

human_df_s_freq <- human_dt_s %>% 
  dplyr::filter(format == "frequency") %>%
  dplyr::mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  dplyr::filter(condition %in% unique(ref_lines_s_MH$condition))

model_df_prob_s_MH <- post_pred_indiv_s_MH %>% 
  dplyr::filter(format == "probability") %>%
  dplyr::mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  dplyr::filter(condition %in% unique(ref_lines_s_MH$condition))

model_df_freq_s_MH <- post_pred_indiv_s_MH %>% 
  dplyr::filter(format == "frequency") %>%
  dplyr::mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  dplyr::filter(condition %in% unique(ref_lines_s_MH$condition))

plot_prob_overlay_s_MH <- plot_overlay(
  human_df    = human_df_s_prob,
  model_df    = model_df_prob_s_MH,
  ref_lines   = ref_lines_s_MH,
  format_label = "Probability Format"
)

plot_freq_overlay_s_MH <- plot_overlay(
  human_df    = human_df_s_freq,
  model_df    = model_df_freq_s_MH,
  ref_lines   = ref_lines_s_MH,
  format_label = "Frequency Format"
)

final_plot_s_MH <- plot_prob_overlay_s_MH + plot_freq_overlay_s_MH +
  patchwork::plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Figure_7.tiff"),
       final_plot_s_MH, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "Figure_7.png"),
       final_plot_s_MH, width = 6, height = 4, dpi = 300)


################################################################################
# 4. Group Level Posterior (Stengard, MH)
################################################################################

human_sum_fmt_s_MH <- compute_all_metrics(
  df         = human_dt_s,
  column     = "response",
  group_vars = c("format")
)

slope_int_fmt_s_MH <- compute_SI_by(
  dt         = human_dt_s,
  group_vars = c("format"),
  predictors = c("BR", "HR", "FAR"),
  column     = "response"
)

human_summary_fmt_s_MH <- human_sum_fmt_s_MH %>%
  dplyr::left_join(slope_int_fmt_s_MH, by = "format")

# ABC by format (Stengard, MH)
abc_group_s_MH <- compute_abc_by_format_MH(
  human_summary = human_summary_fmt_s_MH,
  sim_params    = parameter_dt_s,
  simulations   = simulations_s_MH,
  tol           = 0.05
)

posterior_means_s_MH <- lapply(names(abc_group_s_MH$post_draws), function(fmt) {
  abc_group_s_MH$post_draws[[fmt]] %>%
    dplyr::summarise(dplyr::across(all_of(all_params), ~ mean(.x, na.rm = TRUE))) %>%
    dplyr::mutate(format = fmt)
}) %>%
  dplyr::bind_rows()

plot_list_group_s_MH <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_group(
    param_name = param,
    sim_params = parameter_dt_s,
    abc_fmt    = abc_group_s_MH,
    x_label    = label_map[[param]]
  )
})

final_plot_group_s_MH <- patchwork::wrap_plots(
  plot_list_group_s_MH, ncol = 3, guides = "collect"
) & theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Appendix_E1.tiff"),
       final_plot_group_s_MH, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Appendix_E1.png"),
       final_plot_group_s_MH, width = 8, height = 8, dpi = 300)


################################################################################
# 5. SUBJECT-LEVEL Posterior Check (Experiment, MH)
################################################################################

exp_subject_MH_results <- run_subject_level_MH_validation(
  dataset_name       = "Experiment (MH, individual)",
  human_data         = human_dt,
  sim_params         = parameter_dt,
  simulations        = simulations_MH,
  subject_col        = "subject",  
  abc_tol            = 0.05,
  ppp_alpha          = 0.05,
  calculate_variance = TRUE
)

# Unpack for convenience
abc_indiv_exp_MH        <- exp_subject_MH_results$abc_indiv
observed_indiv_exp_MH   <- exp_subject_MH_results$observed_indiv
post_pred_indiv_exp_MH  <- exp_subject_MH_results$post_pred_indiv
posterior_indiv_exp_MH  <- exp_subject_MH_results$posterior_indiv
stats_long_indiv_exp_MH <- exp_subject_MH_results$stats_long_indiv
ppp_by_subject_exp_MH   <- exp_subject_MH_results$ppp_by_subject
ppp_summary_exp_MH      <- exp_subject_MH_results$ppp_summary_overall


alpha <- 0.05

ppp_hist_exp_MH <- ggplot(ppp_by_subject_exp_MH, aes(x = ppp)) +
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
    title = "PPP distribution (Experiment, individual level, MH)",
    x     = "Posterior predictive p-value (PPP)",
    y     = "Count"
  ) +
  theme_classic(base_size = 10)

ggsave(file.path(plot_dir_tiff, "Figure_14.tiff"),
       ppp_hist_exp_MH, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "Figure_14.png"),
       ppp_hist_exp_MH, width = 6, height = 4, dpi = 300)

# ---- Pooled human vs posterior prediction (Experiment, MH) -----------------

ref_points_exp <- tibble(BR = c(0.01,  0.4,  0.97),
                         HR = c( 0.72 , 0.53 ,0.53 ), 
                         FAR = c(0.42,0.31, 0.11))

ref_lines_exp_MH <- get_reference_lines(human_dt, ref_points_exp)

human_df_exp_prob <- human_dt %>% 
  dplyr::filter(format == "probability") %>%
  dplyr::mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  dplyr::filter(condition %in% unique(ref_lines_exp_MH$condition))   

human_df_exp_freq <- human_dt %>% 
  dplyr::filter(format == "frequency") %>%
  dplyr::mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  dplyr::filter(condition %in% unique(ref_lines_exp_MH$condition))

model_df_prob_exp_MH <- post_pred_indiv_exp_MH %>% 
  dplyr::filter(format == "probability") %>%
  dplyr::mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  dplyr::filter(condition %in% unique(ref_lines_exp_MH$condition))

model_df_freq_exp_MH <- post_pred_indiv_exp_MH %>% 
  dplyr::filter(format == "frequency") %>%
  dplyr::mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  dplyr::filter(condition %in% unique(ref_lines_exp_MH$condition))

plot_prob_overlay_exp_MH <- plot_overlay(
  human_df    = human_df_exp_prob,
  model_df    = model_df_prob_exp_MH,
  ref_lines   = ref_lines_exp_MH,
  format_label = "Probability Format"
)

plot_freq_overlay_exp_MH <- plot_overlay(
  human_df    = human_df_exp_freq,
  model_df    = model_df_freq_exp_MH,
  ref_lines   = ref_lines_exp_MH,
  format_label = "Frequency Format"
)

final_plot_exp_MH <- plot_prob_overlay_exp_MH + plot_freq_overlay_exp_MH +
  patchwork::plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "PPC_Experiment_MH.tiff"),
       final_plot_exp_MH, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPC_Experiment_MH.png"),
       final_plot_exp_MH, width = 6, height = 4, dpi = 300)


################################################################################
# 6. Group Level Posterior (Experiment, MH)
################################################################################

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
  dplyr::left_join(slope_int_fmt_exp,    by = "format") %>%
  dplyr::left_join(variance_int_fmt_exp, by = "format")

# ABC by format (Experiment, MH)
abc_group_exp_MH <- compute_abc_by_format_MH(
  human_summary = human_summary_fmt_exp,
  sim_params    = parameter_dt,
  simulations   = simulations_MH,
  tol           = 0.05
)

posterior_means_exp_MH <- lapply(names(abc_group_exp_MH$post_draws), function(fmt) {
  abc_group_exp_MH$post_draws[[fmt]] %>%
    dplyr::summarise(dplyr::across(all_of(all_params), ~ mean(.x, na.rm = TRUE))) %>%
    dplyr::mutate(format = fmt)
}) %>%
  dplyr::bind_rows()

plot_list_group_exp_MH <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_group(
    param_name = param,
    sim_params = parameter_dt,
    abc_fmt    = abc_group_exp_MH,
    x_label    = label_map[[param]]
  )
})

final_plot_group_exp_MH <- patchwork::wrap_plots(
  plot_list_group_exp_MH, ncol = 3, guides = "collect"
) & theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Appendix_E3.tiff"),
       final_plot_group_exp_MH, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Appendix_E3.png"),
       final_plot_group_exp_MH, width = 8, height = 8, dpi = 300)


################################################################################
# END OF SCRIPT
################################################################################

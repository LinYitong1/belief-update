################################################################################
# INDIVIDUAL-LEVEL ABC POSTERIOR INFERENCE: HIERARCHICAL HEURISTIC-ANALYTICAL
################################################################################
#
# DESCRIPTION:
# Subject-level Bayesian parameter estimation for belief updating models using
# Approximate Bayesian Computation (ABC).
#
# ANALYSIS WORKFLOW (INDIVIDUAL-LEVEL ONLY):
# 1. Environment setup: packages, helper functions
# 2. Data loading (Stengard replication + Experimental data)
# 3. Subject-level ABC + PPC for Stengard:
#      - ABC per subject
#      - Posterior predictive checks (PPP, overlap)
#      - Population-level prior vs posterior (subject means)
#      - Pooled human vs posterior prediction
# 4. Subject-level ABC + PPC for Experimental:
#      - Same steps as above
#
#
# Author: Yitong Lin 
################################################################################


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
simulations    <- as.data.table(readRDS("data/Simulate_Summary_dt.rds")) %>%
  dplyr::filter(model == "M_L")
parameter_dt   <- readRDS("data/parameter_dt.rds")

# ---- 2.3 Simulation data & priors: Stengard --------------------------------
simulations_s  <- as.data.table(readRDS("data/Simulate_Summary_dts.rds")) %>%
  dplyr::filter(model == "M_L")
parameter_dt_s <- readRDS("data/parameter_dt_s.rds")

# ---- 2.4 Model comparison results (optional, not used directly here) -------
posterior_long   <- readRDS("data/posterior_long.rds")
posterior_long_s <- readRDS("data/posterior_long_s.rds")

# ---- 2.5 Parameter specification -------------------------------------------
all_params <- c(paste0("p8_", 1:3), paste0("q_", 1:6), "N9", "v9")

old_names <- all_params

label_values <- c(
  "P(Integrated)", "P(Heuristic)", "P(Random Noise)",
  "P(BO)", "P(HO)", "P(FO)",
  "P(JO)", "P(LS)", "P(50%)",
  "Sample_Size", "Prior_Size"
)
label_map <- setNames(label_values, old_names)


################################################################################
# 3. STENGARD REPLICATION: SUBJECT-LEVEL ABC & VALIDATION
################################################################################

# ---- 3.1 Run subject-level ABC + PPC for Stengard --------------------------
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


# Unpack for convenience
abc_indiv_s       <- stengard_subject_results$abc_indiv
observed_indiv_s  <- stengard_subject_results$observed_indiv
post_pred_indiv_s <- stengard_subject_results$post_pred_indiv
posterior_indiv_s <- stengard_subject_results$posterior_indiv
stats_long_indiv_s <- stengard_subject_results$stats_long_indiv
ppp_by_subject_s   <- stengard_subject_results$ppp_by_subject
ppp_summary_s      <- stengard_subject_results$ppp_summary_overall

print(ppp_summary_s)
#n_total median_ppp q25_ppp q75_ppp prop_extreme
#<int>   <dbl>      <dbl>   <dbl>      <dbl>
#7826       0.57   0.226   0.862       0.181


# ---- 3.2 PPP distribution plot (Stengard, individual-level) ----------------
alpha <- 0.05

ppp_hist_s <- ggplot(ppp_by_subject_s, aes(x = ppp)) +
  geom_histogram(
    bins     = 30,
    color    = "white",
    fill = "#3C5488FF",
    alpha    = 0.9,
    position = "identity"
  ) +
  geom_vline(
    xintercept = c(alpha/2, 1 - alpha/2),
    linetype   = "dashed"
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "PPP distribution (Stengard, individual level)",
    x     = "Posterior predictive p-value (PPP)",
    y     = "Count",
    fill  = ""
  ) +
  theme_classic(base_size = 10)

ggsave(file.path(plot_dir_tiff, "PPP_hist_Stengard_indiv.tiff"), ppp_hist_s,
       width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPP_hist_Stengard_indiv.png"),  ppp_hist_s,
       width = 6, height = 4, dpi = 300)


# ---- 3.3 human vs posterior prediction (Stengard) -------------------
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

final_plot_s

ggsave(file.path(plot_dir_tiff, "PPC_Stengard.tiff"),
       final_plot_s, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPC_Stengard.png"),
       final_plot_s, width = 6, height = 4, dpi = 300)


################################################################################
# 4. EXPERIMENTAL DATA: SUBJECT-LEVEL ABC & VALIDATION
################################################################################

# ---- 4.1 Run subject-level ABC + PPC for Experimental ----------------------
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


# Unpack for convenience
abc_indiv_exp       <- experimental_subject_results$abc_indiv
observed_indiv_exp  <- experimental_subject_results$observed_indiv
post_pred_indiv_exp <- experimental_subject_results$post_pred_indiv
posterior_indiv_exp <- experimental_subject_results$posterior_indiv
stats_long_indiv_exp <- experimental_subject_results$stats_long_indiv
ppp_by_subject_exp   <- experimental_subject_results$ppp_by_subject
ppp_summary_exp      <- experimental_subject_results$ppp_summary_overall

print(ppp_summary_exp)
# n_total median_ppp q25_ppp q75_ppp prop_extreme(two side)
# <int>      <dbl>   <dbl>   <dbl>        <dbl>
# 7080      0.596   0.228   0.926        0.268

# ---- 4.2 PPP distribution plot (Experimental, individual-level) ------------
ppp_hist_exp <-ggplot(ppp_by_subject_exp, aes(x = ppp)) +
  geom_histogram(
    bins     = 30,
    color    = "white",
    fill = "#3C5488FF",
    alpha    = 0.9,
    position = "identity"
  ) +
  geom_vline(
    xintercept = c(alpha/2, 1 - alpha/2),
    linetype   = "dashed"
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "PPP distribution (Experimental, individual level)",
    x     = "Posterior predictive p-value (PPP)",
    y     = "Count",
    fill  = ""
  ) +
  theme_classic(base_size = 10)

ggsave(file.path(plot_dir_tiff, "PPP_hist_Experimental_indiv.tiff"),
       ppp_hist_exp, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPP_hist_Experimental_indiv.png"),
       ppp_hist_exp, width = 6, height = 4, dpi = 300)


# ---- 4.3 Pooled human vs posterior prediction (Experimental) ---------------
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

ggsave(file.path(plot_dir_tiff, "PPC_Experimental.tiff"),
       final_plot, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPC_Experimental.png"),
       final_plot, width = 6, height = 4, dpi = 300)


################################################################################
# 5. Group Level Posterior 
################################################################################
# 1. Compute human summary by format
human_sum_fmt_exp <- compute_all_metrics(
  df         = human_dt,   # your S dataset
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


# 2. Run ABC by format
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

#      p8_1      p8_2      p8_3       q_1       q_2       q_3       q_4       q_5       q_6    N9        v9      format
#1 0.2584371 0.3745074 0.3670555 0.1649831 0.2405736 0.1427853 0.1433336 0.1559266 0.1523978 4.632 0.7904572 probability
#2 0.4287978 0.2647282 0.3064740 0.1919803 0.1648246 0.1372889 0.1994717 0.1581302 0.1483050 5.816 0.6251510   frequency

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

ggsave(file.path(plot_dir_tiff, "Prior_vs_Posterior_Experimental_population.tiff"),
       final_plot_group_exp, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Prior_vs_Posterior_Experimental_population.png"),
       final_plot_group_exp, width = 8, height = 8, dpi = 300)

# 1. Compute human summary by format
human_sum_fmt_s <- compute_all_metrics(
  df         = human_dt_s,   # your S dataset
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


# 2. Run ABC by format
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
#     p8_1      p8_2      p8_3       q_1       q_2       q_3       q_4       q_5       q_6    N9        v9      format
#1 0.2734783 0.3525548 0.3739673 0.1601951 0.1653180 0.1690286 0.1663143 0.1557877 0.1833566 5.078 0.6794729   frequency
#2 0.2396405 0.4033169 0.3570426 0.1677990 0.1648174 0.1469340 0.1613092 0.1671601 0.1919803 5.236 0.7450428 probability
plot_list_group_s <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_param (
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

ggsave(file.path(plot_dir_tiff, "Prior_vs_Posterior_Stengard_population.tiff"),
       final_plot_group_s, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Prior_vs_Posterior_Stengard_population.png"),
       final_plot_group_s, width = 8, height = 8, dpi = 300)


################################################################################
# END OF SCRIPT
################################################################################

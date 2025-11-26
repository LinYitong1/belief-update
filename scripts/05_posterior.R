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

human_data_format <- human_data_format %>%
  distinct(subject, format)%>%
  dplyr::mutate(
    subject = as.double(subject), 
    format  = as.character(format)  
  )


human_dt_s <- as.data.table(load_clean_data("data/s.csv",  type = "stengard"))


human_data_s_format <- human_data_s_format %>%
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
#7826    0.57       0.228    0.86       0.0832


# ---- 3.2 PPP distribution plot (Stengard, individual-level) ----------------
ppp_hist_s <- ggplot(ppp_by_subject_s, aes(x = ppp, fill = bad_fit)) +
  geom_histogram(
    bins     = 30,
    color    = "white",
    alpha    = 0.9,
    position = "identity"
  ) +
  scale_fill_manual(
    values = c(`TRUE` = "#F39B7FFF", `FALSE` = "#3C5488FF"),
    labels = c(`TRUE` = "PPP < .05 (bad fit)", `FALSE` = "PPP ≥ .05")
  ) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
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

# ---- 3.3 Population-level prior vs posterior (subject means, Stengard) -----
# Combine all posterior draws per subject and compute posterior mean per subject
post_draws_s <- data.table::rbindlist(
  abc_indiv_s$abc_out$post_draws,
  use.names = TRUE,
  fill      = TRUE
)

post_subject_means_s <- post_draws_s %>%
  dplyr::group_by(subject) %>%
  dplyr::summarise(
    dplyr::across(all_of(all_params), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Plot prior vs posterior (distribution of subject-level posterior means)
plot_list_pop_s <- lapply(all_params, function(param) {
  plot_param_population(
    param_name         = param,
    sim_params         = parameter_dt_s,
    post_subject_means = post_subject_means_s
  ) +
    labs(x = label_map[[param]])
})

final_plot_pop_s <- wrap_plots(plot_list_pop_s, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Prior_vs_Posterior_Stengard_population.tiff"),
       final_plot_pop_s, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Prior_vs_Posterior_Stengard_population.png"),
       final_plot_pop_s, width = 8, height = 8, dpi = 300)

# ---- 3.4 Pooled human vs posterior prediction (Stengard) -------------------
plot_overlay <- function(human_df, model_df, ref_lines, format_label) {
  combined_df <- bind_rows(
    mutate(human_df, type = "Human"),
    mutate(model_df, type = "Model")
  )
  
  ggplot() +
    geom_histogram(data = combined_df,
                   aes(x = response_pct, y = after_stat(density), fill = type, color = type),
                   binwidth = 3, alpha = 0.7, position = "identity") +
    geom_vline(data = ref_lines, aes(xintercept = value, colour = heuristic),
               linetype = "dashed", size = 0.5) +
    scale_fill_manual(values = c("Human" = "#3c5488", "Model" = "#f39b7f")) +
    scale_color_manual(values = c("Human" = "#3c5488", "Model" = "#f39b7f")) +
    facet_wrap(~ condition, ncol = 1, scales = "free_y") +
    add_okabe_color() +
    coord_cartesian(ylim = c(0,0.2)) +
    labs(title = format_label, x = "Estimates (%)", y = "Density") +
    theme_bw(base_size = 10) +                
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom",
      legend.key.height  = unit(3, "mm"),
      legend.text        = element_text(size = 7), 
      axis.title.x       = element_text(size = 8), 
      axis.title.y       = element_text(size = 8),  
      axis.text.x        = element_text(size = 7), 
      axis.text.y        = element_text(size = 7),  
      strip.text         = element_text(size = 8, face = "bold")  
    )
}
pooled_ppc_plot_s <- plot_human_vs_posterior_pooled(
  human_data      = human_dt_s,
  post_pred_indiv = post_pred_indiv_s,
  binwidth        = 0.05
)

ggsave(file.path(plot_dir_tiff, "Pooled_PPC_Stengard.tiff"),
       pooled_ppc_plot_s, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "Pooled_PPC_Stengard.png"),
       pooled_ppc_plot_s, width = 6, height = 4, dpi = 300)


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

# n_total median_ppp q25_ppp q75_ppp prop_extreme(<0.05)
# <int>      <dbl>   <dbl>   <dbl>        <dbl>
# 7080      0.596   0.226   0.928        0.111

# ---- 4.2 PPP distribution plot (Experimental, individual-level) ------------
ppp_hist_exp <- ggplot(ppp_by_subject_exp, aes(x = ppp, fill = bad_fit)) +
  geom_histogram(
    bins     = 30,
    color    = "white",
    alpha    = 0.9,
    position = "identity"
  ) +
  scale_fill_manual(
    values = c(`TRUE` = "#F39B7FFF", `FALSE` = "#3C5488FF"),
    labels = c(`TRUE` = "PPP < .05 (bad fit)", `FALSE` = "PPP ≥ .05")
  ) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
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

# ---- 4.3 Population-level prior vs posterior (subject means, Experimental) -
post_draws_exp <- data.table::rbindlist(
  abc_indiv_exp$abc_out$post_draws,
  use.names = TRUE,
  fill      = TRUE
)

post_subject_means_exp <- post_draws_exp %>%
  dplyr::group_by(subject) %>%
  dplyr::summarise(
    dplyr::across(all_of(all_params), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

plot_list_pop_exp <- lapply(all_params, function(param) {
  plot_param_population(
    param_name         = param,
    sim_params         = parameter_dt,
    post_subject_means = post_subject_means_exp
  ) +
    labs(x = label_map[[param]])
})

final_plot_pop_exp <- wrap_plots(plot_list_pop_exp, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Prior_vs_Posterior_Experimental_population.tiff"),
       final_plot_pop_exp, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Prior_vs_Posterior_Experimental_population.png"),
       final_plot_pop_exp, width = 8, height = 8, dpi = 300)

# ---- 4.4 Pooled human vs posterior prediction (Experimental) ---------------
pooled_ppc_plot_exp <- plot_human_vs_posterior_pooled(
  human_data      = human_dt,
  post_pred_indiv = post_pred_indiv_exp,
  binwidth        = 0.05
)

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
# 5. OPTIONAL: SUMMARY OF SUBJECT-LEVEL POSTERIOR MEANS (POPULATION VIEW)
################################################################################
# This keeps everything individual-level: we first compute posterior mean
# per subject, then optionally average across subjects per dataset.
human_data_s_format_clean <- human_data_s_format %>%
  # Convert the subject column to character type
  dplyr::mutate(
    subject = as.character(subject)
  )
posterior_means_pop_s <- post_subject_means_s %>%
  left_join(human_data_s_format_clean, by = "subject") %>%  
  group_by(format) %>%  
  dplyr::summarise(
    dplyr::across(all_of(all_params), median),
    .groups = "drop"
  ) %>%
  dplyr::mutate(dataset = "Stengard")

# 1. Create a clean version of the format data with subject as character
human_data_format_clean <- human_data_format %>%
  # Convert the subject column to character type
  dplyr::mutate(
    subject = as.character(subject)
  )

# 2. Perform the join and summary using the cleaned data
posterior_means_pop_exp <- post_subject_means_exp %>%
  # Use the cleaned format data for the join
  left_join(human_data_format_clean, by = "subject") %>%  
  group_by(format) %>%  
  dplyr::summarise(
    dplyr::across(all_of(all_params), median),
    .groups = "drop"
  ) %>%
  dplyr::mutate(dataset = "Experimental")


posterior_means_pop <- dplyr::bind_rows(
  posterior_means_pop_s,
  posterior_means_pop_exp
) %>%
  dplyr::rename_with(~ label_map[.], dplyr::all_of(all_params))

print(posterior_means_pop)
#format      `P(Integrated)` `P(Heuristic)` `P(Random Noise)` `P(BO)` `P(HO)` `P(FO)` `P(JO)` `P(LS)` `P(50%)` Sample_Size Prior_Size
#<chr>              <dbl>          <dbl>             <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>       <dbl>      <dbl>
#1 frequency        0.327          0.315             0.272   0.174   0.172   0.162   0.155   0.151    0.168        4.97      0.651
#2 probability      0.281          0.367             0.287   0.155   0.161   0.152   0.165   0.164    0.167        4.94      0.648
#3 frequency        0.348          0.289             0.273   0.169   0.169   0.149   0.173   0.156    0.160        5.00      0.643
#4 probability      0.244          0.372             0.319   0.140   0.178   0.152   0.141   0.157    0.155        4.86      0.648
################################################################################
# END OF SCRIPT
################################################################################

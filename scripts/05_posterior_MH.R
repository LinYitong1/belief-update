
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
human_dt_s <- as.data.table(load_clean_data("data/s.csv",  type = "stengard"))
human_data_s_format <- human_dt_s %>%
  distinct(subject, format)%>%
  dplyr::mutate(
    subject = as.character(subject), 
    format  = as.character(format)  
  )

# ---- 2.2 Simulation data & priors: Stengard --------------------------------
simulations_s_MH  <- as.data.table(readRDS("data/Simulate_Summary_dts.rds")) %>%
  dplyr::filter(model == "MH")
parameter_dt_s <- readRDS("data/parameter_dt_s.rds")

# ---- 2.3 Model comparison results (optional, not used directly here) -------
posterior_long_s <- readRDS("data/posterior_long_s.rds")

# ---- 2.4 Parameter specification -------------------------------------------
all_params <- c(paste0("p7_", 1:7))

old_names <- all_params

label_values <- c("P(BO)", "P(HO)", "P(FO)",
                  "P(JO)", "P(LS)", "P(50%)", "P(Random)")
label_map <- setNames(label_values, old_names)
################################################################################
# 3. SUBJECT-LEVEL Posterior Check
################################################################################
stengard_subject_MH_results <- run_subject_level_MH_validation(
  dataset_name      = "Stengard (MH, individual)",
  human_data        = human_dt_s,
  sim_params        = parameter_dt_s,
  simulations       = simulations_s_MH,
  subject_col       = "subject",  
  abc_tol           = 0.05,
  ppp_alpha         = 0.05,
  calculate_variance = FALSE
)

stengard_subject_MH_results$ppp_summary_overall

# Unpack for convenience
abc_indiv_s_MH       <- stengard_subject_MH_results$abc_indiv
observed_indiv_s_MH  <- stengard_subject_MH_results$observed_indiv
post_pred_indiv_s_MH <- stengard_subject_MH_results$post_pred_indiv
posterior_indiv_s_MH <- stengard_subject_MH_results$posterior_indiv
stats_long_indiv_s_MH <- stengard_subject_MH_results$stats_long_indiv
ppp_by_subject_s_MH   <- stengard_subject_MH_results$ppp_by_subject
ppp_summary_s_MH      <- stengard_subject_MH_results$ppp_summary_overall

print(ppp_summary_s_MH)
#n_total median_ppp q25_ppp q75_ppp prop_extreme
#<int>   <dbl>      <dbl>   <dbl>      <dbl>
#7826      0.548   0.146   0.858        0.174


# ---- 3.2 PPP distribution plot (Stengard, individual-level) ----------------
ppp_hist_s_MH <- ggplot(ppp_by_subject_s_MH, aes(x = ppp, fill = bad_fit)) +
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

ggsave(file.path(plot_dir_tiff, "PPP_hist_Stengard_indiv_MH.tiff"), ppp_hist_s_MH,
       width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPP_hist_Stengard_indiv_MH.png"),  ppp_hist_s_MH,
       width = 6, height = 4, dpi = 300)

# ---- 3.3 Population-level prior vs posterior (subject means, Stengard) -----
# Combine all posterior draws per subject and compute posterior mean per subject
post_draws_s_MH <- data.table::rbindlist(
  abc_indiv_s_MH$abc_out$post_draws,
  use.names = TRUE,
  fill      = TRUE
)

post_subject_means_s_MH <- post_draws_s_MH %>%
  dplyr::group_by(subject) %>%
  dplyr::summarise(
    dplyr::across(all_of(all_params), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Plot prior vs posterior (distribution of subject-level posterior means)
plot_list_pop_s_MH <- lapply(all_params, function(param) {
  plot_param_population(
    param_name         = param,
    sim_params         = parameter_dt_s,
    post_subject_means = post_subject_means_s_MH
  ) +
    labs(x = label_map[[param]])
})


final_plot_pop_s_MH <- wrap_plots(plot_list_pop_s_MH, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "Prior_vs_Posterior_Stengard_population_MH.tiff"),
       final_plot_pop_s_MH, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Prior_vs_Posterior_Stengard_population_MH.png"),
       final_plot_pop_s_MH, width = 8, height = 8, dpi = 300)

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
pooled_ppc_plot_s_MH <- plot_human_vs_posterior_pooled(
  human_data      = human_dt_s,
  post_pred_indiv = post_pred_indiv_s_MH,
  binwidth        = 0.05
)

ggsave(file.path(plot_dir_tiff, "Pooled_PPC_Stengard.tiff"),
       pooled_ppc_plot_s_MH, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "Pooled_PPC_Stengard.png"),
       pooled_ppc_plot_s_MH, width = 6, height = 4, dpi = 300)

ref_points_s <- tibble(
  BR  = c(0.1, 0.7, 0.9), 
  HR  = c(0.5, 0.9, 0.9), 
  FAR = c(0.3, 0.5, 0.1)
)

ref_lines_s_MH <- get_reference_lines(human_dt_s, ref_points_s)

human_df_s_prob <- human_dt_s %>% 
  filter(format == "probability") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  filter(condition %in% unique(ref_lines_s_MH$condition))   

human_df_s_freq <- human_dt_s %>% 
  filter(format == "frequency") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = response * 100
  ) %>%
  filter(condition %in% unique(ref_lines_s_MH$condition))


model_df_prob_s_MH <- post_pred_indiv_s_MH %>% 
  filter(format == "probability") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  filter(condition %in% unique(ref_lines_s_MH$condition))

model_df_freq_s_MH <- post_pred_indiv_s_MH %>% 
  filter(format == "frequency") %>%
  mutate(
    condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
    response_pct = prediction * 100
  ) %>%
  filter(condition %in% unique(ref_lines_s_MH$condition))

plot_prob_overlay_s_MH <- plot_overlay(
  human_df  = human_df_s_prob,
  model_df  = model_df_prob_s_MH,
  ref_lines = ref_lines_s_MH,
  format_label = "Probability Format"
)

plot_freq_overlay_s_MH <- plot_overlay(
  human_df  = human_df_s_freq,
  model_df  = model_df_freq_s_MH,
  ref_lines = ref_lines_s_MH,
  format_label = "Frequency Format"
)

final_plot_s_MH <- plot_prob_overlay_s_MH + plot_freq_overlay_s_MH +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir_tiff, "PPC_Stengard_MH.tiff"),
       final_plot_s_MH, width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPC_Stengard_MH.png"),
       final_plot_s_MH, width = 6, height = 4, dpi = 300)


################################################################################
# END OF SCRIPT
################################################################################

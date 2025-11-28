
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
#7826      0.546   0.146   0.856        0.276

alpha<-0.05
# ---- 3.2 PPP distribution plot (Stengard, individual-level) ----------------
ppp_hist_s_MH <-ggplot(ppp_by_subject_s_MH, aes(x = ppp)) +
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
  


ggsave(file.path(plot_dir_tiff, "PPP_hist_Stengard_indiv_MH.tiff"), ppp_hist_s_MH,
       width = 6, height = 4, dpi = 600)
ggsave(file.path(plot_dir_png,  "PPP_hist_Stengard_indiv_MH.png"),  ppp_hist_s_MH,
       width = 6, height = 4, dpi = 300)
# ---- 3.4 Pooled human vs posterior prediction (Stengard) -------------------

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
# 5. Group Level Posterior 
################################################################################
human_sum_fmt_s_MH  <- compute_all_metrics(
  df         = human_dt_s,   # your S dataset
  column     = "response",
  group_vars = c("format")
)

slope_int_fmt_s_MH  <- compute_SI_by(
  dt         = human_dt_s,
  group_vars = c("format"),
  predictors = c("BR", "HR", "FAR"),
  column     = "response"
)

human_summary_fmt_s_MH  <- human_sum_fmt_s_MH  %>%
  dplyr::left_join(slope_int_fmt_s_MH,   by = "format") 


# 2. Run ABC by format
abc_group_s_MH  <- compute_abc_by_format_MH(
  human_summary = human_summary_fmt_s_MH ,
  sim_params    = parameter_dt_s,
  simulations   = simulations_s_MH ,
  tol           = 0.05
)

posterior_means_s_MH  <- lapply(names(abc_group_s_MH $post_draws), function(fmt) {
  abc_group_s_MH $post_draws[[fmt]] %>%
    dplyr::summarise(dplyr::across(all_of(all_params), ~ mean(.x, na.rm = TRUE))) %>%
    dplyr::mutate(format = fmt)
}) %>%
  dplyr::bind_rows()
#      p7_1      p7_2      p7_3      p7_4      p7_5      p7_6      p7_7      format
#1  0.1518931 0.1496845 0.1326903 0.1218902 0.1127533 0.1640404 0.1670487   frequency
#20.1683713 0.1456248 0.1169706 0.1343552 0.1179140 0.1315043 0.1852600 probability

plot_list_group_s_MH  <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_group(
    param_name = param,
    sim_params = parameter_dt_s,
    abc_fmt    = abc_group_s_MH ,
    x_label    = label_map[[param]]
  )
})

final_plot_group_s_MH  <- patchwork::wrap_plots(plot_list_group_s_MH , ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")


ggsave(file.path(plot_dir_tiff, "Prior_vs_Posterior_Stengard_population_MH.tiff"),
       final_plot_group_s_MH, width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png,  "Prior_vs_Posterior_Stengard_population_MH.png"),
       final_plot_group_s_MH, width = 8, height = 8, dpi = 300)

################################################################################
# END OF SCRIPT
################################################################################

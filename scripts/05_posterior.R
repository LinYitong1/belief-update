################################################################################
# ABC POSTERIOR INFERENCE: HIERARCHICAL HEURISTIC-ANALYTICAL MODEL (11 PARAMS)
################################################################################
#
# DESCRIPTION:
# Bayesian parameter estimation for belief updating models using Approximate
# Bayesian Computation (ABC). Analyzes both Stengard replication data and
# novel experimental data across probability and frequency presentation formats.
#
# ANALYSIS WORKFLOW:
# 1. Environment Setup: Load packages, functions, configurations
# 2. Data Loading: Import empirical data, simulations, priors
# 3. Stengard Analysis: ABC inference and validation
# 4. Experimental Analysis: Parallel analysis on new data
# 5. Cross-Dataset Synthesis: Comparative parameter estimates
#
# METHODS:
# - ABC with neural network adjustment (tolerance = 0.05)
# - Posterior predictive checks via overlap coefficients
# - Parameter recovery assessment through prior-posterior comparison
#
# Author: Yitong Lin
################################################################################


################################################################################
# 1. ENVIRONMENT SETUP
################################################################################

# ---- 1.1 Package Management ------------------------------------------------
# Install and load required packages for data processing, ABC inference,
# and visualization

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

# ---- 1.2 Source Custom Functions -------------------------------------------
# Load modular helper functions for metrics, models, visualization, and
# posterior analysis

source("R/metrics.R")
source("R/models.R")
source("R/plot_helpers.R")
source("R/table_helper.R") 
source("R/posterior_helper.R")

# ---- 1.3 Output Configuration ----------------------------------------------
# Configure paths and parameters for figure export

plot_dir_tiff    <- "fig/tiff"
plot_dir_png    <- "fig/png"
default_dpi <- 300
plot_width  <- 6
plot_height <- 4.5


################################################################################
# 2. DATA LOADING
################################################################################

# ---- 2.1 Empirical Data ----------------------------------------------------
# Trial-level responses from participants in both probability and frequency
# presentation formats

human_dt <- as.data.table(load_clean_data("data/df.csv", type = "experiment"))
human_dt_s  <- as.data.table(load_clean_data("data/s.csv",   type = "stengard"))

# ---- 2.2 Simulation Data: Experimental Study -------------------------------
# Model predictions from prior parameter samples for the M_L (likelihood-based)
# model variant

simulations <- as.data.table(readRDS("data/Simulate_Summary_dt.rds"))%>%
  filter(model=="M_L")
observed    <-  as.data.table(readRDS("data/Human_Summary_dt.rds"))
parameter_dt <- readRDS("data/parameter_dt.rds")

# ---- 2.3 Simulation Data: Stengard Replication -----------------------------
# Corresponding simulations for replication dataset

simulations_s <- as.data.table(readRDS("data/Simulate_Summary_dts.rds"))%>%
  filter(model=="M_L")
observed_s    <-  as.data.table(readRDS("data/Human_Summary_dts.rds"))
parameter_dt_s <- readRDS("data/parameter_dt_s.rds") 

# ---- 2.4 Individual-Level Model Fits ---------------------------------------
# Posterior model probabilities from participant-level model comparison
# (used for subgroup identification)

posterior_long<-readRDS("data/posterior_long.rds") 
posterior_long_s<-readRDS("data/posterior_long_s.rds") 

# ---- 2.5 Parameter Specification -------------------------------------------
# Define 11-parameter model structure

all_params <- c(paste0("p8_", 1:3), paste0("q_", 1:6), "N9", "v9")


################################################################################
# 3. STENGARD REPLICATION: ABC INFERENCE
################################################################################

# ---- 3.1 Run ABC Posterior Inference ---------------------------------------
# Estimate posterior distributions using ABC rejection with neural network
# adjustment. Tolerance of 0.05 retains top 5% of simulations based on
# summary statistic distance.

stengard_results <- run_analysis_pipeline_overlap(
  dataset_name       = "Stengard",
  human_data         = human_dt_s,
  sim_params         = parameter_dt_s,
  simulations        = simulations_s,
  ppc_subject_id     = "1_1",
  subject_col        = "subject",
  ref_points         = tibble(BR = c(0.1, 0.7, 0.9), 
                              HR = c(0.5, 0.9, 0.9), 
                              FAR = c(0.3, 0.5, 0.1)),
  calculate_variance = FALSE,
  abc_tol =0.05
)

# ---- 3.2 Export Posterior Predictive Visualizations -----------------------
# Overlayed histograms showing empirical data vs posterior predictions

ggsave(file.path(plot_dir_tiff, "Figure_8.tiff"), stengard_results$plot,
       width = 8, height = 5, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_8.png"), stengard_results$plot,
       width = 8, height = 5, dpi = 300)

# ---- 3.3 Extract Posterior Predictions -------------------------------------

post_s<-as.data.table(stengard_results$posteriors_prediction)

# ---- 3.4 Compute Summary Statistics from Posteriors -----------------------
# Aggregate posterior predictions for model validation

posterior_s <- summarise_posterior(
  post_pred= post_s,
  calculate_variance = FALSE)  

stats_long_s <- prepare_stats_long(
  posterior_all =  posterior_s,
  observed_df   = observed_s,
  human_dt      = human_dt_s,
  id_var        = "subject_s")  

# ---- 3.5 Distribution Overlap Coefficient ----------------------------------
# Quantify agreement between observed and predicted distributions using
# histogram overlap (Inman & Bradley, 1989)

overlap_tbl_s<-get_overlap_tbl(stats_long_s, bins = 30)

# ---- 3.6 Posterior Predictive P-values -------------------------------------
# One-sided PPP values (Schmidt et al., 2023) to assess model adequacy
# Values near 0 or 1 indicate systematic misfit

ppp_results_s <- calculate_ppp(stats_long_s, alpha = 0.05)
min(ppp_results_s$ppp)
#  0.148

# ---- 3.7 Parameter Recovery: Prior vs Posterior ----------------------------
# Assess parameter identifiability through posterior contraction

sim_params_agg_s <- aggregate_sim_params(parameter_dt_s)

old_names <- c("p8_1", "p8_2", "p8_3",
               "q_1", "q_2", "q_3", "q_4", "q_5", "q_6",
               "N9", "v9")

# Compute posterior means by format
posterior_means_s_f <- stengard_results$posteriors_parameter$frequency %>%
  summarise(across(all_of(old_names), mean))%>%
  mutate(dataset = "S_F")

posterior_means_s_p <- stengard_results$posteriors_parameter$probability %>%
  summarise(across(all_of(old_names), mean))%>%
  mutate(dataset = "S_P")

posterior_means_s<-rbind(posterior_means_s_f,posterior_means_s_p)

# Define interpretable parameter labels
label_values <- c("P(Integrated)", "P(Heuristic)", "P(Random Noise)",
               "P(BO)", "P(HO)", "P(FO)",
               "P(JO)", "P(LS)", "P(50%)",
               "Sample_Size", "Prior_Size")
label_map <- setNames(label_values, old_names)

# ---- 3.8 Visualize Prior-Posterior Shift -----------------------------------
# Generate density plots comparing prior and posterior for each parameter

plot_list_s <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_param(
    param_name = param,
    prior_sampler =  function(n) sim_params_agg_s[[param]],
    data = stengard_results,
    x_label = label_map[[param]],
  )
})

final_plot_s <- wrap_plots(plot_list_s, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom") 

ggsave(file.path(plot_dir_tiff, "Apendix_E2.tiff"),final_plot_s ,
       width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png, "Apendix_E2.png"), final_plot_s ,
       width = 8, height = 8, dpi = 300)


################################################################################
# 4. EXPERIMENTAL DATA: ABC INFERENCE
################################################################################

# ---- 4.1 Run ABC Posterior Inference ---------------------------------------
# Apply identical methodology to experimental dataset
# Note: calculate_variance=TRUE for variance decomposition analysis

experimental_results <- run_analysis_pipeline_overlap(
  dataset_name       = "Experimental",
  human_data         = human_dt,
  sim_params         = parameter_dt,
  simulations        = simulations,
  ppc_subject_id     = "1",
  subject_col        = "subject",
  ref_points         = tibble(BR = c(0.01,  0.4,  0.97),
                              HR = c( 0.72 , 0.53 ,0.53 ), 
                              FAR = c(0.42,0.31, 0.11)),
  calculate_variance = TRUE,
  abc_tol =0.05)

ggsave(file.path(plot_dir_tiff, "Figure_14.tiff"), experimental_results$plot,
       width = 8, height = 5, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_14.png"), experimental_results$plot,
       width = 8, height = 5, dpi = 300)

# ---- 4.2 Extract Posterior Means by Format ---------------------------------

posterior_means_e_f <- experimental_results$posteriors_parameter$frequency %>%
  summarise(across(all_of(old_names), mean)) %>%
  mutate(dataset = "E_F")

posterior_means_e_p <- experimental_results$posteriors_parameter$probability %>%
  summarise(across(all_of(old_names), mean))%>%
  mutate(dataset = "E_P")

posterior_means_e<-rbind(posterior_means_e_f,posterior_means_e_p)

# ---- 4.3 Secondary Figure Output -------------------------------------------

ggsave(file.path(plot_dir, "Figure_24.png"), experimental_results$plot,
       width = plot_width, height = plot_height, dpi = default_dpi)

# ---- 4.4 Extract Posterior Predictions -------------------------------------

post_exp<-as.data.table(experimental_results$posteriors_prediction)

# ---- 4.5 Posterior Predictive Validation -----------------------------------

posterior_exp <- summarise_posterior(
  post_pred= post_exp,
  calculate_variance = TRUE)  

stats_long_exp <- prepare_stats_long(
  posterior_all =  posterior_exp,
  observed_df   = observed,
  human_dt      = human_dt,
  id_var        = "subject")  

# ---- 4.6 Distribution Overlap Analysis -------------------------------------

overlap_tbl_exp<-get_overlap_tbl( stats_long_exp, bins =30)

# ---- 4.7 Posterior Predictive P-values -------------------------------------

ppp_results_exp <- calculate_ppp(stats_long_exp, alpha = 0.05)
mean(ppp_results_exp$ppp)

# ---- 4.8 Parameter Recovery Analysis ---------------------------------------

sim_params_agg <- aggregate_sim_params(parameter_dt)

# ---- 4.9 Visualize Prior-Posterior Shift -----------------------------------

plot_list <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_param (
    param_name = param,
    prior_sampler = function(n) sim_params_agg[[param]],
    data = experimental_results,
   x_label = label_map[[param]],
  )
})

final_plot <- wrap_plots(plot_list, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom") 

final_plot
ggsave(file.path(plot_dir_tiff, "Apendix_E3.tiff"),final_plot ,
       width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png, "Apendix_E3.png"), final_plot ,
       width = 8, height = 8, dpi = 300)


################################################################################
# 5. CROSS-DATASET SYNTHESIS
################################################################################

# ---- 5.1 Aggregate Posterior Means -----------------------------------------
# Combine results from both datasets and formats

posterior_means<-rbind(posterior_means_s_f,posterior_means_s_p,
                       posterior_means_e_f,posterior_means_e_p)

posterior_means <- posterior_means %>%
  rename_with(~ label_map[.], all_of(old_names))
print(posterior_means)

# ---- 5.2 Generate Publication Tables ---------------------------------------

knitr::kable(
  posterior_means,
  caption = "Posterior means of HAmix model parameters for each experiment and format.",
  digits = 3
)

library(kableExtra)

latex_tbl <- knitr::kable(
  posterior_means,
  format  = "latex",
  booktabs = TRUE,
  caption = "Posterior means of HAmix model parameters for each experiment and format.",
  label   = "tab:hamix-means",
  digits  = 3,
  escape  = TRUE
) |>
  kable_styling(latex_options = c("hold_position", "scale_down"))

# ---- 5.3 Prepare Combined Posterior Distributions -------------------------
# Merge all posterior samples for cross-condition comparison

combined_posteriors_all <- bind_rows(
  "E_F" = experimental_results$posteriors_parameter$frequency %>% dplyr::select(all_of(old_names)),
  "E_P" = experimental_results$posteriors_parameter$probability %>% dplyr::select(all_of(old_names)),
  "S_F" = stengard_results$posteriors_parameter$frequency %>% dplyr::select(all_of(old_names)),
  "S_P" = stengard_results$posteriors_parameter$probability %>% dplyr::select(all_of(old_names)),
  .id = "dataset"
) %>%
  pivot_longer(
    cols = all_of(old_names),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    parameter = factor(parameter, levels = old_names, labels = label_map[old_names]),
    dataset = factor(dataset, levels = c("S_F", "S_P", "E_F", "E_P"),
                     labels = c("Stengard (Freq)", "Stengard (Prob)",
                                "Exp (Freq)", "Exp (Prob)"))
  )


################################################################################
# END OF SCRIPT
################################################################################

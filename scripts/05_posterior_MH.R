################################################################################
# ABC POSTERIOR INFERENCE: MATCHING HEURISTIC MODEL (7 PARAMETERS)
################################################################################
#
# DESCRIPTION:
# Bayesian parameter estimation for the Matching Heuristic (MH) model variant
# using Approximate Bayesian Computation. This simplified 7-parameter model
# excludes sample size and prior strength effects, focusing on strategy mixture
# weights for pure heuristic-based belief updating.
#
#
# ANALYSIS WORKFLOW:
# 1. Environment Setup: Load dependencies and configure output
# 2. Data Loading: Import MH simulations and empirical data
# 3. Stengard Analysis: ABC inference and posterior validation
# 4. Posterior Diagnostics: Overlap, PPP, parameter recovery
# 5. Visualization: Prior-posterior comparison
#
# METHODS:
# - ABC rejection sampling (tolerance = 0.05)
# - Posterior predictive checking
# - Distribution overlap coefficients
# - Format-specific parameter estimation
#
# Author: Yitong Lin
################################################################################


################################################################################
# 1. ENVIRONMENT SETUP
################################################################################

# ---- 1.1 Package Management ------------------------------------------------
# Load statistical, visualization, and data manipulation packages

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

# ---- 1.2 Source Helper Modules ---------------------------------------------
# Import modular functions for model evaluation and posterior analysis

source("R/metrics.R")
source("R/models.R")
source("R/plot_helpers.R")
source("R/posterior_helper.R")

# ---- 1.3 Configure Output Directories --------------------------------------
# Set paths for high-resolution figures

plot_dir_tiff    <- "fig/tiff"
plot_dir_png    <- "fig/png"


################################################################################
# 2. DATA LOADING
################################################################################

# ---- 2.1 Load Empirical Data -----------------------------------------------
# Trial-level participant responses from experimental and replication studies

human_dt <- as.data.table(load_clean_data("data/df.csv", type = "experiment"))
human_dt_s  <- as.data.table(load_clean_data("data/s.csv",   type = "stengard"))

# ---- 2.2 Load MH Model Simulations: Experimental ---------------------------
# Simulated summary statistics from MH model prior parameter samples

simulations_MH <- as.data.table(readRDS("data/Simulate_Summary_dt.rds"))%>%
  filter(model=="MH")
observed    <-  as.data.table(readRDS("data/Human_Summary_dt.rds"))
parameter_dt <- readRDS("data/parameter_dt.rds")

# ---- 2.3 Load MH Model Simulations: Stengard ------------------------------
# Corresponding simulations for replication dataset

simulations_s_MH <- as.data.table(readRDS("data/Simulate_Summary_dts.rds"))%>%
  filter(model=="MH")
observed_s    <-  as.data.table(readRDS("data/Human_Summary_dts.rds"))
parameter_dt_s <- readRDS("data/parameter_dt_s.rds") 

# ---- 2.4 Load Individual-Level Model Fits ----------------------------------
# Pre-computed posterior model probabilities (for reference)

posterior_long<-readRDS("data/posterior_long.rds") 
posterior_long_s<-readRDS("data/posterior_long_s.rds") 

# ---- 2.5 Define MH Model Parameters ----------------------------------------
# 7 mixture weights for different cognitive strategies

all_params <- c(paste0("p7_", 1:7))

################################################################################
# 3. STENGARD REPLICATION: ABC INFERENCE
################################################################################

# ---- 3.1 Run ABC Posterior Estimation --------------------------------------
# Execute complete inference pipeline:
#   - Summary statistics calculation
#   - ABC rejection sampling (5% tolerance)
#   - Posterior predictive simulations
#   - Overlayed histogram visualization

stengard_results_MH <- run_analysis_pipeline_overlap_MH(
  dataset_name       = "Stengard",
  human_data         = human_dt_s,
  sim_params         = parameter_dt_s,
  simulations        = simulations_s_MH,
  ppc_subject_id     = "1_1",
  subject_col        = "subject_s",
  ref_points         = tibble(BR = c(0.1, 0.7, 0.9), 
                              HR = c(0.5, 0.9, 0.9), 
                              FAR = c(0.3, 0.5, 0.1)),
  calculate_variance = FALSE,
  abc_tol =0.05
)

# ---- 3.2 Export Main Results Figure ----------------------------------------
# Save overlayed histograms showing data vs posterior predictions

ggsave(file.path(plot_dir_tiff, "Figure_7.tiff"), stengard_results_MH$plot,
       width = 8, height = 5, dpi = 600)
ggsave(file.path(plot_dir_png, "Figure_7.png"), stengard_results_MH$plot,
       width = 8, height = 5, dpi = 300)


################################################################################
# 4. POSTERIOR VALIDATION
################################################################################

# ---- 4.1 Extract Posterior Predictions -------------------------------------

post_s_MH<-as.data.table(stengard_results_MH$posteriors_prediction)

# ---- 4.2 Compute Summary Statistics ----------------------------------------
# Aggregate posterior samples across trials and formats

posterior_s_MH <- summarise_posterior(
  post_pred= post_s_MH,
  calculate_variance = FALSE)  

# ---- 4.3 Prepare Long-Format Data for Comparison ---------------------------
# Combine observed and model-predicted statistics

stats_long_s_MH <- prepare_stats_long(
  posterior_all =  posterior_s_MH,
  observed_df   = observed_s,
  human_dt      = human_dt_s,
  id_var        = "subject_s")  

# ---- 4.4 Calculate Distribution Overlap ------------------------------------
# Histogram-based overlap coefficient (30 bins)

overlap_tbl_s_MH<-get_overlap_tbl(stats_long_s_MH, bins = 30)

# ---- 4.5 Compute Posterior Predictive P-values ----------------------------
# Bayesian goodness-of-fit test (α = 0.05)

ppp_results_s_MH <- calculate_ppp(stats_long_s_MH, alpha = 0.05)
min(ppp_results_s_MH$ppp)
# 0.098


################################################################################
# 5. PARAMETER POSTERIOR ANALYSIS
################################################################################

# ---- 5.1 Aggregate Prior Samples -------------------------------------------

sim_params_agg_s_MH <- aggregate_sim_params_MH(parameter_dt_s)
old_names <- c("p7_1", "p7_2", "p7_3","p7_4","p7_5","p7_6","p7_7")

# ---- 5.2 Compute Posterior Means by Format ---------------------------------

# Frequency format
posterior_means_s_MH_f <- stengard_results_MH$posteriors_parameter$frequency %>%
  summarise(across(all_of(old_names), mean))%>%
  mutate(dataset = "S_F")

# Probability format
posterior_means_s_MH_p <- stengard_results_MH$posteriors_parameter$probability %>%
  summarise(across(all_of(old_names), mean))%>%
  mutate(dataset = "S_P")

posterior_means_s_MH<-rbind(posterior_means_s_MH_f,posterior_means_s_MH_p)

# ---- 5.3 Define Interpretable Labels ---------------------------------------
# Map parameter names to cognitive strategy descriptions

label_values <- c("P(BO)", "P(HO)", "P(FO)",
                  "P(JO)", "P(LS)", "P(50%)", "P(Random)")
label_map <- setNames(label_values, old_names)

# ---- 5.4 Generate Prior-Posterior Comparison Plots -------------------------
# Visualize parameter identifiability through posterior contraction

plot_list_s_MH <- lapply(all_params, function(param) {
  plot_prior_vs_posterior_param(
    param_name = param,
    prior_sampler =  function(n) sim_params_agg_s_MH [[param]],
    data = stengard_results_MH,
    x_label = label_map[[param]],
  )
})

# ---- 5.5 Combine into Multi-Panel Figure -----------------------------------

final_plot_s_MH<- wrap_plots(plot_list_s_MH, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom") 

# ---- 5.6 Export Parameter Recovery Figure ----------------------------------

ggsave(file.path(plot_dir_tiff, "Apendix_E1.tiff"),final_plot_s_MH,
       width = 8, height = 8, dpi = 600)
ggsave(file.path(plot_dir_png, "Apendix_E1.png"), final_plot_s_MH,
       width = 8, height = 8, dpi = 300)


################################################################################
# END OF SCRIPT
################################################################################

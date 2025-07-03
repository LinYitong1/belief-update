# Belief Updating Project: Refactored Posterior Analysis Pipeline
# ================================================================
# This script wraps repeated code into functions, clarifies the data flow
# and follows a clear top‑down pipeline organisation.  Each STEP from the
# original notebook now maps onto a named function that can be unit‑tested
# or reused across probability / frequency formats and different datasets.
# ----------------------------------------------------------------------

# ---- 1. LIBRARIES -----------------------------------------------------
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
source("R/metrics.R")
source("R/models.R")
source("R/plot_helpers.R")
source("R/posterior_helper.R")
plot_dir    <- "fig"
default_dpi <- 300
plot_width  <- 6
plot_height <- 4.5
# ---- 2. DATA ------------------------------------------------------
human_dt <- as.data.table(load_clean_data("data/df.xlsx", type = "experiment"))
human_dt_s  <- as.data.table(load_clean_data("data/s.csv",   type = "stengard"))

simulations <- as.data.table(readRDS("data/Simulate_Summary_dt.rds"))%>%
  filter(model=="M_L")
observed    <-  as.data.table(readRDS("data/Human_Summary_dt.rds"))
parameter_dt <- readRDS("data/parameter_dt.rds")

simulations_s <- as.data.table(readRDS("data/Simulate_Summary_dts.rds"))%>%
  filter(model=="M_L")
observed_s    <-  as.data.table(readRDS("data/Human_Summary_dts.rds"))
parameter_dt_s <- readRDS("data/parameter_dt_s.rds") 

###############################################################################
#Stengard                                                            #
###############################################################################
stengard_results <- run_analysis_pipeline(
  dataset_name       = "Stengard",
  human_data         = human_dt_s,
  sim_params         = parameter_dt_s,
  simulations        = simulations_s,
  ppc_subject_id     = "1_1",
  subject_col        = "subject_s",
  ref_points         = tibble(BR = c(0.1, 0.7, 0.9), 
                              HR = c(0.5, 0.9, 0.9), 
                              FAR = c(0.3, 0.5, 0.1)),
  calculate_variance = FALSE,
  abc_tol =0.01
)
ggsave(file.path(plot_dir, "Figure_12.png"), stengard_results$plot,
       width = plot_width, height = plot_height, dpi = default_dpi)
###############################################################################
#Experiment                                                           #
###############################################################################
experimental_results <- run_analysis_pipeline(
  dataset_name       = "Experimental",
  human_data         = human_dt,
  sim_params         = parameter_dt,
  simulations        = simulations,
  ppc_subject_id     = "1",
  subject_col        = "subject",
  ref_points         = tibble(BR = c(0.01,  0.4,  0.97),
                              HR = c( 0.72 , 0.9 ,0.53 ), 
                              FAR = c(0.42,0.42, 0.11)),
  calculate_variance = TRUE,
  abc_tol =0.01)
  ggsave(file.path(plot_dir, "Figure_13.png"), experimental_results $plot,
       width = plot_width, height = plot_height, dpi = default_dpi)


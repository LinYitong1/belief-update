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

post_s<-as.data.table(stengard_results$posteriors_prediction)

posterior_s <- summarise_posterior(
  post_pred= post_s,
  calculate_variance = FALSE)  

stats_long_s <- prepare_stats_long(
  posterior_all =  posterior_s,
  observed_df   = observed_s,
  human_dt      = human_dt_s,
  id_var        = "subject_s")  

density_stat_s <- plot_density_stats(
  stats_long_s,
  dens_ncol   = 4)
ggsave(file.path(plot_dir, "Figure_14.png"), density_stat_s,
       width = plot_width, height = plot_height, dpi = default_dpi)

overlap_tbl_s<-get_overlap_tbl(stats_long_s, bins = 30)

hist_overlap_s<-plot_hist_overlap(overlap_tbl_s)
ggsave(file.path(plot_dir, "Figure_16.png"), hist_overlap_s,
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

post_exp<-as.data.table(experimental_results$posteriors_prediction)

 posterior_exp <- summarise_posterior(
   post_pred= post_exp,
   calculate_variance = TRUE)  
 
 stats_long_exp <- prepare_stats_long(
   posterior_all =  posterior_exp,
   observed_df   = observed,
   human_dt      = human_dt,
   id_var        = "subject")  
 
 density_stat_exp<-plot_density_stats(
   stats_long_exp ,
   dens_ncol   = 4)
 ggsave(file.path(plot_dir, "Figure_15.png"), density_stat_exp,
        width = plot_width, height = plot_height, dpi = default_dpi)
 
 overlap_tbl_exp<-get_overlap_tbl( stats_long_exp, bins =30)
 
 hist_overlap_exp<-plot_hist_overlap(overlap_tbl_exp)
 ggsave(file.path(plot_dir, "Figure_17.png"),  hist_overlap_exp,
        width = plot_width, height = plot_height, dpi = default_dpi)
 
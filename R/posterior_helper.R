################################################################################
# POSTERIOR HELPER FUNCTIONS FOR BELIEF UPDATING MODELS
################################################################################
#
# This module implements Approximate Bayesian Computation (ABC) inference and
# posterior predictive checking for cognitive models of belief updating.
#
# Core functionality:
#   - ABC rejection sampling with neural network adjustment
#   - Posterior predictive simulations
#   - Model validation metrics (overlap coefficients, PPP values)
#   - Visualization utilities for posterior diagnostics
#
# Author: Yitong Lin
################################################################################


################################################################################
# 0. DEPENDENCIES
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(abc)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})


################################################################################
# 1. ABC POSTERIOR INFERENCE
################################################################################

# ---- 1.1 ABC Rejection Sampling for 11-Parameter Model --------------------
#' Estimate posterior distributions via ABC rejection sampling
#'
#' Implements neural network-adjusted ABC (Blum & François, 2010) with logit
#' transformation for bounded parameters and log transformation for unbounded.
#'
#' @param obs_stats   Observed summary statistics (numeric vector)
#' @param sim_params  Simulated parameter draws (data.frame with p8_1:3, q_1:6, N9, v9)
#' @param simulations Simulated summary statistics matching obs_stats structure
#' @param tol         Tolerance threshold for rejection (default: 0.05)
#' @return            abc object containing posterior samples

abc_posterior <- function(obs_stats,
                          sim_params,
                          simulations,
                          tol = 0.05) {
  PARAM_COLS <- c(paste0("p8_", 1:3), paste0("q_", 1:6), "N9", "v9")
  SUM_COLS   <- setdiff(names(simulations), c("Iteration", "model"))
  
  # Aggregate to iteration level
  sim_params <- as.data.table(sim_params)
  sim_params_agg <- unique(sim_params[, c("Iteration", PARAM_COLS), with = FALSE])
  setDT(sim_params_agg) 
  
  # Constrain probabilities to (ε, 1-ε) to avoid boundary issues
  prob_cols <- PARAM_COLS[1:9]              
  fix01     <- function(x, eps = 1e-4) pmin(pmax(x, eps), 1 - eps)
  sim_params_agg[, (prob_cols) := lapply(.SD, fix01), .SDcols = prob_cols]
  
  # Remove incomplete cases
  simulations   <- as.data.table(simulations)
  keep_rows     <- complete.cases(simulations[, ..SUM_COLS])
  sim_params_agg <- sim_params_agg[keep_rows]
  simulations    <- simulations[keep_rows]
  stopifnot(nrow(sim_params_agg) == nrow(simulations))
  
  # Transformation specifications
  trans_vec <- c(rep("logit", 9), "log", "log")
  bounds    <- cbind(rep(1e-4, 9), rep(1 - 1e-4, 9))
  
  # ABC with neural network adjustment
  abc::abc(
    target        = obs_stats,
    param         = sim_params_agg[, ..PARAM_COLS],
    sumstat       = simulations[, ..SUM_COLS],
    method        = "neuralnet",
    tol           = tol,
    transf        = trans_vec,
    logit.bounds  = bounds              
  )
}

# ---- 1.2 ABC Rejection Sampling for 7-Parameter MH Model ------------------
#' ABC inference for Matching Heuristic (MH) model variant
#'
#' Specialized version for 7-parameter model without sample size and prior strength.
#'
#' @inheritParams abc_posterior
#' @return abc object with 7 posterior parameter samples

abc_posterior_MH <- function(obs_stats,
                          sim_params,
                          simulations,
                          tol = 0.05) {
  PARAM_COLS <- c(paste0("p7_", 1:7))
  SUM_COLS   <- setdiff(names(simulations), c("Iteration", "model"))
  
  sim_params <- as.data.table(sim_params)
  sim_params_agg <- unique(sim_params[, c("Iteration", PARAM_COLS), with = FALSE])
  setDT(sim_params_agg) 
  
  prob_cols <- PARAM_COLS[1:7]              
  fix01     <- function(x, eps = 1e-4) pmin(pmax(x, eps), 1 - eps)
  sim_params_agg[, (prob_cols) := lapply(.SD, fix01), .SDcols = prob_cols]
  
  simulations   <- as.data.table(simulations)
  keep_rows     <- complete.cases(simulations[, ..SUM_COLS])
  sim_params_agg <- sim_params_agg[keep_rows]
  simulations    <- simulations[keep_rows]
  stopifnot(nrow(sim_params_agg) == nrow(simulations))
  
  trans_vec <- c(rep("logit", 7))
  bounds    <- cbind(rep(1e-4, 7), rep(1 - 1e-4, 7))  
  
  abc::abc(
    target        = obs_stats,
    param         = sim_params_agg[, ..PARAM_COLS],
    sumstat       = simulations[, ..SUM_COLS],
    method        = "neuralnet",
    tol           = tol,
    transf        = trans_vec,
    logit.bounds  = bounds              
  )
}


################################################################################
# 2. FORMAT-SPECIFIC ABC INFERENCE
################################################################################

# ---- 2.1 ABC by Presentation Format (11-Parameter Model) ------------------
#' Compute ABC posteriors separately for probability and frequency formats
#'
#' @param human_summary Summary statistics from empirical data with format column
#' @param sim_params    Parameter samples from prior
#' @param simulations   Simulated summary statistics
#' @param tol           ABC tolerance threshold
#' @return List containing abc results and posterior draws for each format

compute_abc_by_format <- function(human_summary,
                                  sim_params,
                                  simulations,
                                  tol = 0.05) {
  human_summary <- as_tibble(human_summary)
  sim_params    <- as_tibble(sim_params)
  simulations   <- as_tibble(simulations)
  
  stopifnot("format" %in% names(human_summary))
  formats <- unique(human_summary$format)
  
  SUM_COLS <- setdiff(names(simulations), c("Iteration","model"))
  
  abc_res    <- setNames(vector("list", length(formats)), formats)
  post_draws <- setNames(vector("list", length(formats)), formats)
  
  for (fmt in formats) {
    obs_stats <- human_summary %>%
      filter(format == fmt) %>%
      dplyr::select(all_of(SUM_COLS)) %>%
      unlist(use.names = FALSE)
    
    abc_obj <- abc_posterior(obs_stats, sim_params, simulations, tol)
    draws   <- as.data.frame(abc_obj$unadj.values)
    draws$subject <- paste0("global_", fmt)
    
    abc_res[[fmt]]    <- abc_obj
    post_draws[[fmt]] <- draws
  }
  
  list(abc_res = abc_res, post_draws = post_draws, tol = tol)
}

# ---- 2.2 ABC by Presentation Format (7-Parameter MH Model) ----------------
#' Format-specific ABC for MH model variant
#'
#' @inheritParams compute_abc_by_format
#' @return List of ABC results for MH model by format

compute_abc_by_format_MH <- function(human_summary,
                                  sim_params,
                                  simulations,
                                  tol = 0.05) {
  human_summary <- as_tibble(human_summary)
  sim_params    <- as_tibble(sim_params)
  simulations   <- as_tibble(simulations)
  
  stopifnot("format" %in% names(human_summary))
  formats <- unique(human_summary$format)
  
  SUM_COLS <- setdiff(names(simulations), c("Iteration","model"))
  
  abc_res    <- setNames(vector("list", length(formats)), formats)
  post_draws <- setNames(vector("list", length(formats)), formats)
  
  for (fmt in formats) {
    obs_stats <- human_summary %>%
      filter(format == fmt) %>%
      dplyr::select(all_of(SUM_COLS)) %>%
      unlist(use.names = FALSE)
    
    abc_obj <- abc_posterior_MH(obs_stats, sim_params, simulations, tol)
    draws   <- as.data.frame(abc_obj$unadj.values)
    draws$subject <- paste0("global_", fmt)
    
    abc_res[[fmt]]    <- abc_obj
    post_draws[[fmt]] <- draws
  }
  
  list(abc_res = abc_res, post_draws = post_draws, tol = tol)
}


################################################################################
# 3. POSTERIOR PREDICTIVE SIMULATIONS
################################################################################

# ---- 3.1 Generate Predictions (11-Parameter Model) ------------------------
#' Generate posterior predictive samples for trial-level predictions
#'
#' For each posterior parameter draw, simulate model predictions across all
#' trials including stochastic sample size effects.
#'
#' @param post_params Posterior parameter samples
#' @param df_raw      Trial-level data with BR, HR, FAR, true_posterior, response
#' @return Tibble with predictions for each posterior sample × trial combination

generate_ppc_global <- function(post_params, df_raw) {
  n_samples <- nrow(post_params)
  n_trials  <- nrow(df_raw)
  
  trial_fixed <- df_raw %>% dplyr::select(BR, HR, FAR, true_posterior, response)
  rep_param <- function(x) rep(as.numeric(x), n_trials)
  
  results <- map_dfr(seq_len(n_samples), function(i) {
    pars <- post_params[i, ]
    rf_vec <- simulate_and_mutate(df_raw, pars$N9)$relative_frequency
    
    pred_df <- tibble(
      p8_1 = rep_param(pars$p8_1), p8_2 = rep_param(pars$p8_2), p8_3 = rep_param(pars$p8_3),
      q_1   = rep_param(pars$q_1),   q_2   = rep_param(pars$q_2),   q_3   = rep_param(pars$q_3),
      q_4   = rep_param(pars$q_4),   q_5   = rep_param(pars$q_5),   q_6   = rep_param(pars$q_6),
      N9    = rep_param(pars$N9),    v9    = rep_param(pars$v9),
      relative_frequency = rf_vec
    ) %>% bind_cols(trial_fixed)
    
    preds <- mapply(
      large_mixed_model,
      p_1 = pred_df$p8_1, p_2 = pred_df$p8_2, p_3 = pred_df$p8_3,
      q_1 = pred_df$q_1, q_2 = pred_df$q_2, q_3 = pred_df$q_3,
      q_4 = pred_df$q_4, q_5 = pred_df$q_5, q_6 = pred_df$q_6,
      N   = pred_df$N9,  v   = pred_df$v9,
      relative_frequency = pred_df$relative_frequency,
      BR  = pred_df$BR,  HR  = pred_df$HR,  FAR = pred_df$FAR,
      SIMPLIFY = TRUE
    )
    
    values <- if (is.matrix(preds)) preds["value", ] else vapply(preds, `[[`, numeric(1), "value")
    
    tibble(
      sample         = i,
      trial          = seq_len(n_trials),
      prediction     = as.numeric(values),
      BR             = pred_df$BR,
      HR             = pred_df$HR,
      FAR            = pred_df$FAR,
      true_posterior = pred_df$true_posterior
    )
  })
  
  results
}

# ---- 3.2 Generate Predictions (7-Parameter MH Model) ----------------------
#' Posterior predictive simulations for MH model
#'
#' @inheritParams generate_ppc_global
#' @param seed Optional random seed for reproducibility
#' @return Posterior predictive samples for MH model

generate_ppc_global_MH <- function(post_params, df_raw, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n_samples <- nrow(post_params)
  n_trials  <- nrow(df_raw)
  
  trial_fixed <- df_raw %>% dplyr::select(BR, HR, FAR, true_posterior, response)
  rep_param <- function(x) rep(as.numeric(x), n_trials)
  
  purrr::map_dfr(seq_len(n_samples), function(i) {
    pars <- post_params[i, ]
    
    pred_df <- tibble::tibble(
      p7_1 = rep_param(pars$p7_1), p7_2 = rep_param(pars$p7_2), p7_3 = rep_param(pars$p7_3),
      p7_4 = rep_param(pars$p7_4), p7_5 = rep_param(pars$p7_5), p7_6 = rep_param(pars$p7_6),
      p7_7 = rep_param(pars$p7_7)
    ) %>% dplyr::bind_cols(trial_fixed)
    
    preds <- mapply(
      FUN  = mixed_heuristic_model,
      BR   = pred_df$BR,
      HR   = pred_df$HR,
      FAR  = pred_df$FAR,
      p_1  = pred_df$p7_1,
      p_2  = pred_df$p7_2,
      p_3  = pred_df$p7_3,
      p_4  = pred_df$p7_4,
      p_5  = pred_df$p7_5,
      p_6  = pred_df$p7_6,
      p_7  = pred_df$p7_7,
      SIMPLIFY = FALSE
    )
    
    values <- vapply(preds, function(df) df[["value"]][1], numeric(1))
    
    tibble::tibble(
      sample         = i,
      trial          = seq_len(n_trials),
      prediction     = as.numeric(values),
      BR             = pred_df$BR,
      HR             = pred_df$HR,
      FAR            = pred_df$FAR,
      true_posterior = pred_df$true_posterior
    )
  })
}


################################################################################
# 4. VISUALIZATION UTILITIES
################################################################################

# ---- 4.1 Reference Lines for Heuristic Benchmarks -------------------------
#' Compute heuristic benchmark values for reference in plots
#'
#' Calculates values for common heuristics: Bayesian Observer (BO),
#' False Complement (FC), Joint Occurrence (JO), Limited Stochasticity (LS),
#' Representativeness (REP), 50% bias, and normative Bayes.
#'
#' @param df          Data frame with BR, HR, FAR, true_posterior
#' @param ref_points  Conditions to compute benchmarks for
#' @param heuristics  Names of heuristics to include
#' @return Long-format tibble with heuristic values per condition

get_reference_lines <- function(df, ref_points,
                                heuristics = c("REP","BO","FC","JO","LS","50%","Bayes")) {
  df %>%
    semi_join(ref_points, by = c("BR","HR","FAR")) %>%
    mutate(
      Bayes = true_posterior,
      REP   = HR,
      BO    = BR,
      FC    = 1 - FAR,
      JO    = BR * HR,
      LS    = HR - FAR,
      `50%` = 0.5,
      across(c(REP,BO,FC,JO,LS,`50%`,Bayes), ~ .x * 100),
      condition = paste0("BR=",BR,", HR=",HR,", FAR=",FAR)
    ) %>%
    distinct(condition, BO, FC, JO, LS, REP, `50%`, Bayes) %>%
    pivot_longer(-condition, names_to = "heuristic", values_to = "value") %>%
    mutate(heuristic = factor(heuristic, levels = heuristics))
}

# ---- 4.2 Annotate Posterior Predictions with Heuristics -------------------
#' Add heuristic values to posterior predictive samples
#'
#' @param ppc_df PPC predictions with BR, HR, FAR
#' @param vlines Reference lines defining conditions to include
#' @return Annotated tibble ready for visualization

clean_ppc <- function(ppc_df, vlines) {
  ppc_df %>%
    mutate(
      Bayes        = true_posterior,
      REP          = HR,
      BO           = BR,
      FC           = 1 - FAR,
      JO           = BR * HR,
      LS           = HR - FAR,
      `50%`        = 0.5,
      across(c(REP,BO,FC,JO,LS,`50%`,Bayes), ~ .x * 100),
      response_pct = prediction * 100,
      condition    = paste0("BR=",BR,", HR=",HR,", FAR=",FAR)
    ) %>%
    filter(condition %in% unique(vlines$condition))
}

# ---- 4.3 Color Palette and Scaling ----------------------------------------
# Okabe-Ito colorblind-friendly palette for heuristic reference lines
okabe_ito <- c(
  BO    = "#E69F00", FC   = "#56B4E9", JO = "#009E73",
  LS    = "#F0E442", REP  = "#0072B2", `50%` = "#CC79A7", Bayes = "#D55E00"
)

heur_labels <- c(
  BO    = "BO", FC = "FC", JO = "JO",
  LS    = "LS", REP = "REP", `50%` = "50 %", Bayes = "Bayes"
)

add_okabe_color <- function() {
  scale_color_manual(name = "Heuristic", values = okabe_ito, labels = heur_labels)
}


################################################################################
# 5. POSTERIOR VALIDATION METRICS
################################################################################

# ---- 5.1 Summary Statistics from Posterior Samples ------------------------
#' Compute summary statistics and regression slopes from posterior predictions
#'
#' Aggregates posterior predictive samples to compute mean, variance, quantiles,
#' and regression slopes against design variables (BR, HR, FAR).
#'
#' @param post_pred           Posterior predictive samples
#' @param column              Response column name
#' @param group_vars          Grouping variables (sample, format)
#' @param predictors          Design variables for regression
#' @param calculate_variance  Include variance decomposition
#' @param variance_column     Column for variance calculation
#' @param var_group_vars      Grouping for variance decomposition
#' @param var_summary_vars    Summary level for variance
#' @return Tibble with summary statistics per posterior sample

summarise_posterior <- function(
    post_pred,
    column          = "prediction",
    group_vars      = c("sample", "format"),
    predictors      = c("BR", "HR", "FAR"),
    calculate_variance = FALSE,
    variance_column    = column,
    var_group_vars     = c("sample", "format", "BR", "HR", "FAR"),
    var_summary_vars   = c("sample", "format")
) {
  
  dt <- data.table::as.data.table(post_pred)
  
  metrics <- compute_all_metrics(
    df         = dt,
    column     = column,
    group_vars = group_vars
  )
  
  slopes  <- compute_SI_by(
    dt         = dt,
    group_vars = group_vars,
    predictors = predictors,
    column     = column
  )
  
  final_tbl <- dplyr::left_join(metrics, slopes, by = group_vars)
  
  if (calculate_variance) {
    var_tbl <- compute_variance_summary(
      post_pred,
      column       = variance_column,
      group_vars   = var_group_vars,
      summary_vars = var_summary_vars
    )
    
    join_keys <- intersect(names(final_tbl), names(var_tbl))
    final_tbl <- dplyr::left_join(final_tbl, var_tbl, by = join_keys)
  }
  
  final_tbl
}

# ---- 5.2 Combine Observed and Model Statistics ----------------------------
#' Prepare long-format table combining observed and model-predicted statistics
#'
#' @param posterior_all  Summary statistics from posterior
#' @param observed_df    Observed summary statistics
#' @param human_dt       Raw human data (unused, for compatibility)
#' @param group_vars     Grouping variables
#' @param id_var         Subject identifier column
#' @return Long-format tibble with type indicator (Observed vs Model)

prepare_stats_long <- function(
    posterior_all,
    observed_df,
    human_dt,
    group_vars = c("sample", "format"),
    id_var = NULL) {
  
  id_var <- id_var %||%
    (c("subject_s", "subject") %>% intersect(names(observed_df)) %>% first())
  if (is.na(id_var)) stop("ID column not found.", call. = FALSE)
  
  indiv_obs <- observed_df |>
    dplyr::select(-dplyr::all_of(id_var))
  
  dplyr::bind_rows(
    dplyr::mutate(indiv_obs ,      type = "Observed"),
    dplyr::mutate(posterior_all,  type = "Model")
  ) |>
    tidyr::pivot_longer(-c(type, dplyr::all_of(group_vars)),
                        names_to = "stat", values_to = "value") |>
    dplyr::select(-dplyr::all_of(group_vars))
}

# ---- 5.3 Distribution Overlap Coefficient ----------------------------------
#' Compute histogram-based overlap coefficient between two distributions
#'
#' Overlap coefficient = sum of minima of normalized histogram bins.
#' Values range from 0 (no overlap) to 1 (perfect overlap).
#'
#' @param x     First distribution
#' @param y     Second distribution
#' @param bins  Number of histogram bins
#' @return Overlap coefficient

compute_hist_overlap <- function(x, y, bins = 30) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  rng <- range(c(x, y), finite = TRUE)
  if (diff(rng) == 0) return(1)
  
  breaks <- seq(rng[1], rng[2], length.out = bins + 1)
  
  hx <- hist(x, breaks = breaks, plot = FALSE, right = FALSE)
  hy <- hist(y, breaks = breaks, plot = FALSE, right = FALSE)
  
  px <- hx$counts / sum(hx$counts)     
  py <- hy$counts / sum(hy$counts)
  
  sum(pmin(px, py))
}

#' Calculate overlap coefficients for all summary statistics
#'
#' @param stats_long  Long-format table with Observed and Model statistics
#' @param bins        Number of bins for histogram overlap
#' @return Table of overlap coefficients per statistic

get_overlap_tbl <- function(stats_long, bins = 30) {
  stats_long |>
    dplyr::filter(!grepl("subject_s", stat)) |>
    dplyr::group_by(stat) |>
    tidyr::nest() |>
    dplyr::mutate(
      overlap = purrr::map_dbl(data, \(df) {
        x <- df |> dplyr::filter(type == "Model")    |> dplyr::pull(value)
        y <- df |> dplyr::filter(type == "Observed") |> dplyr::pull(value)
        compute_hist_overlap(x, y, bins = bins)
      })
    ) |>
    dplyr::select(stat, overlap) |>
    dplyr::arrange(overlap) |>
    data.table::as.data.table()
}

# ---- 5.4 Posterior Predictive P-values ------------------------------------
#' Calculate posterior predictive p-values (Schmidt et al., 2023)
#'
#' Computes one-sided PPP: proportion of posterior predictive samples
#' where statistic ≥ observed value. Extreme values (close to 0 or 1)
#' indicate model misfit.
#'
#' @param df     Combined observed and model statistics
#' @param alpha  Threshold for flagging poor fit
#' @return Table with PPP values and fit indicators per statistic

calculate_ppp <- function(df, alpha = 0.05) {
  df %>%
    group_by(stat) %>%
    summarise(
      T_obs = value[type == "Observed"][1],
      T_rep = list(value[type == "Model"]),
      ppp   = mean(T_rep[[1]] >= T_obs),      
      .groups = "drop"
    ) %>%
    mutate(
      bad_fit   = ppp < alpha,
      direction = if_else(bad_fit, "under", "ok")   
    )
}


################################################################################
# 6. PARAMETER UTILITIES
################################################################################

# ---- 6.1 Aggregate Parameters by Iteration (11-Parameter Model) -----------
#' Extract unique parameter sets from trial-level simulations
#'
#' @param sim_params  Trial-level parameter table with Iteration column
#' @param param_cols  Parameter names to extract
#' @return Iteration-level parameter table

aggregate_sim_params <- function(sim_params, 
                                 param_cols = c(paste0("p8_", 1:3), 
                                                paste0("q_", 1:6), 
                                                "N9", "v9")) {
  sim_params <- as.data.table(sim_params)
  sim_params_agg <- unique(sim_params[, c("Iteration", param_cols), with = FALSE])
  setDT(sim_params_agg)
  return(sim_params_agg)
}

# ---- 6.2 Aggregate Parameters by Iteration (7-Parameter MH Model) ---------
#' Extract unique parameter sets for MH model
#'
#' @inheritParams aggregate_sim_params
#' @return Iteration-level MH parameter table

aggregate_sim_params_MH <- function(sim_params, 
                                 param_cols = c(paste0("p7_", 1:7))) {
  sim_params <- as.data.table(sim_params)
  sim_params_agg <- unique(sim_params[, c("Iteration", param_cols), with = FALSE])
  setDT(sim_params_agg)
  return(sim_params_agg)
}

# ---- 6.3 Prior vs Posterior Visualization ---------------------------------
#' Plot prior and posterior distributions for a single parameter
#'
#' Generates density plots comparing prior distribution (from simulations)
#' against posterior distributions for frequency and probability formats.
#'
#' @param param_name     Parameter to visualize
#' @param prior_sampler  Function to sample from prior
#' @param data           Analysis results containing posteriors_parameter
#' @param x_label        Axis label for parameter
#' @param title          Plot title (optional)
#' @return ggplot object

plot_prior_vs_posterior_param <- function(param_name,
                                          prior_sampler,
                                          data,
                                          x_label = NULL,
                                          title = NULL) {
  if (is.null(x_label)) x_label <- param_name
  
  posterior_freq <- data$posteriors_parameter$frequency[[param_name]]
  posterior_prob <- data$posteriors_parameter$probability[[param_name]]
  prior <- prior_sampler(n_samples)
  
  df <- data.frame(
    value = c(prior, posterior_freq, posterior_prob),
    source = factor(c(
      rep("Prior", length(prior)),
      rep("Posterior (Freq)", length(posterior_freq)),
      rep("Posterior (Prob)", length(posterior_prob))
    ), levels = c("Prior", "Posterior (Freq)", "Posterior (Prob)"))
  )
  
  ggplot(df, aes(x = value, color = source)) +
    geom_density(size = 0.6, adjust = 1.2) +
    scale_color_manual(values = c(
      "Prior" = "#333333",
      "Posterior (Freq)" = "#D55E00",
      "Posterior (Prob)" = "#0072B2"
    )) +
    labs(title = title, x = x_label, y = "Density") +
    theme_classic(base_size = 10) +             
    theme(
      legend.title  = element_blank(),
      legend.position = "right",
      legend.key.height = unit(3, "mm"),
      legend.text   = element_text(size = 8),    
      plot.title    = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title    = element_text(size = 9),
      axis.text     = element_text(size = 8)
    )
}


################################################################################
# 7. INTEGRATED ANALYSIS PIPELINES
################################################################################

# ---- 7.1 Full Pipeline with Overlayed Histograms (11-Parameter Model) ----
#' Complete ABC inference and posterior predictive checking workflow
#'
#' Executes full analysis pipeline:
#'   1. Compute summary statistics from empirical data by format
#'   2. Run ABC rejection sampling to estimate posteriors
#'   3. Generate posterior predictive samples
#'   4. Create overlayed histograms comparing data and predictions
#'
#' @param dataset_name        Label for output messages
#' @param human_data          Trial-level empirical data
#' @param sim_params          Prior parameter samples
#' @param simulations         Simulated summary statistics from prior
#' @param ppc_subject_id      Subject ID for posterior predictive plot
#' @param subject_col         Column name for subject identifier
#' @param ref_points          Stimulus conditions for reference lines
#' @param calculate_variance  Include variance decomposition in summaries
#' @param abc_tol             ABC rejection tolerance
#' @return List with final plot, posterior samples, predictions, summaries

run_analysis_pipeline_overlap <- function(
    dataset_name,
    human_data,
    sim_params,
    simulations,
    ppc_subject_id,
    subject_col,
    ref_points,
    calculate_variance = FALSE,
    abc_tol = 0.01
) {
  message("--- Running Analysis for: ", dataset_name, " ---")
  
  # Compute empirical summary statistics
  human_sum <- compute_all_metrics(
    df = human_data,
    column = "response",
    group_vars = c("format")
  )
  slope_int <- compute_SI_by(
    dt = human_data,
    group_vars = c("format"),
    predictors = c("BR","HR","FAR"),
    column = "response"
  )
  final_human <- left_join(human_sum, slope_int, by = "format")
  
  if (calculate_variance) {
    var_summary <- compute_variance_summary(
      human_data, column = "response",
      group_vars = c("format","BR","HR","FAR"),
      summary_vars = c("format")
    )
    final_human <- left_join(final_human, var_summary, by = "format")
  }
  
  # ABC posterior estimation
  abc_out <- compute_abc_by_format(
    human_summary = final_human,
    sim_params     = sim_params,
    simulations    = simulations,
    tol            = abc_tol
  )
  
  # Generate posterior predictive samples
  subject_filter <- rlang::expr(!!rlang::sym(subject_col) == ppc_subject_id)
  target_data    <- filter(human_data, !!subject_filter)
  
  ppc_prob <- generate_ppc_global(abc_out$post_draws[["probability"]], target_data) %>%
    mutate(format = "probability")
  ppc_freq <- generate_ppc_global(abc_out$post_draws[["frequency"]],   target_data) %>%
    mutate(format = "frequency")
  ppc_all  <- bind_rows(ppc_prob, ppc_freq)
  
  # Prepare visualization data
  ref_lines <- get_reference_lines(human_data, ref_points)

  ppc_clean_p <- clean_ppc(ppc_prob, ref_lines)
  ppc_clean_f <- clean_ppc(ppc_freq, ref_lines)
  vlines_all<-ref_lines
  human_df_prob <- human_data %>% 
    filter(format == "probability") %>%
    mutate(condition = paste0("BR=",BR,", HR=",HR,", FAR=",FAR),
           response_pct = response * 100)%>%
    filter(condition %in% vlines_all$condition)
  
  human_df_freq <- human_data %>% 
    filter(format == "frequency") %>%
    mutate(condition = paste0("BR=",BR,", HR=",HR,", FAR=",FAR),
           response_pct = response * 100)%>%
    filter(condition %in% vlines_all$condition)
  
  model_df_prob <- ppc_clean_p %>% 
    mutate(response_pct = prediction * 100,
           condition = paste0("BR=",BR,", HR=",HR,", FAR=",FAR))%>%
    filter(condition %in% vlines_all$condition)
  
  model_df_freq <- ppc_clean_f  %>% 
    mutate(response_pct = prediction * 100,
           condition = paste0("BR=",BR,", HR=",HR,", FAR=",FAR))%>%
    filter(condition %in% vlines_all$condition)
  
  # Create overlayed histogram plots
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
      coord_cartesian(ylim = c(0,0.15)) +
      labs(title = format_label, x = "Estimates (%)", y = "Density") +
      theme_bw(base_size = 6) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(3, "mm"),
        legend.text = element_text(size = 5),
        axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.x  = element_text(size = 5),
        axis.text.y  = element_text(size = 5),
        strip.text   = element_text(size = 6, face = "bold")
      )
  }
  
  plot_prob_overlay <- plot_overlay(human_df_prob, model_df_prob, ref_lines, "Probability Format")
  plot_freq_overlay <- plot_overlay(human_df_freq, model_df_freq, ref_lines, "Frequency Format")
  
  final_plot <- plot_prob_overlay + plot_freq_overlay +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  message("--- Analysis for ", dataset_name, " complete. ---")
  
  list(plot = final_plot, 
       posteriors_parameter = abc_out$post_draws, 
       posteriors_prediction = ppc_all,
       human_summary = final_human)
}

# ---- 7.2 Full Pipeline with Overlayed Histograms (7-Parameter MH Model) --
#' Complete analysis pipeline for MH model variant
#'
#' @inheritParams run_analysis_pipeline_overlap
#' @return Analysis results for MH model

run_analysis_pipeline_overlap_MH <- function(
    dataset_name,
    human_data,
    sim_params,
    simulations,
    ppc_subject_id,
    subject_col,
    ref_points,
    calculate_variance = FALSE,
    abc_tol = 0.01
) {
  message("--- Running Analysis for: ", dataset_name, " ---")
  
  # Compute empirical summary statistics
  human_sum <- compute_all_metrics(
    df = human_data,
    column = "response",
    group_vars = c("format")
  )
  slope_int <- compute_SI_by(
    dt = human_data,
    group_vars = c("format"),
    predictors = c("BR","HR","FAR"),
    column = "response"
  )
  final_human <- left_join(human_sum, slope_int, by = "format")
  
  if (calculate_variance) {
    var_summary <- compute_variance_summary(
      human_data, column = "response",
      group_vars = c("format","BR","HR","FAR"),
      summary_vars = c("format")
    )
    final_human <- left_join(final_human, var_summary, by = "format")
  }
  
  # ABC posterior estimation for MH model
  abc_out <- compute_abc_by_format_MH(
    human_summary = final_human,
    sim_params     = sim_params,
    simulations    = simulations,
    tol            = abc_tol
  )
  
  # Generate posterior predictive samples
  subject_filter <- rlang::expr(!!rlang::sym(subject_col) == ppc_subject_id)
  target_data    <- filter(human_data, !!subject_filter)
  
  ppc_prob <- generate_ppc_global_MH(abc_out$post_draws[["probability"]], target_data) %>%
    mutate(format = "probability")
  ppc_freq <- generate_ppc_global_MH(abc_out$post_draws[["frequency"]],   target_data) %>%
    mutate(format = "frequency")
  ppc_all  <- bind_rows(ppc_prob, ppc_freq)
  
  # Prepare visualization data
  ref_lines <- get_reference_lines(human_data, ref_points)
  
  ppc_clean_p <- clean_ppc(ppc_prob, ref_lines)
  ppc_clean_f <- clean_ppc(ppc_freq, ref_lines)
  vlines_all<-ref_lines
  human_df_prob <- human_data %>% 
    filter(format == "probability") %>%
    mutate(condition = paste0("BR=",BR,", HR=",HR,", FAR=",FAR),
           response_pct = response * 100)%>%
    filter(condition %in% vlines_all$condition)
  
  human_df_freq <- human_data %>% 
    filter(format == "frequency") %>%
    mutate(condition = paste0("BR=",BR,", HR=",HR,", FAR=",FAR),
           response_pct = response * 100)%>%
    filter(condition %in% vlines_all$condition)
  
  model_df_prob <- ppc_clean_p %>% 
    mutate(response_pct = prediction * 100,
           condition = paste0("BR=",BR,", HR=",HR,", FAR=",FAR))%>%
    filter(condition %in% vlines_all$condition)
  
  model_df_freq <- ppc_clean_f  %>% 
    mutate(response_pct = prediction * 100,
           condition = paste0("BR=",BR,", HR=",HR,", FAR=",FAR))%>%
    filter(condition %in% vlines_all$condition)
  
  # Create overlayed histogram plots
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
  
  plot_prob_overlay <- plot_overlay(human_df_prob, model_df_prob, ref_lines, "Probability Format")
  plot_freq_overlay <- plot_overlay(human_df_freq, model_df_freq, ref_lines, "Frequency Format")
  
  final_plot <- plot_prob_overlay + plot_freq_overlay +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  message("--- Analysis for ", dataset_name, " complete. ---")
  
  list(plot = final_plot, 
       posteriors_parameter = abc_out$post_draws, 
       posteriors_prediction = ppc_all,
       human_summary = final_human)
}


################################################################################
# END OF MODULE
################################################################################

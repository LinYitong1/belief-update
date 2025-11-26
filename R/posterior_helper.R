################################################################################
# POSTERIOR HELPER FUNCTIONS FOR BELIEF UPDATING MODELS (SUBJECT-LEVEL ONLY)
################################################################################
#
# This module implements Approximate Bayesian Computation (ABC) inference and
# posterior predictive checking for cognitive models of belief updating.
#
# This version focuses on SUBJECT-LEVEL analyses:
#   - ABC rejection sampling with neural network adjustment per subject
#   - Posterior predictive simulations per subject
#   - Model validation metrics (overlap coefficients, PPP values) per subject
#   - Integrated subject-level analysis pipeline
#
# Author: Yitong Lin (modified for subject-level only)
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
# 1. ABC POSTERIOR INFERENCE (MODEL-LEVEL, USED FOR SUBJECT-LEVEL)
################################################################################

# ---- 1.1 ABC Rejection Sampling for 11-Parameter Model --------------------
#' Estimate posterior distributions via ABC rejection sampling
#'
#' Implements neural network-adjusted ABC (Blum & François, 2010) with logit
#' transformation for bounded parameters and log transformation for unbounded.
#' This function is used both for group-level and subject-level ABC, depending
#' on which summary statistics are passed in.
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
#' Can be used for subject-level ABC by passing subject-specific summary statistics.
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
# 2. POSTERIOR PREDICTIVE SIMULATIONS
################################################################################

# ---- 2.1 Generate Predictions (11-Parameter Model) ------------------------
#' Generate posterior predictive samples for trial-level predictions
#'
#' For each posterior parameter draw, simulate model predictions across all
#' trials including stochastic sample size effects. This is used both at the
#' group-level and subject-level depending on the input df_raw.
#'
#' @param post_params Posterior parameter samples (rows = draws)
#' @param df_raw      Trial-level data with BR, HR, FAR, true_posterior, response
#' @return Tibble with predictions for each posterior sample × trial combination

generate_ppc_global <- function(post_params, df_raw) {
  n_samples <- nrow(post_params)
  n_trials  <- nrow(df_raw)
  
  trial_fixed <- df_raw %>% dplyr::select(BR, HR, FAR, true_posterior, response)
  rep_param <- function(x) rep(as.numeric(x), n_trials)
  
  results <- purrr::map_dfr(seq_len(n_samples), function(i) {
    pars <- post_params[i, ]
    rf_vec <- simulate_and_mutate(df_raw, pars$N9)$relative_frequency
    
    pred_df <- tibble::tibble(
      p8_1 = rep_param(pars$p8_1), p8_2 = rep_param(pars$p8_2), p8_3 = rep_param(pars$p8_3),
      q_1   = rep_param(pars$q_1),   q_2   = rep_param(pars$q_2),   q_3   = rep_param(pars$q_3),
      q_4   = rep_param(pars$q_4),   q_5   = rep_param(pars$q_5),   q_6   = rep_param(pars$q_6),
      N9    = rep_param(pars$N9),    v9    = rep_param(pars$v9),
      relative_frequency = rf_vec
    ) %>% dplyr::bind_cols(trial_fixed)
    
    preds <- mapply(
      large_mixed_model,
      p_1 = pred_df$p8_1, p_2 = pred_df$p8_2, p_3 = pred_df$p8_3,
      q_1 = pred_df$q_1,  q_2 = pred_df$q_2,  q_3 = pred_df$q_3,
      q_4 = pred_df$q_4,  q_5 = pred_df$q_5,  q_6 = pred_df$q_6,
      N   = pred_df$N9,   v   = pred_df$v9,
      relative_frequency = pred_df$relative_frequency,
      BR  = pred_df$BR,   HR  = pred_df$HR,   FAR = pred_df$FAR,
      SIMPLIFY = TRUE
    )
    
    values <- if (is.matrix(preds)) preds["value", ] else vapply(preds, `[[`, numeric(1), "value")
    
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
  
  results
}

# ---- 2.2 Generate Predictions (7-Parameter MH Model) ----------------------
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
# 3. POSTERIOR VALIDATION METRICS (GENERIC, USED AT SUBJECT LEVEL)
################################################################################

# ---- 3.1 Summary Statistics from Posterior Samples ------------------------
#' Compute summary statistics and regression slopes from posterior predictions
#'
#' Aggregates posterior predictive samples to compute mean, variance, quantiles,
#' and regression slopes against design variables (BR, HR, FAR).
#'
#' This function is generic: by changing `group_vars` you can use it at
#' subject-level (e.g., group_vars = c("sample", "subject")).
#'
#' @param post_pred           Posterior predictive samples
#' @param column              Response column name
#' @param group_vars          Grouping variables (e.g., sample, subject)
#' @param predictors          Design variables for regression
#' @param calculate_variance  Include variance decomposition
#' @param variance_column     Column for variance calculation
#' @param var_group_vars      Grouping for variance decomposition
#' @param var_summary_vars    Summary level for variance
#' @return Tibble with summary statistics per posterior sample

summarise_posterior <- function(
    post_pred,
    column          = "prediction",
    group_vars      = c("sample", "subject"),
    predictors      = c("BR", "HR", "FAR"),
    calculate_variance = FALSE,
    variance_column    = column,
    var_group_vars     = c("sample", "subject", "BR", "HR", "FAR"),
    var_summary_vars   = c("sample", "subject")
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
      data.table::as.data.table(post_pred),  # ← 改这里
      column       = variance_column,
      group_vars   = var_group_vars,
      summary_vars = var_summary_vars
    )
    
    join_keys <- intersect(names(final_tbl), names(var_tbl))
    final_tbl <- dplyr::left_join(final_tbl, var_tbl, by = join_keys)
  }
  
  
  final_tbl
}


# ---- 3.2 Distribution Overlap Coefficient ---------------------------------
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

# ---- 3.3 Posterior Predictive P-values ------------------------------------
#' Calculate posterior predictive p-values (Schmidt et al., 2023)
#'
#' Computes one-sided PPP: proportion of posterior predictive samples
#' where statistic ≥ observed value. Extreme values (close to 0 or 1)
#' indicate model misfit.
#'
#' This is a generic function; at the subject level, it is typically applied
#' within each subject (see compute_ppp_by_subject()).
#'
#' @param df     Combined observed and model statistics
#' @param alpha  Threshold for flagging poor fit
#' @return Table with PPP values and fit indicators per statistic

calculate_ppp <- function(df, alpha = 0.05) {
  df %>%
    dplyr::group_by(stat) %>%
    dplyr::summarise(
      T_obs = unique(value[type == "Observed"]),
      ppp   = mean(value[type == "Model"] >= T_obs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      bad_fit   = ppp < alpha,
      direction = dplyr::if_else(bad_fit, "under", "ok")
    ) 
}

# ppp<alpha/2 or ppp>1-a/pha/2
# 
################################################################################
# 4. SUBJECT-LEVEL ABC AND PPP
################################################################################

# ---- 4.1 ABC by Subject (11-Parameter Model) ------------------------------
#' Compute ABC posteriors separately for each subject
#'
#' @param human_summary Summary statistics from empirical data with subject column
#' @param sim_params    Parameter samples from prior
#' @param simulations   Simulated summary statistics
#' @param tol           ABC tolerance threshold
#' @param subject_col   Name of subject ID column (default: "subject")
#' @return List containing abc results and posterior draws for each subject

compute_abc_by_subject <- function(human_summary,
                                   sim_params,
                                   simulations,
                                   tol = 0.05,
                                   subject_col = "subject") {
  human_summary <- tibble::as_tibble(human_summary)
  sim_params    <- tibble::as_tibble(sim_params)
  simulations   <- tibble::as_tibble(simulations)
  
  stopifnot(subject_col %in% names(human_summary))
  
  subjects <- as.character(unique(human_summary[[subject_col]]))
  

  SUM_COLS_full <- setdiff(names(simulations), c("Iteration", "model"))
  SUM_COLS      <- intersect(SUM_COLS_full, names(human_summary))
  

  simulations_sub <- simulations[, c("Iteration", "model", SUM_COLS), drop = FALSE]
  
  abc_res    <- setNames(vector("list", length(subjects)), subjects)
  post_draws <- setNames(vector("list", length(subjects)), subjects)
  
  for (sbj in subjects) {
    obs_stats <- human_summary %>%
      dplyr::filter(.data[[subject_col]] == sbj) %>%
      dplyr::select(dplyr::all_of(SUM_COLS)) %>%
      unlist(use.names = FALSE)
    

    abc_obj <- abc_posterior(obs_stats, sim_params, simulations_sub, tol)
    draws   <- as.data.frame(abc_obj$unadj.values)
    draws[[subject_col]] <- sbj
    
    abc_res[[sbj]]    <- abc_obj
    post_draws[[sbj]] <- draws
  }
  
  list(abc_res = abc_res, post_draws = post_draws, tol = tol)
}


# ---- 4.2 Run ABC for All Subjects -----------------------------------------
#' Full subject-level ABC workflow (summary + ABC by subject)
#'
#' @param human_data         Trial-level empirical data
#' @param sim_params         Prior parameter samples
#' @param simulations        Simulated summary statistics
#' @param subject_col        Subject ID column name
#' @param calculate_variance Whether to compute variance summaries
#' @param abc_tol            ABC tolerance
#' @return List with ABC results and human summary statistics

run_abc_all_subjects <- function(human_data,
                                 sim_params,
                                 simulations,
                                 subject_col = "subject",
                                 calculate_variance = FALSE,
                                 abc_tol = 0.05) {
  # 1. Empirical summary statistics per subject
  human_sum <- compute_all_metrics(
    df         = human_data,
    column     = "response",
    group_vars = c(subject_col)
  )
  
  slope_int <- compute_SI_by(
    dt         = human_data,
    group_vars = c(subject_col),
    predictors = c("BR", "HR", "FAR"),
    column     = "response"
  )
  
  final_human <- dplyr::left_join(human_sum, slope_int, by = subject_col)
  
  if (calculate_variance) {
    var_summary <- compute_variance_summary(
      data.table::as.data.table(human_data),           # 确保是 data.table
      column       = "response",
      group_vars   = c(subject_col, "BR", "HR", "FAR"),
      summary_vars = c(subject_col)
    )
    final_human <- dplyr::left_join(final_human, var_summary, by = subject_col)
  }
  
  # 2. ABC per subject
  abc_out <- compute_abc_by_subject(
    human_summary = final_human,
    sim_params    = sim_params,
    simulations   = simulations,
    tol           = abc_tol,
    subject_col   = subject_col
  )
  
  list(
    abc_out       = abc_out,
    human_summary = final_human
  )
}

# ---- 4.3 Long-format stats for subject-level PPP --------------------------
#' Prepare long-format table combining observed and model stats by subject
#'
#' @param posterior_all Summary statistics from posterior (per sample × subject)
#' @param observed_df   Observed summary statistics per subject
#' @param subject_col   Subject ID column name
#' @return Long-format tibble with columns: subject, stat, value, type

prepare_stats_long_subject <- function(posterior_all,
                                       observed_df,
                                       subject_col = "subject") {
  posterior_long <- posterior_all %>%
    tidyr::pivot_longer(
      cols      = -c(sample, !!rlang::sym(subject_col)),
      names_to  = "stat",
      values_to = "value"
    ) %>%
    dplyr::mutate(type = "Model")
  
  observed_long <- observed_df %>%
    tidyr::pivot_longer(
      cols      = -!!rlang::sym(subject_col),
      names_to  = "stat",
      values_to = "value"
    ) %>%
    dplyr::mutate(type = "Observed")
  
  dplyr::bind_rows(observed_long, posterior_long)
}

# ---- 4.4 Subject-wise PPP -------------------------------------------------
#' Compute PPP per subject and summary statistic
#'
#' Wrapper around calculate_ppp(), applied within each subject.
#'
#' @param posterior_summary Posterior summary statistics (per sample × subject)
#' @param observed_df       Observed summary statistics per subject
#' @param subject_col       Subject ID column name
#' @param alpha             PPP threshold
#' @return Tibble with PPP per subject × stat

compute_ppp_by_subject <- function(posterior_summary,
                                   observed_df,
                                   subject_col = "subject",
                                   alpha = 0.05) {
  stats_long_subj <- prepare_stats_long_subject(
    posterior_all = posterior_summary,
    observed_df   = observed_df,
    subject_col   = subject_col
  )
  
  ppp_by_subject <- stats_long_subj %>%
    dplyr::group_by(.data[[subject_col]]) %>%
    dplyr::group_modify(~ calculate_ppp(.x, alpha = alpha)) %>%
    dplyr::ungroup()
  
  ppp_by_subject
}


################################################################################
# 5. INTEGRATED SUBJECT-LEVEL ANALYSIS PIPELINE
################################################################################

# ---- 5.1 Full Subject-Level Validation Pipeline ---------------------------
#' Subject-level ABC + posterior predictive checks + overlap + PPP
#'
#' This function performs a complete subject-level analysis:
#'   1. Run ABC per subject to obtain individual posterior parameters
#'   2. Generate trial-level posterior predictive samples per subject
#'   3. Summarise PPC into sample × subject summary statistics
#'   4. Compute subject-wise overlap coefficients (Model vs Observed)
#'   5. Compute subject-wise PPP values for each summary statistic
#'   6. Provide overall summaries of PPP and overlap distributions
#'
#' @param dataset_name  Label for logging
#' @param human_data    Trial-level empirical data
#' @param sim_params    Prior parameter samples
#' @param simulations   Simulated summary statistics
#' @param subject_col   Subject ID column (default: "subject")
#' @param abc_tol       ABC tolerance
#' @param bins_overlap  Number of bins for histogram overlap
#' @param ppp_alpha     PPP threshold for bad fit
#' @return List containing ABC results, PPC, summaries, PPP, and overlap metrics

run_subject_level_validation <- function(
    dataset_name,
    human_data,
    sim_params,
    simulations,
    subject_col   = "subject",
    abc_tol       = 0.05,
    bins_overlap  = 30,
    ppp_alpha     = 0.05,
    calculate_variance = FALSE
) {
  message("--- Subject-level ABC & validation: ", dataset_name, " ---")
  
  # 1. Run ABC per subject
  abc_indiv <- run_abc_all_subjects(
    human_data         = human_data,
    sim_params         = sim_params,
    simulations        = simulations,
    subject_col        = subject_col,
    calculate_variance = calculate_variance,   # ← 传下去
    abc_tol            = abc_tol
  )
  observed_indiv <- abc_indiv$human_summary
  
  # 2. Generate subject-level posterior predictive samples (trial-level)
  subjects <- sort(unique(human_data[[subject_col]]))
  
  post_pred_indiv <- purrr::map_dfr(subjects, function(sbj) {
    target_data   <- dplyr::filter(human_data, .data[[subject_col]] == sbj)
    draws_subject <- abc_indiv$abc_out$post_draws[[as.character(sbj)]]
    if (is.null(draws_subject) || nrow(draws_subject) == 0) return(NULL)
    
    generate_ppc_global(draws_subject, target_data) %>%
      dplyr::mutate(!!subject_col := sbj)
  })
  
  # 3. Summarise PPC into sample × subject summary statistics
  posterior_indiv <- summarise_posterior(
    post_pred          = post_pred_indiv,
    column             = "prediction",
    group_vars         = c("sample", subject_col),
    predictors         = c("BR", "HR", "FAR"),
    calculate_variance = calculate_variance 
  )
  
  # 4. Combine observed vs model statistics (long format)
  stats_long_indiv <- prepare_stats_long_subject(
    posterior_all = posterior_indiv,
    observed_df   = observed_indiv,
    subject_col   = subject_col
  )
  
  # 5. Subject-level PPP
  ppp_by_subject <- stats_long_indiv %>%
    dplyr::group_by(!!rlang::sym(subject_col)) %>%
    dplyr::group_modify(~ calculate_ppp(.x, alpha = ppp_alpha)) %>%
    dplyr::ungroup()
  
  # 6. Overall PPP summary (across all subject × stat combinations)
  ppp_summary_overall <- ppp_by_subject %>%
    dplyr::summarise(
      n_total      = dplyr::n(),
      median_ppp   = median(ppp, na.rm = TRUE),
      q25_ppp      = quantile(ppp, 0.25, na.rm = TRUE),
      q75_ppp      = quantile(ppp, 0.75, na.rm = TRUE),
      prop_extreme = mean(ppp < ppp_alpha, na.rm = TRUE)
    )
  
  
  message("--- Done subject-level ABC & validation: ", dataset_name, " ---")
  
  list(
    abc_indiv               = abc_indiv,
    observed_indiv          = observed_indiv,
    post_pred_indiv         = post_pred_indiv,
    posterior_indiv         = posterior_indiv,
    stats_long_indiv        = stats_long_indiv,
    ppp_by_subject          = ppp_by_subject,
    ppp_summary_overall     = ppp_summary_overall
  )
}

# -----------------------------------------------------------------------------
# Plot prior vs posterior for one parameter (pooled across subjects)
# -----------------------------------------------------------------------------
plot_param_population <- function(param_name, sim_params, post_subject_means) {
  prior_vals    <- sim_params[[param_name]]
  post_means    <- post_subject_means[[param_name]]
  
  df <- data.frame(
    value  = c(prior_vals, post_means),
    source = factor(
      c(rep("Prior", length(prior_vals)),
        rep("Posterior (subject means)", length(post_means))),
      levels = c("Prior", "Posterior (subject means)")
    )
  )
  
  ggplot(df, aes(x = value, color = source)) +
    geom_density(size = 0.6) +
    scale_color_manual(values = c(
      "Prior"                      = "black",
      "Posterior (subject means)"  = "#D55E00"
    )) +
    labs(
      x = param_name,
      y = "Density",
      color = ""
    ) +
    theme_classic()
}

plot_human_vs_posterior_pooled <- function(
    human_data,
    post_pred_indiv,
    binwidth = 0.05
) {
  human_all <- human_data %>%
    dplyr::mutate(source = "Human")
  
  model_all <- post_pred_indiv %>%
    dplyr::mutate(
      source   = "Model",
      response = prediction
    )
  
  combined <- dplyr::bind_rows(
    human_all %>% dplyr::select(response, source),
    model_all %>% dplyr::select(response, source)
  )
  
  ggplot(combined, aes(x = response, fill = source, color = source)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = binwidth,
      alpha    = 0.6,
      position = "identity"
    ) +
    scale_fill_manual(values = c("Human" = "#3c5488", "Model" = "#f39b7f")) +
    scale_color_manual(values = c("Human" = "#3c5488", "Model" = "#f39b7f")) +
    labs(
      title = "Human vs model posterior predictive (pooled across subjects)",
      x     = "Response / prediction",
      y     = "Density",
      fill  = "",
      color = ""
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position   = "bottom",
      legend.key.height = unit(3, "mm")
    )
}

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

################################################################################
# 8. SUBJECT-LEVEL ABC & PPP FOR MH MODEL (7 PARAMS)
################################################################################

# ---- 8.1 Subject-level posterior summaries (MH) ----------------------------

summarise_posterior_subject_MH <- function(
    post_pred,
    column          = "prediction",
    subject_col     = "subject",
    predictors      = c("BR", "HR", "FAR"),
    calculate_variance = FALSE,
    variance_column    = column
) {
  dt <- data.table::as.data.table(post_pred)
  
  group_vars <- c("sample", subject_col)
  
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
      dt,
      column       = variance_column,
      group_vars   = c("sample", subject_col, "BR", "HR", "FAR"),
      summary_vars = c("sample", subject_col)
    )
    
    join_keys <- intersect(names(final_tbl), names(var_tbl))
    final_tbl <- dplyr::left_join(final_tbl, var_tbl, by = join_keys)
  }
  
  final_tbl
}

# ---- 8.2 Prepare long-format stats by subject (for PPP) --------------------
prepare_stats_long_subject_MH <- function(posterior_all,
                                          observed_df,
                                          subject_col = "subject") {
  posterior_long <- posterior_all %>%
    tidyr::pivot_longer(
      cols      = -c(sample, !!rlang::sym(subject_col)),
      names_to  = "stat",
      values_to = "value"
    ) %>%
    dplyr::mutate(type = "Model")
  
  observed_long <- observed_df %>%
    tidyr::pivot_longer(
      cols      = -!!rlang::sym(subject_col),
      names_to  = "stat",
      values_to = "value"
    ) %>%
    dplyr::mutate(type = "Observed")
  
  dplyr::bind_rows(observed_long, posterior_long)
}

# ---- 8.3 PPP by subject ---------------------------
compute_ppp_by_subject_MH <- function(posterior_summary,
                                      observed_df,
                                      subject_col = "subject",
                                      alpha = 0.05) {
  stats_long_subj <- prepare_stats_long_subject_MH(
    posterior_all = posterior_summary,
    observed_df   = observed_df,
    subject_col   = subject_col
  )
  
  stats_long_subj %>%
    dplyr::group_by(.data[[subject_col]]) %>%
    dplyr::group_modify(~ calculate_ppp(.x, alpha = alpha)) %>%
    dplyr::ungroup()
}

# ---- 8.4 ABC by subject for MH  ----------------------
# human_summary
compute_abc_by_subject_MH <- function(human_summary,
                                      sim_params,
                                      simulations,
                                      tol = 0.05,
                                      subject_col = "subject") {
  human_summary <- tibble::as_tibble(human_summary)
  sim_params    <- tibble::as_tibble(sim_params)
  simulations   <- tibble::as_tibble(simulations)
  
  stopifnot(subject_col %in% names(human_summary))
  
  subjects <- as.character(unique(human_summary[[subject_col]]))
  
  SUM_COLS_full <- setdiff(names(simulations), c("Iteration", "model"))
  SUM_COLS      <- intersect(SUM_COLS_full, names(human_summary))
  
  simulations_sub <- simulations[, c("Iteration", "model", SUM_COLS), drop = FALSE]
  
  abc_res    <- setNames(vector("list", length(subjects)), subjects)
  post_draws <- setNames(vector("list", length(subjects)), subjects)
  
  for (sbj in subjects) {
    obs_stats <- human_summary %>%
      dplyr::filter(.data[[subject_col]] == sbj) %>%
      dplyr::select(dplyr::all_of(SUM_COLS)) %>%
      unlist(use.names = FALSE)
    
    abc_obj <- abc_posterior_MH(obs_stats, sim_params, simulations_sub, tol)
    draws   <- as.data.frame(abc_obj$unadj.values)
    draws[[subject_col]] <- sbj
    
    abc_res[[sbj]]    <- abc_obj
    post_draws[[sbj]] <- draws
  }
  
  list(abc_res = abc_res, post_draws = post_draws, tol = tol)
}

# ---- 8.5 Convenience wrapper: run ABC for all subjects (MH) ----------------
run_abc_all_subjects_MH <- function(human_data,
                                    sim_params,
                                    simulations,
                                    subject_col        = "subject",
                                    calculate_variance = FALSE,
                                    abc_tol            = 0.05) {
  # 1) empirical summary stats per subject
  human_sum <- compute_all_metrics(
    df         = human_data,
    column     = "response",
    group_vars = c(subject_col)
  )
  
  slope_int <- compute_SI_by(
    dt         = human_data,
    group_vars = c(subject_col),
    predictors = c("BR", "HR", "FAR"),
    column     = "response"
  )
  
  final_human <- dplyr::left_join(human_sum, slope_int, by = subject_col)
  
  if (calculate_variance) {
    var_summary <- compute_variance_summary(
      human_data,
      column       = "response",
      group_vars   = c(subject_col, "BR", "HR", "FAR"),
      summary_vars = c(subject_col)
    )
    final_human <- dplyr::left_join(final_human, var_summary, by = subject_col)
  }
  
  # 2) ABC per subject (MH)
  abc_out <- compute_abc_by_subject_MH(
    human_summary = final_human,
    sim_params    = sim_params,
    simulations   = simulations,
    tol           = abc_tol,
    subject_col   = subject_col
  )
  
  list(
    abc_out       = abc_out,
    human_summary = final_human
  )
}

# ---- 8.6 Full subject-level MH validation pipeline -------------------------
run_subject_level_MH_validation <- function(
    dataset_name,
    human_data,
    sim_params,
    simulations,
    subject_col   = "subject",
    abc_tol       = 0.05,
    ppp_alpha     = 0.05,
    calculate_variance = FALSE
) {
  message("--- Subject-level MH ABC & validation: ", dataset_name, " ---")

  abc_indiv <- run_abc_all_subjects_MH(
    human_data         = human_data,
    sim_params         = sim_params,
    simulations        = simulations,
    subject_col        = subject_col,
    calculate_variance = calculate_variance,
    abc_tol            = abc_tol
  )
  observed_indiv <- abc_indiv$human_summary

  subjects <- sort(unique(human_data[[subject_col]]))
  
  post_pred_indiv <- purrr::map_dfr(subjects, function(sbj) {
    target_data   <- dplyr::filter(human_data, .data[[subject_col]] == sbj)
    draws_subject <- abc_indiv$abc_out$post_draws[[as.character(sbj)]]
    if (is.null(draws_subject) || nrow(draws_subject) == 0) return(NULL)
    
    ppc <- generate_ppc_global_MH(draws_subject, target_data)
    
    n_samples <- nrow(draws_subject)
    ppc %>%
      dplyr::mutate(
        !!subject_col := sbj,
        format        = rep(target_data$format, times = n_samples)
      )
  })
  
  # 3. Summarise PPC into sample × subject stats
  posterior_indiv <- summarise_posterior_subject_MH(
    post_pred          = post_pred_indiv,
    column             = "prediction",
    subject_col        = subject_col,
    predictors         = c("BR", "HR", "FAR"),
    calculate_variance = calculate_variance
  )
  
  # 4. PPP per subject × stat
  ppp_by_subject <- compute_ppp_by_subject_MH(
    posterior_summary = posterior_indiv,
    observed_df       = observed_indiv,
    subject_col       = subject_col,
    alpha             = ppp_alpha
  )
  
  # 5. Overall PPP summary
  ppp_summary_overall <- ppp_by_subject %>%
    dplyr::summarise(
      n_total      = dplyr::n(),
      median_ppp   = median(ppp, na.rm = TRUE),
      q25_ppp      = quantile(ppp, 0.25, na.rm = TRUE),
      q75_ppp      = quantile(ppp, 0.75, na.rm = TRUE),
      prop_extreme = mean(ppp < ppp_alpha, na.rm = TRUE)
    )
  
  message("--- Done subject-level MH ABC & validation: ", dataset_name, " ---")
  
  list(
    abc_indiv       = abc_indiv,
    observed_indiv  = observed_indiv,
    post_pred_indiv = post_pred_indiv,
    posterior_indiv = posterior_indiv,
    ppp_by_subject  = ppp_by_subject,
    ppp_summary_overall = ppp_summary_overall
  )
}

################################################################################
# END OF MODULE
################################################################################

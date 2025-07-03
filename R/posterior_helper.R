###############################################################################
# Bayesian Cognition Helper Functions                                         #
# Approximate Bayesian Computation & Posterior Predictive Checks               #
# Author: <your-name>                                                          #
# Date:   2025-07-01                                                           #
###############################################################################

# -----------------------------------------------------------------------------
# 0. Load Required Packages
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)  # Fast data manipulation
  library(dplyr)       # Data wrangling
  library(purrr)       # Functional programming tools
  library(abc)         # Approximate Bayesian Computation
  library(tidyr)       # Tidy data utilities
  library(ggplot2)  
  library(patchwork)# Data visualization
})

# -----------------------------------------------------------------------------
# 1. ABC Posterior Estimation (Rejection Sampling)
# -----------------------------------------------------------------------------
#' Estimate posterior samples via ABC rejection
#'
#' @param obs_stats   Numeric vector of observed summary statistics (order must match `simulations`).
#' @param sim_params  Data frame of simulated parameters (columns p8_1,p8_2,p8_3,q_1..q_6,N9,v9).
#' @param simulations Data frame of simulated summary statistics (excluding Iteration, model).
#' @param tol         Tolerance fraction for ABC (default 0.01).
#' @return            An 'abc' class object containing posterior draws.
abc_posterior <- function(obs_stats,
                          sim_params,
                          simulations,
                          tol = 0.01) {
  # Parameter and summary-statistic columns
  PARAM_COLS <- c(paste0("p8_",1:3), paste0("q_",1:6), "N9", "v9")
  SUM_COLS   <- setdiff(names(simulations), c("Iteration","model"))
  
  # Input validation
  stopifnot(is.numeric(obs_stats), length(obs_stats) == length(SUM_COLS))
  stopifnot(all(PARAM_COLS %in% names(sim_params)), all(SUM_COLS %in% names(simulations)))
  
  sim_params  <- as.data.table(sim_params)
  simulations <- as.data.table(simulations)
  
  abc::abc(
    target  = obs_stats,
    param   = sim_params[, ..PARAM_COLS],
    sumstat = simulations[, ..SUM_COLS],
    method  = "rejection",
    tol     = tol
  )
}

# -----------------------------------------------------------------------------
# 2. ABC by Stimulus Format
# -----------------------------------------------------------------------------
#' Compute ABC posteriors separately for each format
#'
#' @param human_summary Tibble with observed summary stats and 'format'.
#' @param sim_params    Tibble of parameter draws including 'format'.
#' @param simulations   Tibble of simulated stats including 'format'.
#' @param tol           Tolerance for rejection (default 0.01).
#' @return              List with 'abc_res' (per-format abc objects) and
#'                      'post_draws' (raw posterior samples).
compute_abc_by_format <- function(human_summary,
                                  sim_params,
                                  simulations,
                                  tol = 0.01) {
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
    draws$subject <- "global"
    
    abc_res[[fmt]]    <- abc_obj
    post_draws[[fmt]] <- draws
  }
  
  list(abc_res = abc_res, post_draws = post_draws, tol = tol)
}

# -----------------------------------------------------------------------------
# 3. Posterior Predictive Checks (Global)
# -----------------------------------------------------------------------------
#' Generate posterior predictive distributions for trials
#'
#' @param post_params Data frame of posterior parameter draws.
#' @param df_raw      Raw trial data with BR, HR, FAR, true_posterior, response.
#' @return            Tibble of sample × trial predictions and covariates.
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

# -----------------------------------------------------------------------------
# 4. Reference Lines Generator
# -----------------------------------------------------------------------------
#' Create benchmark heuristic lines for plotting
#'
#' @param df          Data frame of observations or predictions.
#' @param ref_points  Data frame of conditions (BR, HR, FAR).
#' @param heuristics  Character vector of heuristic names.
#' @return            Tibble of heuristic values per condition.
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

# -----------------------------------------------------------------------------
# 5. Clean & Annotate PPC Results
# -----------------------------------------------------------------------------
#' Annotate PPC data with heuristic values
#'
#' @param ppc_df   PPC predictions data frame (with BR, HR, FAR).
#' @param vlines   Reference lines tibble including 'condition'.
#' @return         Filtered and annotated tibble for plotting.
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

# -----------------------------------------------------------------------------
# 6. Plotting Helpers
# -----------------------------------------------------------------------------
# Okabe–Ito colorblind-friendly palette
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

#' Plot distribution with heuristic reference lines
#'
#' @param data Data frame with 'response_pct' and 'condition'.
#' @param vlines Reference lines tibble.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
plot_distribution <- function(data, vlines, xlab, ylab) {
  ggplot(data, aes(x = response_pct)) +
    geom_histogram( aes(y = after_stat(density)), binwidth = 3, fill = "#5b5e6e") +
    geom_vline(data = vlines, aes(xintercept = value, colour = heuristic),
               linetype = "dashed", size = .6) +
    add_okabe_color() +
    facet_wrap(~ condition, ncol = 1, scales = "free_y") +
    coord_cartesian(ylim = c(0,0.13)) +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 5) +
    theme(
      panel.grid.major.y = element_line(size = 0.3, colour = "grey85"),
      panel.grid.minor   = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold"),
      legend.position    = "bottom",
      legend.key.height  = unit(4, "mm"),
      legend.text        = element_text(size = 6),
      axis.title.x       = element_text(size = 5),  
      axis.title.y       = element_text(size = 5),  
      axis.text.x        = element_text(size = 4.5),
      axis.text.y        = element_text(size = 4.5)
    )
}

# -----------------------------------------------------------------------------
# 7. Human Estimates by Format
# -----------------------------------------------------------------------------
#' Plot human responses by format using a custom function
#'
#' @param df_raw     Raw human data with 'format', BR, HR, FAR, true_posterior, response.
#' @param formats    Vector of format labels.
#' @param vlines_all Reference lines tibble.
#' @param plot_fn    Plotting function (e.g., plot_distribution).
#' @return           Named list of ggplot objects.
plot_human_by_format <- function(df_raw, formats, vlines_all, plot_fn) {
  clean_human <- function(df, fmt) {
    df %>%
      filter(format == fmt) %>%
      mutate(
        condition = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
        Bayes = true_posterior, REP = HR, BO = BR, FC = 1 - FAR,
        JO = BR * HR, LS = HR - FAR, `50%` = 0.5,
        across(c(REP,BO,FC,JO,LS,`50%`,Bayes), ~ .x * 100),
        response_pct = response * 100
      ) %>%
      filter(condition %in% vlines_all$condition)
  }
  
  plots <- map(formats, function(fmt) {
    plt <- clean_human(df_raw, fmt) %>%
      plot_fn(vlines_all, paste("Human Estimates (%) -", fmt), "Density") 
    plt
  })
  names(plots) <- formats
  plots
}

# -----------------------------------------------------------------------------
# 8. Full Analysis Pipeline
# -----------------------------------------------------------------------------
#' Run complete pipeline: summaries, ABC, PPC, and plotting
#'
#' @param dataset_name   Label for dataset.
#' @param human_data     Raw human trial data.
#' @param sim_params     Simulation parameter draws.
#' @param simulations    Simulation summary statistics.
#' @param ppc_subject_id ID of subject for PPC.
#' @param subject_col    Column name for subject filter.
#' @param ref_points     Data frame of benchmark conditions.
#' @param calculate_variance Logical to include variance in summaries.
#' @param abc_tol        ABC tolerance (default 0.01).
#' @return               List with final plot, posterior draws, and summaries.
run_analysis_pipeline <- function(
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
  
  # 1. Human summary by format
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
  
  # 2. ABC estimation
  abc_out <- compute_abc_by_format(
    human_summary = final_human,
    sim_params     = sim_params,
    simulations    = simulations,
    tol            = abc_tol
  )
  
  # 3. Posterior predictive checks
  subject_filter <- rlang::expr(!!rlang::sym(subject_col) == ppc_subject_id)
  target_data    <- filter(human_data, !!subject_filter)
  
  ppc_prob <- generate_ppc_global(abc_out$post_draws[["probability"]], target_data) %>%
    mutate(format = "probability")
  ppc_freq <- generate_ppc_global(abc_out$post_draws[["frequency"]],   target_data) %>%
    mutate(format = "frequency")
  ppc_all  <- bind_rows(ppc_prob, ppc_freq)
  
  # 4. Reference lines & cleaning
  ref_lines <- get_reference_lines(human_data, ref_points)
  ppc_clean_p <- clean_ppc(ppc_prob, ref_lines)
  ppc_clean_f <- clean_ppc(ppc_freq, ref_lines)
  
  # 5. Plotting
  human_plots <- plot_human_by_format(
    df_raw  = human_data,
    formats = unique(human_data$format),
    vlines_all = ref_lines,
    plot_fn = plot_distribution
  )
  model_p <- plot_distribution(ppc_clean_p, ref_lines,
                               "Posterior Prediction (%) - Probability","Density")
  model_f <- plot_distribution(ppc_clean_f, ref_lines,
                               "Posterior Prediction (%) - Frequency","Density")
  
  final_plot <- (human_plots[[1]] | human_plots[[2]] | model_p | model_f) +
    plot_layout(ncol = 4, guides = "collect") &
    theme(
      plot.title        = element_text(size = 5, hjust = 0.5),
      legend.position   = "bottom",
      legend.key.height = unit(4, "mm"),
      legend.text       = element_text(size = 5),
      strip.text        = element_text(face = "bold")
    )
  
  message("--- Analysis for ", dataset_name, " complete. ---")
  list(plot = final_plot, 
       posteriors_parameter = abc_out$post_draws, 
       posteriors_prediction = ppc_all,
       human_summary = final_human)
}

# -----------------------------------------------------------------------------
# 9. Individual‑level Posterior vs. Observed Diagnostics
# -----------------------------------------------------------------------------

# The following helpers compute summary statistics per posterior draw, merge
# them with the observed participant‑level statistics, and visualise the extent
# to which the model reproduces individual variability.

# 9.1 Metric aggregation + regression slopes
# ── Summarise posterior metrics, regression slopes, and (optional) variance ──
# ── Summarise posterior metrics, slopes, and (optional) posterior variance ──
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
  
  # 1. metrics + slopes --------------------------------------------------------
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
  
  # 2. posterior variance (on demand) -----------------------------------------
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


# 9.2 Prepare long‑format table (Model + Observed)
prepare_stats_long <- function(
    posterior_all,
    observed_df,
    human_dt,
    group_vars = c("sample", "format"),
    id_var = NULL) {
  
  id_var <- id_var %||%
    (c("subject_s", "subject") %>% intersect(names(observed_df)) %>% first())
  if (is.na(id_var)) stop("ID column not found.", call. = FALSE)
  
  format_tbl <- dplyr::distinct(human_dt, .data[[id_var]], format)
  
  indiv_obs <- observed_df |>
    dplyr::left_join(format_tbl, by = id_var) |>
    dplyr::select(-dplyr::all_of(id_var))
  
  dplyr::bind_rows(
    dplyr::mutate(indiv_obs,      type = "Observed"),
    dplyr::mutate(posterior_all,  type = "Model")
  ) |>
    tidyr::pivot_longer(-c(type, dplyr::all_of(group_vars)),
                        names_to = "stat", values_to = "value") |>
    dplyr::select(-dplyr::all_of(group_vars))
}


# 9.3 Histogram‑overlap statistic
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

# -----------------------------------------------------------------------------
# 10. Visualisation helpers for individual‑level diagnostics
# -----------------------------------------------------------------------------

plot_hist_overlap <- function(overlap_tbl) {
  ggplot2::ggplot(
    overlap_tbl,
    ggplot2::aes(x = reorder(stat, overlap), y = overlap)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = stat, y = 0, yend = overlap),
      colour = "grey60"
    ) +
    ggplot2::geom_point(size = 2, colour = "firebrick") +
    ggplot2::coord_flip() +
    ggplot2::labs(y = "Histogram overlap", x = NULL) +
    ggplot2::theme_bw(base_size = 7)
}

plot_density_stats <- function(
    stats_long,
    dens_ncol   = 4) {
  
  ggplot2::ggplot(
    stats_long,
    ggplot2::aes(value, fill = type, colour = type)
  ) +
    ggplot2::geom_density(alpha = .3, adjust = 1.4) +
    ggplot2::geom_rug(
      data   = subset(stats_long, grepl("intercept_", stat) & abs(value) > 30),
      sides  = "b", colour = "grey50", alpha = .5
    ) +
    ggplot2::facet_wrap(~ factor(stat), scales = "free", ncol = dens_ncol) +
    ggplot2::scale_fill_manual(values = c(Model = "#1f78b4", Observed = "#fb9a99")) +
    ggplot2::scale_colour_manual(values = c(Model = "#1f78b4", Observed = "#fb9a99")) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::labs(
      x       = NULL,
      y       = "Density",
      caption = "Intercept panels trimmed to |x| < 30; dashed rug marks extreme values"
    )
}



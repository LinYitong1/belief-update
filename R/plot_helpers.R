# =============================================================
# Helper Functions for Belief Updating Analysis
# Author: Yitong Lin
# Date: June 2025
# =============================================================

# ==================== Task Utilities ==========================
# Function: Create a unique task identifier string from BR, HR, FAR values
create_task_id <- function(BR, HR, FAR) {
  paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR)
}

# Function: Filter rows that match a specific task (optional utility)
filter_task <- function(df, BR, HR, FAR) {
  df %>% filter(BR == BR, HR == HR, FAR == FAR)
}

# ==================== Heuristic Models ========================
# Function: Add predictions from various heuristic models and Bayes to a dataframe
compute_heuristics <- function(df) {
  df %>% mutate(
    Bayes = true_posterior,      # Normative Bayesian posterior
    REP   = HR,                  # Representative heuristic (uses hit rate)
    BO    = BR,                  # Base-rate-only heuristic
    FC    = 1 - FAR,             # False-conclusion heuristic (1 - false alarm rate)
    JO    = BR * HR,             # Joint occurrence heuristic (base rate × hit rate)
    LS    = HR - FAR,            # Likelihood subtraction heuristic
    `50%` = 0.5,                 # Fixed anchor at 0.5 (pure conservatism)
    across(c(REP, BO, FC, JO, LS, `50%`, Bayes), ~ .x * 100),  # Convert to %
    response_pct = response * 100  # Convert human response to %
  )
}

# Function: Extract reference vertical lines for heuristics from a single task row
extract_reference_lines <- function(data, palette = okabe_ito) {
  data %>%
    slice(1) %>%
    dplyr::select(all_of(names(palette))) %>%
    pivot_longer(everything(), names_to = "heuristic", values_to = "value") %>%
    filter(!is.na(value), between(value, 0, 100)) %>%
    mutate(heuristic = factor(heuristic, levels = names(palette)))
}

# Function: Same as above, but takes a task label directly
extract_task_vlines <- function(df, task_label, heuristics) {
  df %>%
    filter(task_id == task_label) %>%
    slice(1) %>%
    dplyr::select(all_of(heuristics)) %>%
    pivot_longer(everything(), names_to = "heuristic", values_to = "value") %>%
    mutate(heuristic = factor(heuristic, levels = heuristics))
}
# Function:Filter out responses that match any heuristic
filter_nonmatches <- function(df, tol = 0.03) {
  preds <- data.frame(
    Bayes = round(df$true_posterior, 2),
    REP   = df$HR,
    BO    = df$BR,
    FC    = 1 - df$FAR,
    JO    = round(df$BR * df$HR, 2),
    LS    = round(df$HR - df$FAR, 2),
    `50%` = 0.5
  )
  diffs <- abs(preds - df$response)
  df[rowSums(diffs <= tol) == 0, ]
}
# ==================== Data Loading ============================
# Function: Load and clean data from different sources
load_clean_data <- function(path, type = c("experiment", "stengard", "sirota")) {
  type <- match.arg(type)
  
  if (type == "experiment") {
    # Load your own experimental data
    readxl::read_excel(path) %>%
      transmute(
        subject, format, n_trial, rt,
        BR             = br,
        HR             = hr,
        FAR            = far,
        response       = round(inputvalue / 100, 2),
        true_posterior = round(correctAnswer / 100, 2)
      )
    
  } else if (type == "stengard") {
    # Load Stengård et al. data
    readr::read_csv(path) %>%
      transmute(
        subject, format, trial,
        true_posterior, response,
        BR             = br,
        HR             = hr,
        FAR            = far,
        subject_s      = paste0(subject, "_", ifelse(format == "frequency", 1, 2))
      )
    
  } else if (type == "sirota") {
    # Load Sirota et al. data
    readr::read_csv(path) %>%
      mutate(Man = factor(Man, levels = c("Probability", "Frequency"))) %>%
      transmute(
        format       = Man,
        true_posterior,
        response,
        BR, HR, FAR,
        subject_s    = Id
      )
  } else {
    stop("Unknown type")
  }
}

# ==================== Model Simulation ========================
# Function: Simulate predictions using Bayesian Sampler models
simulate_bs_models <- function(tt, n_iter = 91, mean_v = 0.6, mean_N = 5) {
  parameter_list <- vector("list", n_iter)
  
  for (i in seq_len(n_iter)) {
    # Sample parameters from exponential and geometric distributions
    v1 <- rexp(1, rate = 1 / mean_v)
    N1 <- rgeom(1, prob = 1 / mean_N) + 1
    v2 <- rexp(1, rate = 1 / mean_v)
    N2 <- rgeom(1, prob = 1 / mean_N) + 1
    u  <- runif(1)  # asymmetry noise for BS_R
    
    # Simulate relative frequency responses
    rf1 <- simulate_and_mutate(tt, N1)$relative_frequency
    rf2 <- simulate_and_mutate(tt, N2)$relative_frequency
    
    # Store all parameters and simulated outcomes
    parameter_list[[i]] <- data.frame(
      Iteration = i,
      task_id = create_task_id(tt$BR, tt$HR, tt$FAR),
      BR = tt$BR, HR = tt$HR, FAR = tt$FAR,
      true_posterior = round(tt$true_posterior, 3),
      N1, v1, N2, v2, u,
      relative_frequency1 = rf1,
      relative_frequency2 = rf2
    )
  }
  data.table::rbindlist(parameter_list)
}

# Function: Apply BS and BS_R prediction models
predict_bs_models <- function(param_df) {
  param_df %>%
    rowwise() %>%
    mutate(
      BS_NH   = BS(N = N1, v = v1, relative_frequency = relative_frequency1),
      BS_R_NH = BS_R(N = N2, v = v2, relative_frequency = relative_frequency2, u = u)
    ) %>%
    ungroup() %>%
    pivot_longer(cols = c(BS_NH, BS_R_NH), names_to = "model", values_to = "predict") %>%
    mutate(
      predict_pct = as.numeric(predict) * 100,
      model = recode(model,
                     BS_NH   = "Bayesian Sampler",
                     BS_R_NH = "Asymmetric Bayesian Sampler"),
      model_type = "Bayesian Sampler"
    )
}

# ==================== Plotting Utilities =======================

# Okabe–Ito colorblind-friendly palette
okabe_ito <- c(
  BO    = "#E69F00", FC   = "#56B4E9", JO    = "#009E73",
  LS    = "#F0E442", REP  = "#0072B2", `50%` = "#CC79A7", Bayes = "#D55E00"
)

# Human-readable heuristic labels
heur_labels <- c(
  BO    = "BO", FC = "FC", JO = "JO",
  LS    = "LS", REP = "REP", `50%` = "50 %", Bayes = "Bayes"
)

# ggplot helper for custom color palette with legend
add_okabe_color <- function() {
  scale_color_manual(
    name   = "Heuristic",
    values = okabe_ito,
    labels = heur_labels
  )
}

# Function: Plot three-panel comparison of human responses vs. Bayesian Sampler
plot_task_panel <- function(task_label, include_legend = FALSE) {
  human    <- S_clean %>% filter(task_id == task_label)
  vlines   <- extract_task_vlines(human, task_label, heuristics)
  model_df <- pv_all %>% filter(task_id == task_label)
  
  # Panel 1: Probability format
  p_hum_prob <- ggplot(
    human %>% filter(format == "probability"),
    aes(x = response_pct)
  ) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 3, fill = "#5b5e6e") +
    geom_vline(
      data      = vlines,
      aes(xintercept = value, colour = heuristic),
      linetype  = "dashed", size = 0.6,
      show.legend = include_legend
    ) +
    labs(
      title = paste0(task_label, "\nFormat: Probability"),
      x = "Estimated probability (%)", y = "Density"
    ) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.13)) +
    add_okabe_color() +
    theme_bw(base_size = 7) +
    theme(
      legend.position = if(include_legend) "bottom" else "none",
      panel.grid.major.y = element_line(size = 0.3, colour = "grey85"),
      panel.grid.minor   = element_blank()
    )
  
  # Panel 2: Frequency format
  p_hum_freq <- ggplot(
    human %>% filter(format == "frequency"),
    aes(x = response_pct)
  ) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 3, fill = "#5b5e6e") +
    geom_vline(
      data      = vlines,
      aes(xintercept = value, colour = heuristic),
      linetype  = "dashed", size = 0.6
    ) +
    labs(
      title = "Format: Frequency",
      x = "Estimated probability (%)", y = "Density"
    ) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.13)) +
    add_okabe_color() +
    theme_bw(base_size = 7) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(size = 0.3, colour = "grey85"),
      panel.grid.minor   = element_blank()
    )
  
  # Panel 3: Model predictions
  p_mod_bs <- ggplot(
    model_df %>% filter(model_type == "Bayesian Sampler"),
    aes(x = predict_pct, fill = model)
  ) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 2,
                   alpha = 0.4, colour = NA, position = "identity") +
    labs(title = "Bayesian Sampler", x = "Predicted probability (%)", y = "Density") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.13)) +
    scale_fill_manual(
      values = c(
        "Bayesian Sampler"          = "#E41A1C",
        "Asymmetric Bayesian Sampler" = "#377EB8"
      )
    ) +
    theme_bw(base_size = 7) +
    theme(
      legend.position = if(include_legend) "bottom" else "none",
      panel.grid.major.y = element_line(size = 0.3, colour = "grey85"),
      panel.grid.minor   = element_blank()
    )
  
  # Combine panels vertically
  (p_hum_prob / p_hum_freq / p_mod_bs) +
    patchwork::plot_layout(ncol = 1, heights = c(3, 3, 3))
}

# Function: Plot participant mean response vs. true posterior
plot_mean_vs_tp <- function(df, group_vars, color_var, legend_title = NULL, file_name = NULL) {
  summary_df <- df %>%
    group_by(across(all_of(c(group_vars, "true_posterior")))) %>%
    summarise(
      mean_r = mean(response),
      sd_r   = sd(response),
      n      = n(),
      se     = sd_r / sqrt(n()),
      .groups = 'drop'
    )
  
  p <- ggplot(summary_df,
              aes(x = true_posterior, y = mean_r,
                  color = .data[[color_var]],
                  group = .data[[color_var]])) +
    geom_point(size = 1.2) +
    geom_errorbar(aes(ymin = mean_r - se, ymax = mean_r + se), width = 0.03) +
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(values = cb_palette) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "True Posterior", y = "Participant Response", color = legend_title) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  
  if (!is.null(file_name)) {
    ggsave(file.path(plot_dir, file_name), plot = p,
           width = plot_width, height = plot_height, dpi = default_dpi)
  }
  invisible(p)
}

# Function: Compute paired profiles (e.g., for probability vs. frequency format)
compute_paired_profiles <- function(df, value_col, group_col = "format") {
  grouped <- aggregate(
    as.formula(paste(value_col, "~ BR + HR + FAR +", group_col)),
    data = df, FUN = mean
  )
  prob_val <- grouped[[value_col]][grouped[[group_col]] == "probability"]
  freq_val <- grouped[[value_col]][grouped[[group_col]] == "frequency"]
  paired(prob_val, freq_val)
}

# Function: Plot paired profiles (mean or deviation) across two conditions
plot_paired_profile <- function(paired_obj, title, ylab, labels = c("Probability", "Frequency")) {
  plot(paired_obj, type = "profile") +
    theme_classic(base_size = 10) +
    labs(x = NULL, y = ylab) +
    ggtitle(title) +
    scale_x_discrete(labels = labels) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 9)
    )
}

# ==================== Constants ================================
# List of heuristics (used for plotting and extraction)
heuristics <- c("BO", "FC", "JO", "LS", "REP", "50%", "Bayes")

# Color palette for formats
cb_palette <- c("#0072B2", "#D55E00")  # Probability, Frequency

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
    readr::read_csv(path) %>%
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
      title = paste0(task_label),
      x = "Human Estimates (%)-Probability", y = "Density"
    )+
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.13)) +
    add_okabe_color() +
    theme_bw(base_size = 5) +
    theme(
      legend.position = if(include_legend) "bottom" else "none",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
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
      x = "Human Estimates (%)-Frequency", y = "Density"
    ) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.13)) +
    add_okabe_color() +
    theme_bw(base_size = 5) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  # Panel 3: Model predictions
  p_mod_bs <- ggplot(
    model_df %>% filter(model_type == "Bayesian Sampler"),
    aes(x = predict_pct, fill = model)
  ) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 3,
                   alpha = 0.4, colour = NA, position = "identity") +
    labs(x = " Model Prediction (%)", y = "Density") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.13)) +
    scale_fill_manual(
      values = c(
        "Bayesian Sampler"          = "#E41A1C",
        "Asymmetric Bayesian Sampler" = "#377EB8"
      )
    ) +
    theme_bw(base_size = 5) +
    theme(
      legend.position = if(include_legend) "bottom" else "none",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
  
  # Combine panels vertically
  (p_hum_prob / p_hum_freq / p_mod_bs) +
    patchwork::plot_layout(ncol = 1, heights = c(3, 3, 3))
}

# Function: Plot participant mean response vs. true posterior
plot_mean_vs_tp <- function(df, group_vars, color_var,
                            legend_title = NULL,
                            file_stub = NULL) {
  
  summary_df <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, "true_posterior")))) %>%
    dplyr::summarise(
      mean_r = mean(response),
      sd_r   = sd(response),
      n      = dplyr::n(),
      se     = sd_r / sqrt(n),
      .groups = "drop"
    )
  
  p <- ggplot(summary_df,
              aes(x = true_posterior, y = mean_r,
                  colour = .data[[color_var]],
                  group  = .data[[color_var]])) +
    geom_point(size = 1.2) +
    geom_errorbar(aes(ymin = mean_r - se, ymax = mean_r + se), width = 0.03) +
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_colour_manual(values = cb_palette, name = legend_title) +
    scale_x_continuous("True Posterior",
                       limits = c(0, 1), breaks = seq(0, 1, .25)) +
    scale_y_continuous("Participant Response",
                       limits = c(0, 1), breaks = seq(0, 1, .25)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_bw(base_size = 12) +
    theme(
      legend.position   = "bottom",
      legend.key.height = unit(4, "mm"),
      legend.text       = element_text(size = 10)
    )
  
  if (!is.null(file_stub)) {
    dir.create(plot_dir_tiff, showWarnings = FALSE, recursive = TRUE)
    dir.create(plot_dir_png,  showWarnings = FALSE, recursive = TRUE)

    ggsave(
      file.path(plot_dir_tiff, paste0(file_stub, ".tiff")),
      plot  = p,
      width = plot_width,
      height = plot_height,
      dpi   = 600
    )

    ggsave(
      file.path(plot_dir_png, paste0(file_stub, ".png")),
      plot  = p,
      width = plot_width,
      height = plot_height,
      dpi   = 300
    )
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

# ----- Variance -----
plot_model_accuracy <- function(df, file_stub = NULL) {
  df <- df %>%
    mutate(model_clean = model %>%
             str_replace("^H_NH$", "50%") %>% 
             str_replace("^HO_NH$", "Rep") %>% 
             str_replace("^FO_NH$", "FC") %>% 
             str_remove("_NH$") %>%
             str_replace("^BS_R$", "A_BS") %>%                     
             str_replace("^B$", "Bayes") %>%                       
             str_replace("^M_L$", "Hybrid_BS_MH") %>%              
             str_replace("^MIN_BS_", "Hybrid_BS_") %>% 
             str_replace("^Hybrid_BS_H$", "Hybrid_BS_50%")
    )
  
  heuristic_families <- c("BO", "Rep", "FC", "JO", "LS", "50%", "MH")
  
  bs_row <- df %>% filter(model_clean == "BS")
  bs_expanded <- do.call(rbind, lapply(heuristic_families, function(fam) {
    bs_row %>%
      mutate(
        model_clean  = "BS",
        family_group = fam,
        model_type   = "Bayesian Sampler"
      )
  }))
  
  other_models <- df %>%
    filter(model_clean != "BS") %>%
    mutate(
      family_group = case_when(
        model_clean %in% c("BO", "Hybrid_BS_BO") ~ "BO",
        model_clean %in% c("Rep", "Hybrid_BS_HO") ~ "Rep",
        model_clean %in% c("FC", "Hybrid_BS_FO") ~ "FC",
        model_clean %in% c("JO", "Hybrid_BS_JO") ~ "JO",
        model_clean %in% c("LS", "Hybrid_BS_LS") ~ "LS",
        model_clean %in% c("50%", "Hybrid_BS_50%") ~ "50%",
        model_clean %in% c("MH", "Hybrid_BS_MH") ~ "MH",
        model_clean == "Bayes" ~ "Bayes",
        model_clean == "LA" ~ "LA",
        model_clean == "AH" ~ "AH",
        model_clean == "LE" ~ "LE",
        model_clean == "A_BS" ~ "A_BS"
      ),
      model_type = case_when(
        grepl("^Hybrid", model_clean) ~ "Heuristic-Anchored Bayesian Sampler",
        grepl("^Bayes",  model_clean) ~ "Bayes",
        grepl("^A_BS",   model_clean) ~ "Asymmetric Bayesian Sampler",
        TRUE                           ~ "Heuristic"
      )
    )
  
  df_final <- bind_rows(other_models, bs_expanded) %>%
    filter(model_clean != "A_BS") %>%
    mutate(
      family_group = factor(family_group, levels = c(
        "BO", "Rep", "FC", "JO", "LS", "50%", "MH",
        "LA", "AH", "LE", "Bayes"
      )),
      model_type = factor(
        model_type,
        levels = c("Heuristic", "Bayesian Sampler",
                   "Heuristic-Anchored Bayesian Sampler", "Bayes")
      )
    )
  
  fill_cols <- c(
    "Heuristic"                        = cb_palette[1],
    "Bayesian Sampler"                = cb_palette[2],
    "Heuristic-Anchored Bayesian Sampler" = cb_palette[3],
    "Bayes"                           = cb_palette[4]
  )
  
  p <- ggplot(df_final, aes(x = family_group, y = Accuracy, fill = model_type)) +
    geom_bar(position = position_dodge(width = 0.7),
             stat = "identity", width = 0.6) +
    scale_fill_manual(values = fill_cols, name = NULL) +
    scale_y_continuous(
      "Proportion Accurate (< 3%)",
      limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)
    ) +
    labs(
      title = "Model Accuracy Comparison",
      x = "Model"
    ) +
    theme_bw(base_size = 9) +               
    theme(
      plot.title   = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 8),
      legend.position   = "bottom",
      legend.key.height = unit(3, "mm"),
      legend.text       = element_text(size = 8)
    )
  
  
  if (!is.null(file_stub)) {
    dir.create(plot_dir_tiff, showWarnings = FALSE, recursive = TRUE)
    dir.create(plot_dir_png,  showWarnings = FALSE, recursive = TRUE)
    
    ggsave(
      file.path(plot_dir_tiff, paste0(file_stub, ".tiff")),
      plot  = p,
      width = plot_width,
      height = plot_height,
      dpi   = 600
    )
    
    ggsave(
      file.path(plot_dir_png, paste0(file_stub, ".png")),
      plot  = p,
      width = plot_width,
      height = plot_height,
      dpi   = 300
    )
  }
  
  
  invisible(p)
}


# ==================== Constants ================================
# List of heuristics (used for plotting and extraction)
heuristics <- c("BO", "FC", "JO", "LS", "REP", "50%", "Bayes")

# Color palette for formats
cb_palette <- c("#0072B2", "#D55E00")  # Probability, Frequency

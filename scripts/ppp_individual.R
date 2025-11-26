compute_abc_by_subject <- function(human_summary,
                                   sim_params,
                                   simulations,
                                   tol = 0.05,
                                   sub = "subject") {
  human_summary <- as_tibble(human_summary)
  sim_params    <- as_tibble(sim_params)
  simulations   <- as_tibble(simulations)
  stopifnot(sub %in% names(human_summary))
  subjects <- unique(human_summary[[sub]])
  
  SUM_COLS <- setdiff(names(simulations), c("Iteration", "model"))
  
  abc_res    <- setNames(vector("list", length(subjects)), subjects)
  post_draws <- setNames(vector("list", length(subjects)), subjects)
  
  for (sbj in subjects) {
    obs_stats <- human_summary %>%
      dplyr::filter(.data[[sub]] == sbj) %>%
      dplyr::select(all_of(SUM_COLS)) %>%
      unlist(use.names = FALSE)
    
    abc_obj <- abc_posterior(obs_stats, sim_params, simulations, tol)
    draws   <- as.data.frame(abc_obj$unadj.values)
    
    draws$subject <- sbj
    
    abc_res[[as.character(sbj)]]    <- abc_obj
    post_draws[[as.character(sbj)]] <- draws
  }
  
  list(abc_res = abc_res, post_draws = post_draws, tol = tol)
}

run_abc_all_subjects <- function(human_data,
                                 sim_params,
                                 simulations,
                                 subject_col,
                                 calculate_variance = FALSE,
                                 abc_tol = 0.05) {
  human_sum <- compute_all_metrics(
    df = human_data,
    column = "response",
    group_vars = c(subject_col)
  )
  
  slope_int <- compute_SI_by(
    dt = human_data,
    group_vars  = c(subject_col),
    predictors  = c("BR","HR","FAR"),
    column      = "response"
  )
  
  final_human <- dplyr::left_join(human_sum, slope_int, by = subject_col)
  
  if (calculate_variance) {
    var_summary <- compute_variance_summary(
      human_data,
      column       = "response",
      group_vars   = c(subject_col,"BR","HR","FAR"),
      summary_vars = c(subject_col)
    )
    final_human <- dplyr::left_join(final_human, var_summary, by = subject_col)
  }
  
  abc_out <- compute_abc_by_subject(
    human_summary = final_human,
    sim_params    = sim_params,
    simulations   = simulations,
    tol           = abc_tol,
    sub           = subject_col
  )
  
  list(
    abc_out      = abc_out,
    human_summary = final_human
  )
}

abc_res_s <- run_abc_all_subjects(
  human_data        = human_dt_s,
  sim_params        = parameter_dt_s,
  simulations       = simulations_s,
  subject_col       = "subject",
  calculate_variance = FALSE,
  abc_tol           = 0.05
)
head(abc_res_s$abc_out)

plot_group_ppc_overlap <- function(
    dataset_name,
    human_data,
    abc_out,
    subject_col,
    ref_points
) {
  message("--- Running GROUP-LEVEL PPC (individual ABC) for: ", dataset_name, " ---")
  
  subjects <- unique(human_data[[subject_col]])

  ppc_list <- lapply(subjects, function(sbj) {
    target_data <- dplyr::filter(human_data, .data[[subject_col]] == sbj)
    
    draws_subject <- abc_out$post_draws[[as.character(sbj)]]
    if (is.null(draws_subject)) return(NULL)
    
    ppc_prob <- generate_ppc_global(
      draws_subject,
      dplyr::filter(target_data, format == "probability")
    ) %>%
      dplyr::mutate(
        format = "probability",
        !!subject_col := sbj
      )
    
    ppc_freq <- generate_ppc_global(
      draws_subject,
      dplyr::filter(target_data, format == "frequency")
    ) %>%
      dplyr::mutate(
        format = "frequency",
        !!subject_col := sbj
      )
    
    dplyr::bind_rows(ppc_prob, ppc_freq)
  })
  
  ppc_all <- dplyr::bind_rows(ppc_list)

  ref_lines   <- get_reference_lines(human_data, ref_points)
  ppc_prob    <- dplyr::filter(ppc_all, format == "probability")
  ppc_freq    <- dplyr::filter(ppc_all, format == "frequency")
  
  ppc_clean_p <- clean_ppc(ppc_prob, ref_lines)
  ppc_clean_f <- clean_ppc(ppc_freq, ref_lines)
  vlines_all  <- ref_lines

  human_df_prob <- human_data %>% 
    dplyr::filter(format == "probability") %>%
    dplyr::mutate(
      condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
      response_pct = response * 100
    ) %>%
    dplyr::filter(condition %in% vlines_all$condition)
  
  human_df_freq <- human_data %>% 
    dplyr::filter(format == "frequency") %>%
    dplyr::mutate(
      condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR),
      response_pct = response * 100
    ) %>%
    dplyr::filter(condition %in% vlines_all$condition)

  model_df_prob <- ppc_clean_p %>% 
    dplyr::mutate(
      response_pct = prediction * 100,
      condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR)
    ) %>%
    dplyr::filter(condition %in% vlines_all$condition)
  
  model_df_freq <- ppc_clean_f  %>% 
    dplyr::mutate(
      response_pct = prediction * 100,
      condition    = paste0("BR=", BR, ", HR=", HR, ", FAR=", FAR)
    ) %>%
    dplyr::filter(condition %in% vlines_all$condition)
  
  plot_overlay <- function(human_df, model_df, ref_lines, format_label) {
    combined_df <- dplyr::bind_rows(
      dplyr::mutate(human_df, type = "Human"),
      dplyr::mutate(model_df, type = "Model")
    )
    
    ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data = combined_df,
        ggplot2::aes(
          x = response_pct,
          y = after_stat(density),
          fill = type,
          color = type
        ),
        binwidth = 3,
        alpha    = 0.7,
        position = "identity"
      ) +
      ggplot2::geom_vline(
        data = ref_lines,
        ggplot2::aes(xintercept = value, colour = heuristic),
        linetype = "dashed",
        size     = 0.5
      ) +
      ggplot2::scale_fill_manual(values  = c("Human" = "#3c5488", "Model" = "#f39b7f")) +
      ggplot2::scale_color_manual(values = c("Human" = "#3c5488", "Model" = "#f39b7f")) +
      ggplot2::facet_wrap(~ condition, ncol = 1, scales = "free_y") +
      add_okabe_color() +
      ggplot2::coord_cartesian(ylim = c(0,1)) +
      ggplot2::labs(title = format_label, x = "Estimates (%)", y = "Density") +
      ggplot2::theme_bw(base_size = 6) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        legend.position    = "bottom",
        legend.key.height  = grid::unit(3, "mm"),
        legend.text        = ggplot2::element_text(size = 5),
        axis.title.x       = ggplot2::element_text(size = 6),
        axis.title.y       = ggplot2::element_text(size = 6),
        axis.text.x        = ggplot2::element_text(size = 5),
        axis.text.y        = ggplot2::element_text(size = 5),
        strip.text         = ggplot2::element_text(size = 6, face = "bold")
      )
  }
  
  plot_prob_overlay <- plot_overlay(
    human_df_prob, model_df_prob, ref_lines,
    "Probability Format (pooled over subjects)"
  )
  plot_freq_overlay <- plot_overlay(
    human_df_freq, model_df_freq, ref_lines,
    "Frequency Format (pooled over subjects)"
  )
  
  final_plot <- plot_prob_overlay + plot_freq_overlay +
    patchwork::plot_layout(ncol = 2, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  
  message("--- GROUP-LEVEL PPC (individual ABC) for ", dataset_name, " complete. ---")
  
  list(
    plot          = final_plot,
    ppc_all       = ppc_all,
    model_df_prob = model_df_prob,
    model_df_freq = model_df_freq
  )
}
ref_points <- tibble(
  BR  = c(0.1, 0.7, 0.9), 
  HR  = c(0.5, 0.9, 0.9), 
  FAR = c(0.3, 0.5, 0.1)
)
sub20 <- unique(human_dt_s$subject)[1:2]
human_sub20 <- dplyr::filter(human_dt_s, subject %in% sub20)
stengard_group_ppc <- plot_group_ppc_overlap(
  dataset_name = "Stengard",
  human_data   =  human_dt_s ,
  abc_out      = abc_res_s$abc_out, 
  subject_col  = "subject",
  ref_points   = ref_points
)

stengard_group_ppc$plot
stengard_group_ppc$model_df_freq
stengard_group_ppc$model_df_prob

library(ggplot2)
library(patchwork)

plot_y025 <- stengard_group_ppc$plot & 
  ggplot2::coord_cartesian(ylim = c(0, 0.15))

plot_y025

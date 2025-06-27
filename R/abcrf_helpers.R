# R/abcrf_helpers.R
# ================================================================
# Helper functions for ABC-RF model classification and posterior
# probability reconstruction from terminal node matches
# ---------------------------------------------------------------
# Author: Yitong Lin
# Created: 2025-06-27
# ================================================================

#' Compute posterior probabilities for each observed subject
#' using terminal node matches and class distributions.
#'
#' @param observed     Data frame of observed summary statistics.
#' @param simulations  Data frame of simulated statistics with factor `model`.
#' @param abc_model    A trained abcrf object (from abcrf::abcrf).
#'
#' @return A (N_obs × N_classes) matrix of posterior probabilities.
#' @export
compute_posterior <- function(observed, simulations, abc_model) {
  ## ── Dependencies ──────────────────────────────────────────────
  requireNamespace("matrixStats", quietly = TRUE)
  requireNamespace("foreach",     quietly = TRUE)
  
  ## ── Basic information ─────────────────────────────────────────
  classes    <- levels(simulations$model)
  n_classes  <- length(classes)
  
  # Terminal-node IDs for the training set and the observed set
  tn_train <- predict(abc_model$model.rf,
                      data  = simulations,
                      type  = "terminalNodes")$predictions
  tn_obs   <- predict(abc_model$model.rf,
                      data  = observed,
                      type  = "terminalNodes")$predictions
  n_trees  <- ncol(tn_train)
  
  ## ── 1.  For each tree, count how many training samples of each
  ##        class fall into every leaf ────────────────────────────
  training_record <- foreach::foreach(k = seq_len(n_trees),
                                      .packages = "matrixStats") %dopar% {
                                        leaf_train_k <- tn_train[, k]
                                        leaf_obs_k   <- tn_obs[  , k]
                                        
                                        # Largest leaf index we will need for this tree
                                        max_leaf <- max(c(0L, leaf_train_k, leaf_obs_k), na.rm = TRUE)
                                        M <- matrix(0L, nrow = max_leaf, ncol = n_classes)   # May be 0×n_classes
                                        
                                        if (max_leaf > 0L && length(leaf_train_k) > 0L) {
                                          # One-hot-encode the class factor for the training samples
                                          X <- matrix(0L, nrow = length(leaf_train_k), ncol = n_classes)
                                          X[cbind(seq_along(leaf_train_k), as.integer(simulations$model))] <- 1L
                                          
                                          # Sum rows by leaf ID
                                          M_counts <- rowsum(X, group = leaf_train_k, reorder = FALSE)
                                          idx      <- as.integer(rownames(M_counts))
                                          valid    <- idx > 0L & idx <= max_leaf
                                          M[idx[valid], ] <- M_counts[valid, , drop = FALSE]
                                        }
                                        M   # (max_leaf × n_classes) count matrix for this tree
                                      }
  
  ## ── 2.  Convert counts to probabilities for each tree, producing
  ##        an N_obs × n_classes matrix, then sum them ────────────
  posterior <- foreach::foreach(
    k        = seq_len(n_trees),
    .combine = "+",
    .packages = "matrixStats",
    .export   = c("training_record", "n_classes")
  ) %dopar% {
    
    obs_nodes <- tn_obs[, k]          # Leaf IDs for observed samples
    counts    <- training_record[[k]] # Count matrix for this tree
    
    # Start with a uniform distribution for every sample
    prob_mat <- matrix(1 / n_classes,
                       nrow = length(obs_nodes),
                       ncol = n_classes)
    
    # Only customize if the count matrix is non-empty
    if (nrow(counts) > 0L && ncol(counts) == n_classes) {
      in_bounds <- obs_nodes > 0L & obs_nodes <= nrow(counts)
      if (any(in_bounds)) {
        idx         <- obs_nodes[in_bounds]
        leaf_counts <- counts[idx, , drop = FALSE]
        row_sums    <- matrixStats::rowSums2(leaf_counts)
        
        non_empty <- row_sums > 0L
        if (any(non_empty)) {
          prob_mat[in_bounds[non_empty], ] <-
            leaf_counts[non_empty, , drop = FALSE] / row_sums[non_empty]
        }
        # Rows where row_sums == 0 keep the uniform distribution
      }
    }
    prob_mat   # N_obs × n_classes for this tree
  }
  
  ## ── 3.  Average across trees ──────────────────────────────────
  posterior <- posterior / n_trees
  colnames(posterior) <- classes
  posterior
}

#' Convert posterior matrix to long tibble for plotting or merging
#'
#' @param posterior       Matrix (N_obs × N_models) from compute_posterior()
#' @param observed        Observed data (must contain `subject` column or will assign id)
#' @param format_label    Character string for "probability" or "frequency"
#'
#' @return A tibble in long format: id, model, p, format
#' @export
to_long_posterior <- function(posterior, observed) {
  stopifnot(nrow(posterior) == nrow(observed))
  classes <- colnames(posterior)
  
  if ("subject" %in% names(observed)) {
    id_vector <- observed$subject
  } else if ("subject_s" %in% names(observed)) {
    id_vector <- observed$subject_s
  } else {
    warning("No 'subject' or 'subject_s' column found; assigning sequential IDs.")
    id_vector <- seq_len(nrow(observed))
  }
  
  tibble::as_tibble(posterior) %>%
    dplyr::mutate(id = factor(id_vector)) %>%
    tidyr::pivot_longer(cols = !id, names_to = "model", values_to = "p") %>%
    dplyr::mutate(
      model = factor(model, levels = classes),
      model_group = case_when(
        model %in% "B_NH" ~ "Bayesian",
        model %in% c("MH","HO_NH","BO_NH","FO_NH","JO_NH","LS_NH","H_NH") ~ "Single Heuristic",
        model %in% c("AH_NH","LE_NH") ~ "Adaptive Heuristic",
        model == "LA_NH" ~ "Linear Averaging",
        model %in% c("BS_NH","BS_R_NH") ~ "Bayesian Sampler",
        model == "R" ~ "Random",
        model %in% c("MIN_BS_BO","MIN_BS_HO","MIN_BS_FO",
                     "MIN_BS_JO","MIN_BS_LS","MIN_BS_H","M_L") ~ "Mixed Models",
        TRUE ~ "Other" # Default case
      )
    )
}


build_matrix <- function(df, fmt) {
  df %>%
    filter(format == fmt) %>%
    dplyr::select(id, model, p) %>%      
    pivot_wider(names_from  = model,
                values_from =p,
                values_fill = 0) %>%                  
    arrange(id) %>% 
    tibble::column_to_rownames("id") %>% 
    as.matrix()
}

# -------- Plotting Function --------
plot_bms <- function(df, title) {
  ggplot(df, aes(x = reorder(model, r), y = r, fill = group)) +
    geom_col(color = "black", width = 0.7) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.1))) +
    scale_fill_brewer(palette = "Set2") +
    labs(y = "Posterior model frequency",
         x = NULL, title = title, fill = "Model class") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 12, face = "bold")
    )
}

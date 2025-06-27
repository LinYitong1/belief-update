# models.R - Define cognitive models
#
# Description:
# Created on: 2025-06-26
# Author: Yitong Lin
# -----------------------------

# Basic Model
B <- function(BR, HR, FAR) {
  result <- BR * HR / ((BR * HR) + (1 - BR) * FAR)
  return(result)
}
BO <- function(BR) {
  result <- BR
  return(result)
}
HO <- function(HR) {
  result <- HR
  return(result)
}
FO <- function(FAR) {
  result <- (1 - FAR)
  return(result)
}
JO <- function(BR, HR) {
  result <- (BR * HR)
  return(result)
}
LS <- function(HR, FAR) {
  result <-(HR - FAR)
  return(result)
}
H <- function(...) {
  result <- 0.5
  return(result)
}
epsilon <- function() { 
  result<-runif(1, 0, 1)
  return(result)
}

create_mixed_model_1 <- function(model_name, model_fn, param_names, r_param) {
  function(...) {
    args <- list(...)
    
    # Determine num_trials from the longest data parameter
    data_params <- args[param_names]
    if (length(data_params) == 0) {
      # Handle case with no data params, e.g., H()
      num_trials <- if (!is.null(args[[r_param]])) length(args[[r_param]]) else 1
    } else {
      num_trials <- max(sapply(data_params, length))
    }
    
    # Now, handle the noise parameter 'r'
    r <- args[[r_param]]
    if (is.null(r)) {
      stop(sprintf("Missing required parameter '%s'", r_param))
    }
    if (length(r) == 1) {
      r <- rep(r, num_trials)
    } else if (length(r) != num_trials) {
      stop("'r' must be a single value or match the length of the longest parameter")
    }
    
    # Expand scalar data parameters to match num_trials
    for (p in param_names) {
      if (length(args[[p]]) == 1) {
        args[[p]] <- rep(args[[p]], num_trials)
      }
    }
    
    params <- lapply(param_names, function(p) args[[p]])
    names(params) <- param_names
    
    # The rest of the function is good...
    # ... (checking lengths, running the loop)
    if (!all(sapply(params, length) == num_trials)) {
      stop("All parameters must have the same length after expansion")
    }
    
    models <- character(num_trials)
    values <- numeric(num_trials)
    R <- runif(num_trials)
    for (i in seq_len(num_trials)) {
      if (R[i] > r[i]) {
        models[i] <- model_name
        values[i] <- do.call(model_fn, lapply(params, `[[`, i))
      } else {
        models[i] <- "Random"
        values[i] <- epsilon()
      }
    }
    
    data.frame(model = models, value = values)
  }
}


# 1-7
B_NH <- create_mixed_model_1("B", B, c("BR", "HR", "FAR"), "r")
BO_NH <- create_mixed_model_1("BO", BO, c("BR"), "r")
HO_NH <- create_mixed_model_1("HO", HO, c("HR"), "r")
FO_NH <- create_mixed_model_1("FO", FO, c("FAR"), "r")
JO_NH <- create_mixed_model_1("JO", JO, c("BR", "HR"), "r")
LS_NH <- create_mixed_model_1("LS", LS, c("HR", "FAR"), "r")
H_NH <- create_mixed_model_1("H", H, c(), "r")

# Bayesian sampler base functions
base_BS <- function(N, v, relative_frequency, prior_value) {
  result <- (N / (N + v)) * relative_frequency + 
    ((prior_value * v) / (N + v))
}

BS <- function(N, v, relative_frequency) {
  base_BS(N, v, relative_frequency, 0.5)
}

BS_R <- function(N, v, relative_frequency, u) {
  base_BS(N, v, relative_frequency, u)
}

#8-9
BS_NH <- create_mixed_model_1("BS", BS, c("N", "v", "relative_frequency"), "r")
BS_R_NH <- create_mixed_model_1("BS_R", BS_R, c("N", "v", "relative_frequency", "u"), "r")

# 10
LA <- function(BR, HR, FAR, wBR, wHR, wFAR) {
  result <- wBR * BR + wHR * HR + wFAR * FAR
  result <- (result - min(0, result)) / max(1, result)
  return(result)
}

LA_NH <- create_mixed_model_1("LA", LA, c("BR", "HR", "FAR", "wBR", "wHR", "wFAR"), "r")

# 11
select_best_heuristic_NH <-function(BR, HR, FAR, r) {
  num_trials <- length(BR)
  result <- numeric(num_trials)
  R <- runif(num_trials)
  heuristics <- vector("list", length(BR))
  differences <- vector("list", length(BR))
  for (i in seq_along(BR)) {
    if (R[i] > r) {
      d_bayes <- B(BR[i], HR[i], FAR[i])
      heuristics[[i]] <- list(
        BO(BR[i]),
        HO(HR[i]),
        FO(FAR[i]),
        JO(BR[i], HR[i]),
        LS(HR[i], FAR[i]),
        H()
      )
      differences[[i]] <- sapply(heuristics[[i]], function(d_heuristic) {
        abs(d_heuristic - d_bayes)
      })
      best_heuristic_index <- which.min(differences[[i]])
      best_heuristic_value <- heuristics[[i]][[best_heuristic_index]]
      result[i] <- best_heuristic_value
    } else {
      result[i] <- epsilon()
    }
  }
  return(result)
}

# 12
lexicographic_model_NH <- function(BR, HR, FAR, delta_HR, delta_BR, delta_FAR, r) {
  num_trials <- length(BR)
  models <- character(num_trials)
  values <- numeric(num_trials)
  R_noise <- runif(num_trials)
  
  for (i in 1:num_trials) {
    if (R_noise[i] > r) {
      in_HR_range  <- HR[i] >= (0.5 - delta_HR) && HR[i] <= (0.5 + delta_HR)
      in_BR_range  <- BR[i] >= (0.5 - delta_BR) && BR[i] <= (0.5 + delta_BR)
      in_FAR_range <- FAR[i] >= (0.5 - delta_FAR) && FAR[i] <= (0.5 + delta_FAR)
      
      if (!in_HR_range) {
        models[i] <- "HO"
        values[i] <- HR[i]
      } else if (in_HR_range && !in_BR_range) {
        models[i] <- "BO"
        values[i] <- BR[i]
      } else if (in_HR_range && in_BR_range && in_FAR_range) {
        models[i] <- "JO"
        values[i] <- BR[i] * HR[i]
      } else {
        models[i] <- "LS"
        values[i] <- HR[i] - FAR[i]
      }
    } else {
      models[i] <- "Random"
      values[i] <- epsilon()
    }
  }
  
  return(data.frame(model = models, value = values))
}

# Simulation helper
simulate_and_mutate <- function(df, num_samples) {
  df <- df %>%
    mutate(relative_frequency = numeric(n()))  
  for (i in seq_along(df$true_posterior)) {
    true_posterior <- df$true_posterior[i]
    samples <- rbinom(num_samples, 1, true_posterior)
    count_ones <- sum(samples == 1)
    df$relative_frequency[i] <- count_ones / num_samples
  }
  return(df)
}

# Integrated BS functions
In_BS_H <- function(N, v, relative_frequency) {
  base_BS(N, v, relative_frequency, 0.5)
}
In_BS_BO <- function(N, v, relative_frequency, BR) {
  base_BS(N, v, relative_frequency, BR)
}
In_BS_HO <- function(N, v, relative_frequency, HR) {
  base_BS(N, v, relative_frequency, HR)
}
In_BS_FO <- function(N, v, relative_frequency, FAR) {
  base_BS(N, v, relative_frequency, (1 - FAR))
}
In_BS_JO <- function(N, v, relative_frequency, BR, HR) {
  base_BS(N, v, relative_frequency, (BR * HR))
}
In_BS_LS <- function(N, v, relative_frequency, HR, FAR) {
  base_BS(N, v, relative_frequency, (HR - FAR))
}

# Mixed model creator
create_mixed_model_2 <- function(model_name, model_fn, param_names) {
  function(...) {
    args <- list(...)
    p_vector <- as.numeric(c(args$p_1, args$p_2, args$p_3))
    if(any(is.na(p_vector))) {
      stop("Probability values (p_1, p_2, p_3) must be numeric and non-missing.")
    }
    params <- args[param_names]
    num_trials <- if("relative_frequency" %in% names(params)) {
      length(params$relative_frequency)
    } else {
      1
    }
    
    models <- character(num_trials)
    values <- numeric(num_trials)
    
    base_info <- switch(model_name,
                        "MIN_BS_BO" = list(fn = BO, name = "BO"),
                        "MIN_BS_HO" = list(fn = HO, name = "HO"),
                        "MIN_BS_FO" = list(fn = FO, name = "FO"),
                        "MIN_BS_JO" = list(fn = JO, name = "JO"),
                        "MIN_BS_LS" = list(fn = LS, name = "LS"),
                        "MIN_BS_H"  = list(fn = H,  name = "H"),
                        stop(sprintf("Unknown model name: %s", model_name))
    )
    
    for(i in seq_len(num_trials)) {
      model_index <- rmultinom(1, size = 1, prob = p_vector)
      
      if(model_index[1] == 1) {
        models[i] <- model_name
        param_values <- lapply(params, function(x) {
          if(length(x) == 1) x else x[i]
        })
        values[i] <- do.call(model_fn, param_values)
        
      } else if(model_index[2] == 1) {
        models[i] <- base_info$name
        base_params <- switch(base_info$name,
                              "BO" = list(BR = params$BR),
                              "HO" = list(HR = params$HR),
                              "FO" = list(FAR = params$FAR),
                              "JO" = list(BR = params$BR, HR = params$HR),
                              "LS" = list(HR = params$HR, FAR = params$FAR),
                              "H" = list()
        )
        base_params_values <- lapply(base_params, function(x) {
          if(length(x) == 1) x else x[i]
        })
        values[i] <- do.call(base_info$fn, base_params_values)
        
      } else {
        models[i] <- "Random"
        values[i] <- epsilon()
      }
    }
    data.frame(model = models, value = values)
  }
}

# Create mixed BS models
MIN_BS_BO <- create_mixed_model_2("MIN_BS_BO", In_BS_BO, 
                                  param_names = c("N", "v", "relative_frequency", "BR"))
MIN_BS_HO <- create_mixed_model_2("MIN_BS_HO", In_BS_HO, 
                                  param_names = c("N", "v", "relative_frequency", "HR"))
MIN_BS_FO <- create_mixed_model_2("MIN_BS_FO", In_BS_FO, 
                                  param_names = c("N", "v", "relative_frequency", "FAR"))
MIN_BS_JO <- create_mixed_model_2("MIN_BS_JO", In_BS_JO, 
                                  param_names = c("N", "v", "relative_frequency", "BR", "HR"))
MIN_BS_LS <- create_mixed_model_2("MIN_BS_LS", In_BS_LS, 
                                  param_names = c("N", "v", "relative_frequency", "HR", "FAR"))
MIN_BS_H <- create_mixed_model_2("MIN_BS_H", In_BS_H, 
                                 param_names = c("N", "v", "relative_frequency"))


make_bs_integrated <- function(theta_fn) {
  force(theta_fn)  # capture in the closure
  function(N, v, relative_frequency, ...) {
    theta <- theta_fn(...)                             # compute θ
    base_BS(N, v, relative_frequency, theta)           # call base model
  }
}

model_specs <- list(
  BO = list(threshold_fn = function(BR) BR,               params = c("BR")),
  HO = list(threshold_fn = function(HR) HR,               params = c("HR")),
  FO = list(threshold_fn = function(FAR) (1 - FAR),         params = c("FAR")),
  JO = list(threshold_fn = function(BR, HR) (BR * HR),      params = c("BR", "HR")),
  LS = list(threshold_fn = function(HR, FAR) (HR - FAR),    params = c("HR", "FAR")),
  H  = list(threshold_fn = function() 0.5,                params = character(0))
)

# Build one integrated model per spec
integrated_model_fns <- lapply(model_specs, function(sp) {
  make_bs_integrated(sp$threshold_fn)
})

create_mixed_model_3 <- function(integrated_model_fns, model_specs) {
  force(model_specs)           
  force(integrated_model_fns)
  
  # Return callable that the user will actually run
  function(p_1, p_2, p_3, 
           q_1,q_2,q_3,q_4,q_5,q_6,                         
           N = 1L, v = 1, relative_frequency = 1,
           ...,                          
           .seed = NULL) {
    
    normalize_prob <- function(x, name){
      if(any(x < 0) || sum(x) == 0)
        stop(name, " must be non-negative and have positive sum.")
      x / sum(x)
    }
    p <- normalize_prob(c(p_1, p_2, p_3),  "p")
    q <- normalize_prob(c(q_1, q_2, q_3, q_4, q_5, q_6), "q")
    
    # ---- reproducibility ----------------------------------------------------
    if (!is.null(.seed)) set.seed(.seed)
    
    
    # ---- determine number of trials ----------------------------------------
    n_trials <- max(length(relative_frequency), 1L)
    
    
    # ---- draw categories and base model names ------------------------------
    categories <- {
      draws <- rmultinom(n_trials, 1, prob = p)      
      labs  <- c("Integrated", "Base", "Random")   
      labs[ apply(draws, 2, which.max) ]            
    }
    
    base_choices <- {
      draws <- rmultinom(n_trials, 1, prob = q)
      labs  <- names(model_specs)                    
      labs[ apply(draws, 2, which.max) ]
    }
    
    # ---- allocate output vectors -------------------------------------------
    values <- numeric(n_trials)
    labels <- character(n_trials)
    
    # ---- iterate over trials ------------------------------------------------
    for (i in seq_len(n_trials)) {
      cat_i  <- categories[i]
      base_i <- base_choices[i]
      spec   <- model_specs[[base_i]]
      
      if (cat_i == "Random") {
        # Truly random value in [0,1)
        values[i] <- epsilon()
        labels[i] <- "Random"
        
      } else {
        # Pull the extra parameters required by the threshold fn
        extra <- lapply(spec$params, function(pn) list(...)[[pn]][i])
        names(extra) <- spec$params
        
        if (cat_i == "Base") {
          # Plain threshold function
          values[i] <- do.call(spec$threshold_fn, extra)
          labels[i] <- base_i
          
        } else {  # Integrated
          values[i] <- do.call(integrated_model_fns[[base_i]],
                               c(list(N = N[i],
                                      v = v[i],
                                      relative_frequency = relative_frequency[i]),
                                 extra))
          labels[i] <- paste0("In_", base_i)
        }
      }
    }
    data.frame(model = labels, value = values, stringsAsFactors = FALSE)
  }
}


large_mixed_model <- create_mixed_model_3(integrated_model_fns, model_specs)

# 19
library(MCMCpack)

mixed_heuristic_model <- function(BR, HR, FAR, p_1, p_2, p_3, p_4, p_5, p_6, p_7) {
  # Combine the probabilities into a vector
  p_combined <- c(p_1, p_2, p_3, p_4, p_5, p_6, p_7)
  p_vector <- as.numeric(p_combined)
  num_trials <- length(BR)
  
  # Generate heuristic selections for all trials at once
  heuristic_matrix <- rmultinom(num_trials, 1, p_vector)  # 7 x 45 matrix
  
  # Identify selected heuristics
  selected_heuristics <- apply(heuristic_matrix, 2, which.max)  # Vector of length 45
  
  # Initialize vectors to store results
  models <- character(num_trials)
  values <- numeric(num_trials)
  
  # Apply heuristics based on selection
  for (i in 1:num_trials) {
    heuristic <- selected_heuristics[i]
    if (heuristic == 1) {
      models[i] <- "BO"
      values[i] <- BO(BR[i])
    } else if (heuristic == 2) {
      models[i] <- "HO"
      values[i] <- HO(HR[i])
    } else if (heuristic == 3) {
      models[i] <- "FO"
      values[i] <- FO(FAR[i])
    } else if (heuristic == 4) {
      models[i] <- "JO"
      values[i] <- JO(BR[i], HR[i])
    } else if (heuristic == 5) {
      models[i] <- "LS"
      values[i] <- LS(HR[i], FAR[i])
    } else if (heuristic == 6) {
      models[i] <- "H"
      values[i] <- H()
    } else if (heuristic == 7) {
      models[i] <- "Random"
      values[i] <- epsilon()
    } else {
      models[i] <- NA
      values[i] <- NA
    }
  }
  
  # Return the results as a data frame
  return(data.frame(model = models, value = values))
}


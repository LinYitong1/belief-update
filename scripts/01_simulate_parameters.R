# Simulate Parameters
#
# Description:
# Created on: 2025-06-26
# Author: Yitong Lin
# -----------------------------

# scripts/01_simulate_parameters.R ------------------------------------------
pacman::p_load(data.table, dplyr, tidyr, magrittr, MCMCpack, purrr,readxl,parallel)
n_cores <- detectCores() - 2
if (n_cores < 1) n_cores <- 1 
## 0. Load task design from raw experiment data ------------------------------
raw_path <- "data/df.xlsx"
S_path <- "data/s.csv"

df <- readxl::read_excel(raw_path) |>
  transmute(
    subject, format, n_trial, rt,
    BR = br,
    HR = hr,
    FAR = far,
    response = round(inputvalue / 100, 2),
    true_posterior = round(correctAnswer / 100, 2)
  )

s <- read_csv(S_path) |>
  transmute(
    subject, format,trial,
    true_posterior,response,
    BR = br,
    HR = hr,
    FAR = far,
    subject_s = paste0(subject, "_", ifelse(format == "frequency", 1, 2))
  )

# Select task parameters (BR, HR, FAR, posterior) for one subject
tt <- df %>% filter(subject == "1")
tt_s<- s %>% filter(subject_s == "1_1")
## 1. Helper functions -------------------------------------------------------

# Monte Carlo sampler: return relative frequency of success in N Bernoulli trials
simulate_and_mutate <- function(df, num_samples) {
  df <- df %>%
    mutate(relative_frequency = numeric(n()))
  for (i in seq_along(df$true_posterior)) {
    p <- df$true_posterior[i]
    df$relative_frequency[i] <- mean(rbinom(num_samples, 1, p))
  }
  return(df)
}


# Draw a Dirichlet-distributed probability vector of length k
dirichlet_weights <- function(K) as.numeric(MCMCpack::rdirichlet(1, rep(1, K)))

## 2. Parameter generation (one iteration) ----------------------------------

# Generate one complete parameter set for all trials
generate_one_iteration <- function(iter_id, tt) {
  mean_v <- 0.6 
  mean_N <- 5
  
  rate_v <- 1 / mean_v
  rate_N <- 1 / mean_N
  
  # Sample exponential parameters v1–v9 and geometric sample sizes N1–N9
  v <- rexp(9, rate = rate_v)
  N <- rgeom(9, prob = rate_N) + 1
  
  # For each trial, simulate the observed relative frequency for every trial
  rf_list <- lapply(seq_along(N), function(j) {
    simulate_and_mutate(tt, N[j])$relative_frequency
  })
  
  # Random noise (12 scalars between 0.05 and 0.95)
  r_vals <- runif(12, 0.05, 0.95)
  names(r_vals) <- paste0("r", 1:12)
  
  # Dirichlet weight vectors (used by heuristic/mixed models)
  p_list <- list(
    p1 = dirichlet_weights(3),
    p2 = dirichlet_weights(3),
    p3 = dirichlet_weights(3),
    p4 = dirichlet_weights(3),
    p5 = dirichlet_weights(3),
    p6 = dirichlet_weights(3),
    p7 = dirichlet_weights(7),
    p8 = dirichlet_weights(3),
    q  = dirichlet_weights(6)
  )
  
  # Linear weights for additive models
  wBR  <- runif(1, -3, 3)
  wHR  <- runif(1, -3, 3)
  wFAR <- runif(1, -3, 3)
  
  # Lexicographic thresholds
  delta_HR  <- runif(1, 0, 0.5)
  delta_BR  <- runif(1, 0, 0.5)
  delta_FAR <- runif(1, 0, 0.5)
  
  # Additional uniform noise (used in softmax response etc.)
  u <- runif(1)
  
  # Assemble trial-level data table with base info
  trial_dt <- data.table(
    Iteration = iter_id,
    trial_id = seq_len(nrow(tt)),
    BR = tt$BR,
    HR = tt$HR,
    FAR = tt$FAR,
    true_posterior = tt$true_posterior
  )
  
  # Broadcast all sampled scalar values to each trial
  trial_dt[, paste0("N", 1:9) := as.list(N)]
  trial_dt[, paste0("v", 1:9) := as.list(v)]
  trial_dt[, paste0("r", 1:12) := as.list(r_vals)]
  for (k in seq_along(p_list)) {
    trial_dt[, paste0(names(p_list)[k], "_", seq_along(p_list[[k]])) := as.list(p_list[[k]])]
  }
  trial_dt[, c("wBR", "wHR", "wFAR") := .(wBR, wHR, wFAR)]
  trial_dt[, c("delta_HR", "delta_BR", "delta_FAR") := .(delta_HR, delta_BR, delta_FAR)]
  trial_dt[, u := u]
  
  # Add relative frequency columns (each of length = n_trials)
  for (j in seq_along(rf_list)) {
    trial_dt[[paste0("relative_frequency", j)]] <- rf_list[[j]]
  }
  
  return(trial_dt)
}

## 3. Wrapper to run multiple iterations -------------------------------------
# Run Monte Carlo sampling for all trials and return combined parameter table
run_simulations <- function(tt, n_iter = 50, seed = 123) {
  stopifnot(all(c("BR", "HR", "FAR", "true_posterior") %in% names(tt)))
  set.seed(seed)
  
  results <- mclapply(
    seq_len(n_iter),
    function(i) generate_one_iteration(i, tt),
    mc.cores = n_cores
  )
  
  data.table::rbindlist(results)
}



## 4. Run simulation -----------------------------------------------
parameter_dt <- run_simulations(tt, n_iter = 10000)
parameter_dt_s <- run_simulations(tt_s, n_iter = 10000)
saveRDS(parameter_dt, "data/parameter_dt.rds")
saveRDS(parameter_dt_s, "data/parameter_dt_s.rds")

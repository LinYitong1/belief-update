
pacman::p_load(data.table, dplyr, tidyr, magrittr, MCMCpack, purrr, readr, parallel)
source("R/models.R")
source("R/metrics.R")
source("R/mean_variance.R")

n_cores <- detectCores() - 2
if (n_cores < 1) { n_cores <- 1 }

# Task template
raw <- read_csv("data/Exp_2.csv", show_col_types = FALSE)
if ("practice" %in% names(raw)) { raw <- raw |> filter(is.na(practice) | !practice) }
if ("is_check_trial" %in% names(raw)) { raw <- raw |> filter(is.na(is_check_trial) | !is_check_trial) }
if ("invalid_input" %in% names(raw)) { raw <- raw |> filter(!coalesce(invalid_input, FALSE)) }

complete_block_id <- raw |>
  dplyr::count(subject_id, block_num, name = "n_trials") |>
  filter(n_trials == 54) |>
  slice(1)
if (nrow(complete_block_id) == 0) { stop("No complete 54-trial block found.") }
template_subject <- complete_block_id$subject_id[1]
template_block <- complete_block_id$block_num[1]
cat("Task template subject:", template_subject, "\n")
cat("Task template block:", template_block, "\n")

tt <- raw |>
  filter(subject_id == template_subject, block_num == template_block) |>
  transmute(BR = br, HR = hr, FAR = far,
            true_posterior = round(correctAnswer / 100, 2))
stopifnot(nrow(tt) == 54)
item_check <- tt |> dplyr::count(BR, HR, FAR)
stopifnot(nrow(item_check) == 27)
stopifnot(all(item_check$n == 2))
cat("Trials:", nrow(tt), "\n")
cat("Unique problems:", nrow(item_check), "\n")

# Original helper functions
simulate_and_mutate <- function(df, num_samples) {
  df <- df |> mutate(relative_frequency = numeric(n()))
  for (i in seq_along(df$true_posterior)) {
    p <- df$true_posterior[i]
    df$relative_frequency[i] <- mean(rbinom(num_samples, 1, p))
  }
  return(df)
}

dirichlet_weights <- function(K) {
  as.numeric(MCMCpack::rdirichlet(1, rep(1, K)))
}

constrained_p_weights <- function(random_max = 0.05) {
  p_random <- runif(1, min = 0, max = random_max)
  integrated_share <- rbeta(1, shape1 = 1, shape2 = 1)
  c((1 - p_random) * integrated_share,
    (1 - p_random) * (1 - integrated_share), p_random)
}

generate_one_iteration <- function(iter_id, tt) {
  mean_v <- 0.6
  mean_N <- 5
  rate_v <- 1 / mean_v
  rate_N <- 1 / mean_N
  v9 <- rexp(1, rate = rate_v)
  N9 <- rgeom(1, prob = rate_N) + 1
  relative_frequency9 <- simulate_and_mutate(tt, N9)$relative_frequency
  p8 <- constrained_p_weights(random_max = 0.05)
  q <- dirichlet_weights(6)
  trial_dt <- data.table(
    Iteration = iter_id, trial_id = seq_len(nrow(tt)),
    BR = tt$BR, HR = tt$HR, FAR = tt$FAR,
    true_posterior = tt$true_posterior,
    N9 = N9, v9 = v9,
    p8_1 = p8[1], p8_2 = p8[2], p8_3 = p8[3],
    q_1 = q[1], q_2 = q[2], q_3 = q[3],
    q_4 = q[4], q_5 = q[5], q_6 = q[6],
    relative_frequency9 = relative_frequency9
  )
  return(trial_dt)
}

run_simulations <- function(tt, n_iter = 50, seed = 123) {
  stopifnot(all(c("BR", "HR", "FAR", "true_posterior") %in% names(tt)))
  set.seed(seed)
  results <- mclapply(
    seq_len(n_iter),
    function(i) { generate_one_iteration(i, tt) },
    mc.cores = n_cores
  )
  data.table::rbindlist(results)
}

# Simulations and M_L predictions
N_ITER <- 10000
cat("\nRunning", N_ITER, "M_L simulations...\n")
parameter_dt_ML <- run_simulations(tt, n_iter = N_ITER, seed = 123)
cat("Parameter simulation completed.\n")

generate_predictions_ML <- function(dt) {
  setDT(dt)
  pred <- dt[, .(
    M_L = large_mixed_model(
      p_1 = p8_1, p_2 = p8_2, p_3 = p8_3,
      q_1 = q_1, q_2 = q_2, q_3 = q_3,
      q_4 = q_4, q_5 = q_5, q_6 = q_6,
      N = N9, v = v9, relative_frequency = relative_frequency9,
      BR = BR, HR = HR, FAR = FAR
    )$value
  ), by = .I]
  pred[, I := NULL]
  return(cbind(dt, pred))
}

prediction_dt_ML <- generate_predictions_ML(parameter_dt_ML)
prediction_dt_ML <- prediction_dt_ML |>
  as.data.frame() |>
  mutate(model = "M_L", predict = round(as.numeric(M_L), 2),
         true_posterior = round(as.numeric(true_posterior), 2)) |>
  dplyr::select(-M_L)
prediction_dt_ML <- as.data.table(prediction_dt_ML)
cat("Prediction generation completed.\n")
cat("Trial-level simulation rows:", nrow(prediction_dt_ML), "\n")

# Summary statistics
M_Simulate_ML <- compute_all_metrics(prediction_dt_ML)
R_Simulate_ML <- compute_SI_by(prediction_dt_ML)
MeanVar_Simulate_ML <- compute_variance_summary(
  dt = prediction_dt_ML, column = "predict",
  group_vars = c("model", "Iteration", "BR", "HR", "FAR"),
  summary_vars = c("model", "Iteration")
)
V_Simulate_ML <- prediction_dt_ML[, {
  item_counts <- .SD[, .N, by = .(BR, HR, FAR)]
  complete_items <- item_counts[N == 2L, .(BR, HR, FAR)]
  n_complete <- nrow(complete_items)
  complete_dt <- merge(.SD, complete_items, by = c("BR", "HR", "FAR"), all = FALSE)
  mv_slope <- block_mv_slope(
    dt = complete_dt, item_cols = c("BR", "HR", "FAR"),
    response_col = "predict", response_scale = 1, expected_repeats = 2L
  )
  .(mv_slope = mv_slope, mv_n_items = n_complete)
}, by = .(model, Iteration)]


final_Simulate_ML <- M_Simulate_ML |>
  left_join(R_Simulate_ML, by = c("model", "Iteration")) |>
  left_join(as.data.frame(MeanVar_Simulate_ML), by = c("model", "Iteration")) |>
  left_join(as.data.frame(V_Simulate_ML), by = c("model", "Iteration"))
parameter_ML <- parameter_dt_ML |>
  as.data.frame() |>
  distinct(Iteration, N9, v9, p8_1, p8_2, p8_3,
           q_1, q_2, q_3, q_4, q_5, q_6)
saveRDS(parameter_ML, "data/Exp2_ML_Parameters_r05.rds")
ML_reference <- parameter_ML |> left_join(final_Simulate_ML, by = "Iteration")
stopifnot(nrow(ML_reference) == N_ITER)
saveRDS(ML_reference, "data/Exp2_ML_reference_r05.rds")

# Load packages and analysis helpers
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(here)
})
source(here::here("R", "models.R"))
source(here::here("R", "metrics.R"))
source(here::here("R", "mean_variance.R"))

# PPC settings and input data
SEED <- 20261011L
N_PPC_DRAWS <- 500L
set.seed(SEED)

param_cols <- c(paste0("p8_", 1:3), paste0("q_", 1:6), "N9", "v9")
raw <- readr::read_csv(here::here("data", "Exp2", "Exp_2.csv"), show_col_types = FALSE)
post <- tibble::as_tibble(readRDS(here::here("data", "Exp2", "Exp2_ML_Posterior_Draws_r05.rds")))
reference <- as.data.frame(readRDS(here::here("data", "Exp2", "Exp2_ML_reference_r05.rds")))
human_stats <- tibble::as_tibble(readRDS(here::here("data", "Exp2", "Exp2_Human_Stats.rds")))
stopifnot(all(c("subject_block", param_cols) %in% names(post)))
stopifnot(all(c("subject_id", "block_type", "format", "br", "hr", "far", "inputvalue", "correctAnswer") %in% names(raw)))

reference_stats <- setdiff(names(reference), c("Iteration", "model", param_cols, "mv_n_items"))
stopifnot(length(reference_stats) > 0L, !anyDuplicated(reference_stats))
message("Reference-table response statistics: ", length(reference_stats))

separators <- c("__", "_", " ", "|", "-", "::")
score_sep <- function(s) sum(unique(paste(raw$subject_id, raw$block_type, sep = s)) %in% unique(as.character(post$subject_block)))
scores <- vapply(separators, score_sep, integer(1))
chosen_sep <- separators[which.max(scores)]
message("Detected subject-block separator: '", chosen_sep, "'")

human <- raw %>%
  dplyr::mutate(
    subject_id = as.character(subject_id),
    subject_block = paste(subject_id, block_type, sep = chosen_sep),
    time_pressure = as.character(block_type),
    format = as.character(format),
    BR = as.numeric(br), HR = as.numeric(hr), FAR = as.numeric(far),
    true_posterior = round(as.numeric(correctAnswer) / 100, 2),
    response = round(as.numeric(inputvalue) / 100, 2)
  ) %>%
  dplyr::filter(is.finite(response), is.finite(true_posterior),
                is.finite(BR), is.finite(HR), is.finite(FAR)) %>%
  dplyr::select(subject_id, subject_block, format, time_pressure,
                BR, HR, FAR, true_posterior, response)

keys_data <- sort(unique(human$subject_block))
keys_post <- sort(unique(as.character(post$subject_block)))
if (!setequal(keys_data, keys_post)) {
  stop("Participant-block keys do not match. Data-only: ",
       paste(head(setdiff(keys_data, keys_post), 5), collapse = ", "),
       "; posterior-only: ", paste(head(setdiff(keys_post, keys_data), 5), collapse = ", "))
}
stopifnot(!anyNA(human$subject_block))
message(length(keys_data), " matching participant-blocks")

# Compute PPC summary statistics
ppc_statistics <- function(dat, response_col = "response") {
  dt <- data.table::copy(data.table::as.data.table(dat))
  stopifnot(response_col %in% names(dt), nrow(dt) > 0L)
  dt[, `:=`(
    model = "M_L",
    Iteration = 1L,
    predict = round(as.numeric(get(response_col)), 2),
    true_posterior = round(as.numeric(true_posterior), 2)
  )]
  
  M <- compute_all_metrics(dt)
  R <- compute_SI_by(dt)
  MeanVar <- compute_variance_summary(
    dt = dt, column = "predict",
    group_vars = c("model", "Iteration", "BR", "HR", "FAR"),
    summary_vars = c("model", "Iteration")
  )
  V <- dt[, {
    item_counts <- .SD[, .N, by = .(BR, HR, FAR)]
    complete_items <- item_counts[N == 2L, .(BR, HR, FAR)]
    n_complete <- nrow(complete_items)
    complete_dt <- merge(.SD, complete_items,
                         by = c("BR", "HR", "FAR"), all = FALSE)
    mv_slope <- if (n_complete >= 2L) {
      block_mv_slope(
        dt = complete_dt, item_cols = c("BR", "HR", "FAR"),
        response_col = "predict", response_scale = 1,
        expected_repeats = 2L
      )
    } else NA_real_
    .(mv_slope = mv_slope, mv_n_items = n_complete)
  }, by = .(model, Iteration)]
  
  final_stats <- as.data.frame(M) %>%
    dplyr::left_join(as.data.frame(R), by = c("model", "Iteration")) %>%
    dplyr::left_join(as.data.frame(MeanVar), by = c("model", "Iteration")) %>%
    dplyr::left_join(as.data.frame(V), by = c("model", "Iteration"))
  stopifnot(nrow(final_stats) == 1L, !anyDuplicated(names(final_stats)))
  unexpected <- setdiff(reference_stats, names(final_stats))
  if (length(unexpected)) stop("Missing reference statistics: ", paste(unexpected, collapse = ", "))
  numeric_ok <- vapply(final_stats[reference_stats], is.numeric, logical(1))
  if (!all(numeric_ok)) stop("Non-numeric statistics: ", paste(reference_stats[!numeric_ok], collapse = ", "))
  values <- unlist(final_stats[1L, reference_stats, drop = FALSE], use.names = TRUE)
  stopifnot(identical(names(values), reference_stats))
  values
}

# Simulate responses from posterior parameters
simulate_habs_once <- function(template, pars) {
  n <- nrow(template)
  N <- as.numeric(pars$N9)
  rf <- simulate_and_mutate(template, N)$relative_frequency
  
  rep_param <- function(x) rep(as.numeric(x), n)
  preds <- mapply(
    large_mixed_model,
    p_1 = rep_param(pars$p8_1), p_2 = rep_param(pars$p8_2), p_3 = rep_param(pars$p8_3),
    q_1 = rep_param(pars$q_1), q_2 = rep_param(pars$q_2), q_3 = rep_param(pars$q_3),
    q_4 = rep_param(pars$q_4), q_5 = rep_param(pars$q_5), q_6 = rep_param(pars$q_6),
    N = rep_param(pars$N9), v = rep_param(pars$v9),
    relative_frequency = rf,
    BR = template$BR, HR = template$HR, FAR = template$FAR,
    SIMPLIFY = TRUE
  )
  
  values <- if (is.matrix(preds)) preds["value", ] else {
    vapply(preds, `[[`, numeric(1), "value")
  }
  stopifnot(length(values) == n, all(is.finite(as.numeric(values))))
  as.numeric(values)
}

# Match observed statistics to participant blocks
blocks <- split(human, human$subject_block)
post_blocks <- split(post, as.character(post$subject_block))

if (!"subject_block" %in% names(human_stats)) {
  stop("Exp2_Human_Stats.rds lacks 'subject_block'. Available columns: ",
       paste(names(human_stats), collapse = ", "))
}
human_stats <- human_stats %>%
  dplyr::mutate(subject_block = as.character(subject_block))
block_metadata <- human %>%
  dplyr::distinct(subject_block, subject_id, format, time_pressure)
stopifnot(!anyDuplicated(block_metadata$subject_block))

observed <- human_stats %>%
  dplyr::select(subject_block, dplyr::all_of(reference_stats)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(reference_stats),
                      names_to = "stat", values_to = "observed") %>%
  dplyr::left_join(block_metadata, by = "subject_block") %>%
  dplyr::select(subject_block, subject_id, format, time_pressure, stat, observed)
stopifnot(nrow(observed) == length(blocks) * length(reference_stats),
          !anyNA(observed$subject_id))
message("Loaded existing ABC human statistics: ", nrow(observed), " values")

example <- ppc_statistics(blocks[[1L]], "response")
stopifnot(identical(names(example), reference_stats))

# Generate posterior predictive replicates
ppc_rep <- purrr::imap_dfr(blocks, function(dat, key) {
  draws <- post_blocks[[key]]
  if (is.null(draws) || !nrow(draws)) stop("Missing posterior draws for ", key)
  stopifnot(nrow(draws) == N_PPC_DRAWS)
  sampled_rows <- seq_len(nrow(draws))
  message("Simulating ", key, " (", nrow(dat), " trials; ", N_PPC_DRAWS, " replicates)")
  purrr::map_dfr(seq_len(N_PPC_DRAWS), function(i) {
    y_rep <- simulate_habs_once(dat, draws[sampled_rows[i], , drop = FALSE])
    dat_rep <- dplyr::mutate(dat, predicted = y_rep)
    stats <- ppc_statistics(dat_rep, "predicted")
    tibble::tibble(subject_block = key, subject_id = dat$subject_id[1L],
                   format = dat$format[1L], time_pressure = dat$time_pressure[1L],
                   replicate = i, stat = names(stats), replicated = unname(stats))
  })
})

# Compute posterior predictive p-values
ppc_tests <- ppc_rep %>%
  dplyr::inner_join(
    observed,
    by = c("subject_block", "subject_id", "format", "time_pressure", "stat")
  ) %>%
  dplyr::group_by(subject_block, subject_id, format, time_pressure, stat) %>%
  dplyr::summarise(
    observed = dplyr::first(observed),
    n_replicates = sum(is.finite(replicated)),
    n_requested = dplyr::n(),
    median_replicated = if (any(is.finite(replicated))) {
      median(replicated[is.finite(replicated)])
    } else NA_real_,
    q025 = if (any(is.finite(replicated))) {
      unname(quantile(replicated[is.finite(replicated)], .025))
    } else NA_real_,
    q975 = if (any(is.finite(replicated))) {
      unname(quantile(replicated[is.finite(replicated)], .975))
    } else NA_real_,
    upper_tail = if (is.finite(dplyr::first(observed)) &&
                     any(is.finite(replicated))) {
      mean(replicated[is.finite(replicated)] >= dplyr::first(observed))
    } else NA_real_,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ppp = 1 - 2 * abs(upper_tail - 0.5),
    valid_check = is.finite(observed) & n_replicates == N_PPC_DRAWS,
    bad_fit = ifelse(valid_check, ppp < .05, NA)
  )

stopifnot(nrow(ppc_tests) == length(blocks) * length(reference_stats))
invalid_checks <- ppc_tests %>% dplyr::filter(!valid_check)
message("Invalid/partial PPC checks: ", nrow(invalid_checks), " (see saved CSV)")

ppc_by_condition <- ppc_tests %>%
  dplyr::group_by(format, time_pressure, stat) %>%
  dplyr::summarise(
    n = dplyr::n(), n_valid = sum(valid_check),
    median_ppp = if (any(valid_check)) median(ppp[valid_check]) else NA_real_,
    q25_ppp = if (any(valid_check)) unname(quantile(ppp[valid_check], .25)) else NA_real_,
    q75_ppp = if (any(valid_check)) unname(quantile(ppp[valid_check], .75)) else NA_real_,
    proportion_ppp_below_05 = if (any(valid_check)) mean(bad_fit[valid_check]) else NA_real_,
    .groups = "drop"
  )
ppc_overall <- ppc_tests %>%
  dplyr::summarise(n_checks = dplyr::n(), n_valid = sum(valid_check),
                   median_ppp = median(ppp[valid_check], na.rm = TRUE),
                   proportion_ppp_below_05 = mean(bad_fit[valid_check], na.rm = TRUE))

out_dir <- here::here("results", "Exp2_PPC_r05")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  list(
    overall = ppc_overall,
    by_condition = ppc_by_condition,
    by_participant_block = ppc_tests
  ),
  file.path(out_dir, "Exp2_PPC_Results.rds")
)

cat("\n===== Full two-block PPC =====\n")
print(ppc_overall)
print(ppc_by_condition, n = Inf)

# Summarise first-block PPC
first_block_keys <- raw |>
  dplyr::filter(block_num == 1) |>
  dplyr::transmute(
    subject_block = paste(as.character(subject_id), block_type, sep = chosen_sep)
  ) |>
  dplyr::distinct()

ppc_first <- ppc_tests |>
  dplyr::semi_join(first_block_keys, by = "subject_block")

ppc_first_overall <- ppc_first |>
  dplyr::summarise(
    n_checks = dplyr::n(),
    n_valid = sum(valid_check, na.rm = TRUE),
    median_ppp = median(ppp[valid_check], na.rm = TRUE),
    proportion_ppp_below_05 = mean(ppp[valid_check] < .05, na.rm = TRUE)
  )

ppc_first_by_condition <- ppc_first |>
  dplyr::group_by(format, time_pressure, stat) |>
  dplyr::summarise(
    n = dplyr::n(),
    n_valid = sum(valid_check, na.rm = TRUE),
    median_ppp = median(ppp[valid_check], na.rm = TRUE),
    proportion_ppp_below_05 = mean(ppp[valid_check] < .05, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n===== First-block PPC =====\n")
print(ppc_first_overall)
print(ppc_first_by_condition, n = Inf)

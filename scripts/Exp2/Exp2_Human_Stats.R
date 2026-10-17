#Experiment 2 Human Summary Stats
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(readr)
})
source("R/metrics.R")
source("R/mean_variance.R")

# Retain participants represented in two block types; no minimum trial count is applied.
raw <- readr::read_csv("data/Exp2/Exp_2.csv", show_col_types = FALSE)
counts <- raw |>
  dplyr::count(subject_id, block_type, name = "n_valid")

keep_subjects <- counts |>
  dplyr::group_by(subject_id) |>
  dplyr::summarise(
    n_blocks = dplyr::n_distinct(block_type),
    min_valid = min(n_valid),
    .groups = "drop"
  ) |>
  dplyr::filter(n_blocks == 2L) |>
  dplyr::pull(subject_id)

human <- raw |>
  dplyr::filter(subject_id %in% keep_subjects) |>
  dplyr::transmute(
    subject = as.character(subject_id),
    format = factor(format, levels = c("frequency", "probability")),
    time_pressure = factor(block_type, levels = c("Self-paced", "Time-limited")),
    subject_block = paste(subject_id, block_type, sep = "__"),
    BR = br, HR = hr, FAR = far,
    response = round(inputvalue / 100, 2),
    true_posterior = round(correctAnswer / 100, 2)
  )
dt <- data.table::as.data.table(human)

# Compute response statistics using the simulation metric functions.
metrics <- compute_all_metrics(
  df = dt, column = "response", group_vars = "subject_block"
)
slopes <- compute_SI_by(
  dt = dt, group_vars = "subject_block",
  predictors = c("BR", "HR", "FAR"), column = "response"
)

# Use items with exactly two valid repetitions; require at least three items for mv_slope.
variance <- dt[, {
  pair_counts <- .SD[, .N, by = .(BR, HR, FAR)]
  complete_items <- pair_counts[N == 2L, .(BR, HR, FAR)]
  paired_trials <- merge(.SD, complete_items, by = c("BR", "HR", "FAR"))
  .(
    mv_slope = if (nrow(complete_items) >= 3L) {
      block_mv_slope(
        dt = paired_trials, item_cols = c("BR", "HR", "FAR"),
        response_col = "response", response_scale = 1,
        expected_repeats = 2L
      )
    } else NA_real_,
    mv_n_items = nrow(complete_items)
  )
}, by = subject_block]
mean_variance <- compute_variance_summary(
  dt = dt, column = "response",
  group_vars = c("subject_block", "BR", "HR", "FAR"),
  summary_vars = "subject_block"
)
# Retain one summary row per participant-block with fitting metadata.
metadata <- human |>
  dplyr::distinct(subject_block, subject, format, time_pressure)

human_summary <- metadata |>
  dplyr::left_join(as.data.frame(metrics), by = "subject_block") |>
  dplyr::left_join(as.data.frame(slopes), by = "subject_block") |>
  dplyr::left_join(as.data.frame(mean_variance), by = "subject_block") |>
  dplyr::left_join(as.data.frame(variance), by = "subject_block")

stopifnot(nrow(human_summary) == 2L * length(keep_subjects))
saveRDS(human_summary, "data/Exp2/Exp2_Human_Stats.rds")

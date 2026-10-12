pacman::p_load(data.table, dplyr, tidyr, purrr, abc)

ABC_TOL <- 0.05
PARAM_COLS <- c(paste0("p8_", 1:3), paste0("q_", 1:6), "N9", "v9")

# 1. Load and align reference simulations
reference <- as.data.frame(readRDS("data/Exp2_ML_reference_r05.rds"))
SUM_COLS <- setdiff(names(reference), c("Iteration", "model", PARAM_COLS, "mv_n_items"))
sim_params <- reference |> dplyr::select(Iteration, all_of(PARAM_COLS))
simulations <- reference |> dplyr::select(Iteration, model, all_of(SUM_COLS))

human_summary <- as.data.frame(
  readRDS("data/Exp2_human_Stats.rds")
)
cat("Human blocks:", nrow(human_summary), "\n")

stopifnot(nrow(human_summary) > 0)
stopifnot(!anyDuplicated(human_summary$subject_block))

SUM_COLS_COMMON <- intersect(SUM_COLS, names(human_summary))
missing_stats <- setdiff(SUM_COLS, SUM_COLS_COMMON)
cat("Reference rows:", nrow(reference),
    "| Human blocks:", nrow(human_summary),
    "| Shared statistics:", length(SUM_COLS_COMMON), "\n")
if (length(missing_stats)) cat("Missing human statistics:", paste(missing_stats, collapse = ", "), "\n")

# 2. ABC fit: same transforms and unadjusted accepted values as the original
abc_posterior_ML <- function(obs_stats, sim_params, simulations, tol = ABC_TOL) {
  sum_cols <- setdiff(names(simulations), c("Iteration", "model"))
  params <- unique(data.table::as.data.table(sim_params), by = "Iteration")
  sims <- data.table::as.data.table(simulations)
  
  # Ensure the parameter and statistic rows refer to the same simulation.
  stopifnot(identical(as.character(params$Iteration), as.character(sims$Iteration)))
  
  prob_cols <- PARAM_COLS[1:9]
  params[, (prob_cols) := lapply(.SD, function(x) pmin(pmax(x, 1e-4), 1 - 1e-4)),
         .SDcols = prob_cols]
  
  keep <- complete.cases(sims[, ..sum_cols])
  params <- params[keep]
  sims <- sims[keep]
  
  abc::abc(
    target = obs_stats,
    param = params[, ..PARAM_COLS],
    sumstat = sims[, ..sum_cols],
    method = "rejection",
    tol = tol,
    transf = c(rep("logit", 9), "log", "log"),
    logit.bounds = cbind(rep(1e-4, 9), rep(1 - 1e-4, 9))
  )
}

# 3. Fit every participant x block (both blocks)
blocks <- as.character(unique(human_summary$subject_block))
abc_results <- setNames(vector("list", length(blocks)), blocks)
posterior_draws <- setNames(vector("list", length(blocks)), blocks)
n_dropped_stats <- integer(length(blocks))

for (i in seq_along(blocks)) {
  block_id <- blocks[i]
  human_block <- human_summary[human_summary$subject_block == block_id, , drop = FALSE]
  valid <- vapply(human_block[SUM_COLS_COMMON], function(x) is.finite(as.numeric(x)), logical(1))
  block_cols <- SUM_COLS_COMMON[valid]
  n_dropped_stats[i] <- length(SUM_COLS_COMMON) - length(block_cols)
  if (!length(block_cols)) stop("No finite statistics for block: ", block_id)
  
  obs <- unlist(human_block[block_cols], use.names = FALSE)
  sims_block <- simulations |> dplyr::select(Iteration, model, all_of(block_cols))
  fit <- abc_posterior_ML(obs, sim_params, sims_block)
  
  draws <- as.data.frame(fit$unadj.values)
  stopifnot(ncol(draws) == length(PARAM_COLS))
  names(draws) <- PARAM_COLS
  draws$subject_block <- block_id
  
  abc_results[[block_id]] <- fit
  posterior_draws[[block_id]] <- draws
  cat(sprintf("[%d/%d] %s\n", i, length(blocks), block_id))
}
cat("Blocks with dropped non-finite statistics:", sum(n_dropped_stats > 0), "\n")

metadata <- human_summary |> dplyr::select(subject_block, subject, format, time_pressure) |> distinct()
posterior_draws_all <- bind_rows(posterior_draws) |> left_join(metadata, by = "subject_block")
posterior_medians <- posterior_draws_all |>
  group_by(subject_block, subject, format, time_pressure) |>
  summarise(across(all_of(PARAM_COLS), ~ median(.x, na.rm = TRUE)), .groups = "drop")
posterior_summary_long <- posterior_draws_all |>
  pivot_longer(all_of(PARAM_COLS), names_to = "parameter", values_to = "value") |>
  group_by(subject_block, subject, format, time_pressure, parameter) |>
  summarise(mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            lower_95 = quantile(value, .025, na.rm = TRUE),
            upper_95 = quantile(value, .975, na.rm = TRUE),
            n_draws = sum(is.finite(value)), .groups = "drop")

cat("Participants fitted:", n_distinct(posterior_medians$subject),
    "| Blocks fitted:", nrow(posterior_medians),
    "| Accepted draws:", nrow(posterior_draws_all), "\n")
print(posterior_medians)

# 4. Preserve the original outputs
saveRDS(posterior_draws_all, "data/Exp2_ML_Posterior_Draws_r05.rds")

# 5. Paired tests (SP vs TL) for P(Integrated) and N, within each format
compare_parameter <- function(data, parameter) {
  wide <- data |>
    dplyr::select(subject, format, time_pressure, value = all_of(parameter)) |>
    pivot_wider(names_from = time_pressure, values_from = value) |>
    mutate(diff = `Time-limited` - `Self-paced`)
  
  descriptives <- wide |>
    group_by(format) |>
    summarise(n = sum(is.finite(`Self-paced`) & is.finite(`Time-limited`)),
              mean_SP = mean(`Self-paced`, na.rm = TRUE),
              mean_TL = mean(`Time-limited`, na.rm = TRUE),
              mean_diff = mean(diff, na.rm = TRUE),
              sd_diff = sd(diff, na.rm = TRUE),
              dz = mean_diff / sd_diff, .groups = "drop") |>
    mutate(parameter = parameter)
  
  paired_tests <- wide |>
    filter(is.finite(`Self-paced`), is.finite(`Time-limited`)) |>
    group_by(format) |>
    group_modify(~ {
      tt <- t.test(.x$`Time-limited`, .x$`Self-paced`, paired = TRUE)
      tibble(t = unname(tt$statistic), df = unname(tt$parameter),
             p = tt$p.value, ci_low = tt$conf.int[1], ci_high = tt$conf.int[2])
    }) |>
    ungroup() |>
    mutate(parameter = parameter)
  
  format_test <- t.test(diff ~ format, data = wide)
  format_comparison <- tibble(parameter = parameter,
                              t = unname(format_test$statistic),
                              df = unname(format_test$parameter),
                              p = format_test$p.value,
                              ci_low = format_test$conf.int[1],
                              ci_high = format_test$conf.int[2])
  list(descriptives = descriptives, paired_tests = paired_tests,
       format_comparison = format_comparison)
}

res_pint <- compare_parameter(posterior_medians, "p8_1")
res_N <- compare_parameter(posterior_medians, "N9")
labels <- c(p8_1 = "P(Integrated)", N9 = "N")

parameter_descriptives <- bind_rows(res_pint$descriptives, res_N$descriptives) |>
  mutate(parameter = recode(parameter, !!!labels)) |>
  dplyr::select(parameter, format, n, mean_SP, mean_TL, mean_diff, sd_diff, dz)
paired_results <- bind_rows(res_pint$paired_tests, res_N$paired_tests) |>
  mutate(parameter = recode(parameter, !!!labels))
format_results <- bind_rows(res_pint$format_comparison, res_N$format_comparison) |>
  mutate(parameter = recode(parameter, !!!labels))

cat("\nParameter descriptives:\n"); print(parameter_descriptives)
cat("\nPaired tests (BH corrected across four tests):\n"); print(paired_results)
cat("\nFormat comparisons of TL-SP changes:\n"); print(format_results)

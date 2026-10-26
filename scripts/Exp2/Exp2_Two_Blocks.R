# ---------------------------------------------------------------------------
# Experiment 2: two-block analyses
# ---------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = here::here())

pacman::p_load(
  dplyr, tidyr, readr, lme4, lmerTest,
  emmeans, here, purrr, data.table, ggplot2
)

source("R/metrics.R")
source("R/models.R")
source("R/posterior_helper.R")
source("R/plot_helpers.R")
source("R/table_helper.R")

# ── 00. Constants ───────────────────────────────────────────────────────────

PARAM_COLS <- c(paste0("p8_", 1:3), paste0("q_", 1:6), "N9", "v9")
MAIN_PARAMS <- c(
  p8_1 = "Integrated", p8_2 = "Heuristic", p8_3 = "Random", N9 = "N", v9 = "v"
)
Q_LABELS <- c(q_1 = "BO", q_2 = "HO", q_3 = "FO", q_4 = "JO", q_5 = "LS", q_6 = "50%")

# ── 01. Data preparation ────────────────────────────────────────────────────
## 01.1 Behavioural data
raw <- read_csv(here::here("data", "Exp2", "Exp_2.csv"), show_col_types = FALSE)

order_metadata <- raw |>
  filter(block_num == 1) |>
  transmute(
    subject_id = as.character(subject_id),
    block_order = if_else(block_type == "Self-paced", "SP-first", "TL-first")
  ) |>
  distinct()

mf <- raw |>
  transmute(
    subject_id = as.character(subject_id),
    format = factor(format, levels = c("frequency", "probability")),
    time_pressure = factor(block_type, levels = c("Self-paced", "Time-limited")),
    block_num,
    judgment = inputvalue / 100,
    bayes = correctAnswer / 100,
    BO = br, HO = hr, FO = 1 - far, JO = br * hr, LS = hr - far
  ) |>
  left_join(order_metadata, by = "subject_id") |>
  mutate(
    block_order = factor(block_order, levels = c("SP-first", "TL-first")),
    TP = as.integer(time_pressure == "Time-limited"),
    abs_error = abs(judgment - bayes)
  )

cat("Two-block behavioural sample:", n_distinct(mf$subject_id), "participants;",
    nrow(mf), "trials\n")

## 01.2 ABC posterior parameter data
posterior_path <- here::here("data", "Exp2", "Exp2_ML_Posterior_Draws_r05.rds")
posterior_draws <- readRDS(posterior_path)

block_keys <- raw |>
  transmute(
    subject_id = as.character(subject_id),
    subject_block = paste(subject_id, block_type, sep = "__"),
    format = as.character(format),
    time_pressure = as.character(block_type),
    block_num
  ) |>
  distinct() |>
  left_join(order_metadata, by = "subject_id")

posterior_medians <- posterior_draws |>
  semi_join(block_keys, by = "subject_block") |>
  group_by(subject_block) |>
  summarise(across(all_of(PARAM_COLS), ~ median(.x, na.rm = TRUE)), .groups = "drop") |>
  left_join(block_keys, by = "subject_block") |>
  mutate(
    format = factor(format, levels = c("frequency", "probability")),
    time_pressure = factor(time_pressure, levels = c("Self-paced", "Time-limited")),
    block_order = factor(block_order, levels = c("SP-first", "TL-first")),
    TP = as.integer(time_pressure == "Time-limited")
  )

complete_subjects <- posterior_medians |>
  dplyr::count(subject_id) |>
  filter(n == 2) |>
  pull(subject_id)

posterior_medians <- posterior_medians |>
  filter(subject_id %in% complete_subjects)

cat("Two-block ABC sample:", n_distinct(posterior_medians$subject_id), "participants\n")

## 01.3 Exact-match scoring (fractional credit)
match_cols <- c("match_Bayes", "match_BO", "match_HO", "match_FO",
                "match_JO", "match_LS", "match_50")

match_scores <- calculate_match_scores(
  mf |>
    transmute(
      subject = subject_id,
      format,
      true_posterior = bayes,
      HR = HO, BR = BO, FAR = 1 - FO,
      response = round(judgment, 2)
    )
)

match_trial <- bind_cols(
  mf |> dplyr::select(subject_id, format, time_pressure, block_num, block_order, TP),
  match_scores |>
    dplyr::select(-subject, -format) |>
    dplyr::rename(
      match_Bayes = Bayes, match_BO = BO, match_HO = REP,
      match_FO = FC, match_JO = JO, match_LS = LS, match_50 = X50.
    )
)

mf <- mf |>
  mutate(exact_any = rowSums(
    dplyr::select(match_trial, dplyr::all_of(match_cols))
  ) > 0)

# ── 02. Model-fitting helper ────────────────────────────────────────────────
fit_two_block <- function(formula, fallback_corr, fallback_simple, data, label) {
  model <- tryCatch(lmer(formula, data = data), error = function(e) NULL)
  if (is.null(model) || isSingular(model)) {
    message(label, ": failed or singular; trying no-correlation structure.")
    model <- tryCatch(lmer(fallback_corr, data = data), error = function(e) NULL)
  }
  if (is.null(model) || isSingular(model)) {
    message(label, ": still failed or singular; using random-intercept-only.")
    model <- tryCatch(lmer(fallback_simple, data = data), error = function(e) NULL)
  }
  if (is.null(model) || isSingular(model)) {
    message(label, ": reduced model remains singular; inspect before reporting.")
  }
  model
}

# ── 03. Accuracy: absolute error ────────────────────────────────────────────
m_error <- fit_two_block(
  abs_error ~ format * time_pressure + block_order * time_pressure + (1 + TP | subject_id),
  abs_error ~ format * time_pressure + block_order * time_pressure + (1 + TP || subject_id),
  abs_error ~ format * time_pressure + block_order * time_pressure + (1 | subject_id),
  mf, "Absolute error"
)

error_anova <- anova(m_error, type = 3)
error_emm <- emmeans(
  m_error,
  ~ time_pressure | format,
  lmer.df = "satterthwaite",
  lmerTest.limit = 10000  
)
error_contrasts <- contrast(error_emm, method = list("TL - SP" = c(-1, 1)))

cat("\n=== Absolute error: Type III tests ===\n")
print(error_anova)

## Fixed effects: format x time-pressure interaction reported in the manuscript.
cat("\n=== Absolute error: fixed effects ===\n")
print(summary(m_error)$coefficients)

cat("\n=== Absolute error: condition means and TL - SP contrasts ===\n")
print(error_emm)
print(summary(error_contrasts, infer = TRUE))

# ── 04. Conservatism: Bayesian posterior slopes ─────────────────────────────

m_conservatism <- lmer(
  judgment ~ bayes * format * time_pressure +
    block_order * time_pressure +
    (1 + bayes * TP | subject_id),
  data = mf,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

conservatism_slopes <- emmeans::emtrends(
  m_conservatism,
  ~ format * time_pressure,
  var = "bayes",
  lmer.df = "satterthwaite",
  lmerTest.limit = 10000
)

conservatism_contrasts <- contrast(
  conservatism_slopes,
  method = list("TL - SP" = c(-1, 1)),
  by = "format"
)

cat("\n=== Conservatism: fixed effects ===\n")
print(summary(m_conservatism)$coefficients)

cat("\n=== Conservatism: posterior slopes and tests against 1 ===\n")
print(conservatism_slopes)
print(summary(conservatism_slopes, null = 1, infer = c(TRUE, TRUE)))

cat("\n=== Conservatism: TL - SP contrasts ===\n")
print(summary(conservatism_contrasts, infer = TRUE))

# ── 05. Non-exact responses (preregistered) ─────────────────────────────────

## Restricted to trials matching no candidate prediction; distinguishes
## reduced sampling from a pure random-response account.
mf_nonexact <- mf |> filter(!exact_any)

m_nonexact <- lmer(
  judgment ~ bayes * format * time_pressure +
    block_order * time_pressure +
    (1 + bayes * TP | subject_id),
  data = mf_nonexact,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

nonexact_anova <- anova(m_nonexact, type = 3)
nonexact_slopes <- emtrends(
  m_nonexact, ~ format * time_pressure, var = "bayes",
  lmer.df = "satterthwaite"
)
nonexact_contrasts <- contrast(
  nonexact_slopes,
  method = list("TL - SP" = c(-1, 1)),
  by = "format"
)

cat("\n=== Non-exact responses ===\n")
print(nonexact_anova)
print(nonexact_slopes)
print(summary(nonexact_contrasts, infer = TRUE))

# ── 06. Exact Bayesian and heuristic matching ───────────────────────────────

## Participant-level proportions of responses matching the Bayesian posterior
## and, separately, any heuristic anchor.
match_props <- match_trial |>
  dplyr::select(subject_id, format, time_pressure, block_order, TP, all_of(match_cols)) |>
  pivot_longer(all_of(match_cols), names_to = "match_type", values_to = "match") |>
  mutate(match_type = dplyr::recode(
    match_type,
    match_Bayes = "Bayes", match_BO = "BO", match_HO = "HO",
    match_FO = "FO", match_JO = "JO", match_LS = "LS",
    match_50 = "50%"
  )) |>
  group_by(subject_id, format, time_pressure, block_order, TP, match_type) |>
  summarise(match_prop = mean(match), .groups = "drop")

bayes_props <- match_props |> filter(match_type == "Bayes")

m_bayes <- fit_two_block(
  match_prop ~ format * time_pressure + block_order * time_pressure + (1 + TP | subject_id),
  match_prop ~ format * time_pressure + block_order * time_pressure + (1 + TP || subject_id),
  match_prop ~ format * time_pressure + block_order * time_pressure + (1 | subject_id),
  bayes_props, "Exact matches: Bayes"
)

heuristic_props <- match_trial |>
  mutate(
    heuristic_match =
      match_BO + match_HO + match_FO +
      match_JO + match_LS + match_50
  ) |>
  group_by(subject_id, format, time_pressure, block_order, TP) |>
  summarise(match_prop = mean(heuristic_match), .groups = "drop")

m_heuristic_order <- fit_two_block(
  match_prop ~ format * time_pressure +
    block_order * time_pressure +
    (1 + TP | subject_id),
  
  match_prop ~ format * time_pressure +
    block_order * time_pressure +
    (1 + TP || subject_id),
  
  match_prop ~ format * time_pressure +
    block_order * time_pressure +
    (1 | subject_id),
  
  heuristic_props,
  "Overall heuristic matching"
)

cat("\n=== Exact matching: Type III tests ===\n")
print(anova(m_bayes, type = 3))
print(anova(m_heuristic_order, type = 3))

## Fixed effects for the Bayesian model give the time-pressure x block-order
## interaction; condition means and TL - SP contrasts follow.
for (name in c("Bayes", "Heuristic")) {
  model <- if (name == "Bayes") m_bayes else m_heuristic_order
  
  cat("\n=== Exact", name, "matching ===\n")
  if (name == "Bayes") print(summary(model)$coefficients)
  
  emm <- emmeans(
    model,
    ~ time_pressure | format,
    lmer.df = "satterthwaite",
    lmerTest.limit = 10000
  )
  print(emm)
  print(summary(contrast(emm, method = list("TL - SP" = c(-1, 1))), infer = TRUE))
}

# ── 07. Block-order effects and Figure 15 ───────────────────────────────────

## 07.1 Mean judgments
m_judgment <- fit_two_block(
  judgment ~ format * time_pressure + block_order * time_pressure + (1 + TP | subject_id),
  judgment ~ format * time_pressure + block_order * time_pressure + (1 + TP || subject_id),
  judgment ~ format * time_pressure + block_order * time_pressure + (1 | subject_id),
  mf, "Judgment"
)

cat("\n=== Mean judgments: fixed effects ===\n")
print(summary(m_judgment)$coefficients)

## 07.2 Bayesian matching by block order
bayes_order_emm <- emmeans(
  m_bayes,
  ~ time_pressure | block_order,
  lmer.df = "satterthwaite",
  lmerTest.limit = 10000
)

cat("\n=== Bayesian matching by block order ===\n")
print(bayes_order_emm)
print(summary(
  contrast(bayes_order_emm, method = list("TL - SP" = c(-1, 1))),
  infer = TRUE
))

## 07.3 Behavioural figure data
valid_subjects <- unique(as.character(posterior_medians$subject_id))

accuracy_plot <- mf |>
  filter(subject_id %in% valid_subjects) |>
  group_by(subject_id, format, time_pressure, block_num, .drop = TRUE) |>
  summarise(value = mean(abs_error, na.rm = TRUE) * 100, .groups = "drop") |>
  mutate(measure = "Absolute error")

slope_plot <- mf |>
  filter(subject_id %in% valid_subjects) |>
  group_by(subject_id, format, time_pressure, block_num, .drop = TRUE) |>
  group_modify(~ {
    dat <- .x |>
      filter(is.finite(judgment), is.finite(bayes))
    
    slope <- if (nrow(dat) >= 3 &&
                 n_distinct(dat$bayes) >= 2) {
      unname(coef(lm(judgment ~ bayes, data = dat))[2])
    } else {
      NA_real_
    }
    
    tibble(value = slope)
  }) |>
  ungroup() |>
  mutate(measure = "Posterior slope")

matching_plot <- match_trial |>
  filter(subject_id %in% valid_subjects) |>
  mutate(
    heuristic_match =
      match_BO + match_HO + match_FO +
      match_JO + match_LS + match_50
  ) |>
  group_by(subject_id, format, time_pressure, block_num, .drop = TRUE) |>
  summarise(
    `Exact Bayesian matching` = mean(match_Bayes, na.rm = TRUE),
    `Exact heuristic matching` = mean(heuristic_match, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = c(`Exact Bayesian matching`, `Exact heuristic matching`),
    names_to = "measure",
    values_to = "value"
  )

behaviour_plot_data <- bind_rows(
  accuracy_plot,
  slope_plot,
  matching_plot
) |>
  mutate(
    measure = factor(
      measure,
      levels = c(
        "Absolute error",
        "Posterior slope",
        "Exact Bayesian matching",
        "Exact heuristic matching"
      )
    ),
    format = factor(
      format,
      levels = c("frequency", "probability"),
      labels = c("Frequency", "Probability")
    ),
    time_pressure = factor(
      time_pressure,
      levels = c("Self-paced", "Time-limited")
    )
  )

behaviour_plot_summary <- behaviour_plot_data |>
  group_by(format, measure, block_num, time_pressure, .drop = TRUE) |>
  summarise(
    n = sum(is.finite(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

## 07.4 Figure 15
dir.create(here::here("fig"), showWarnings = FALSE, recursive = TRUE)
p_behaviour <- plot_exp2_blocks(behaviour_plot_summary)
save_exp2_figure(p_behaviour, "Figure_15")

# ── 08. HABS_mixed parameter estimates ──────────────────────────────────────

## 08.1 Preregistered paired comparisons across conditions
parameter_long <- posterior_medians |>
  pivot_longer(all_of(PARAM_COLS), names_to = "parameter", values_to = "value")

paired_test <- function(data, group_col = "parameter") {
  data |>
    group_by(across(all_of(c(group_col, "format")))) |>
    group_modify(~ {
      wide <- .x |> dplyr::select(subject_id, time_pressure, value) |>
        pivot_wider(names_from = time_pressure, values_from = value)
      tt <- t.test(wide$`Time-limited`, wide$`Self-paced`, paired = TRUE)
      tibble(
        n = nrow(wide),
        mean_SP = mean(wide$`Self-paced`),
        mean_TL = mean(wide$`Time-limited`),
        difference_TL_minus_SP = mean(wide$`Time-limited` - wide$`Self-paced`),
        t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
        CI_low = tt$conf.int[1], CI_high = tt$conf.int[2]
      )
    }) |>
    ungroup()
}

main_long <- parameter_long |>
  filter(parameter %in% names(MAIN_PARAMS)) |>
  mutate(parameter = dplyr::recode(parameter, !!!MAIN_PARAMS))

q_long <- parameter_long |>
  filter(parameter %in% names(Q_LABELS)) |>
  mutate(heuristic = dplyr::recode(parameter, !!!Q_LABELS))

cat("\n=== Paired parameter comparisons (TL - SP) ===\n")
print(paired_test(main_long), n = Inf)
print(paired_test(q_long, "heuristic"), n = Inf)

## 08.2 Condition means of posterior medians
PARAM_LABELS <- c(
  p8_1 = "P(Integrated)",
  p8_2 = "P(Heuristic)",
  p8_3 = "P(Random)",
  N9   = "N",
  v9   = "v",
  q_1  = "q_BO",
  q_2  = "q_HO",
  q_3  = "q_FO",
  q_4  = "q_JO",
  q_5  = "q_LS",
  q_6  = "q_50%"
)

full_parameter_table <- posterior_medians |>
  dplyr::select(format, time_pressure, dplyr::all_of(names(PARAM_LABELS))) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(names(PARAM_LABELS)),
    names_to = "parameter",
    values_to = "value"
  ) |>
  dplyr::group_by(parameter, format, time_pressure) |>
  dplyr::summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(
    parameter = factor(
      parameter,
      levels = names(PARAM_LABELS),
      labels = unname(PARAM_LABELS)
    ),
    condition = paste(
      as.character(format),
      ifelse(time_pressure == "Self-paced", "SP", "TL"),
      sep = "_"
    ),
    estimate = dplyr::case_when(
      parameter == "N" ~ sprintf("%.2f", mean),
      parameter == "P(Random)" ~ sprintf("%.4f", mean),
      TRUE ~ sprintf("%.3f", mean)
    )
  ) |>
  dplyr::select(parameter, condition, estimate) |>
  tidyr::pivot_wider(names_from = condition, values_from = estimate) |>
  dplyr::arrange(parameter)

cat("\n=== Parameter condition means ===\n")
print(full_parameter_table, n = Inf)

# ── 09. Parameter mixed-effects models ──────────────────────────────────────

## One fit per parameter, sum-to-zero contrasts so that Type III tests and
## emmeans contrasts come from the same model.
parameter_models <- lapply(c("p8_1", "p8_2", "p8_3", "N9"), function(param) {
  data <- posterior_medians |>
    mutate(
      value = .data[[param]],
      format = factor(format),
      time_pressure = factor(time_pressure),
      block_order = factor(block_order)
    )
  contrasts(data$format) <- contr.sum(2)
  contrasts(data$time_pressure) <- contr.sum(2)
  contrasts(data$block_order) <- contr.sum(2)
  lmer(
    value ~ format * time_pressure +
      block_order * time_pressure +
      (1 | subject_id),
    data = data
  )
})
names(parameter_models) <- c("p8_1", "p8_2", "p8_3", "N9")

for (param in c("p8_1", "p8_2", "p8_3", "N9")) {
  cat("\n=== Parameter:", MAIN_PARAMS[[param]], "===\n")
  model <- parameter_models[[param]]
  print(anova(model, type = 3))
  
  emm <- emmeans(model, ~ time_pressure | format, lmer.df = "satterthwaite")
  print(summary(contrast(emm, method = list("TL - SP" = c(-1, 1))), infer = TRUE))
  
  emm_order <- emmeans(model, ~ time_pressure | block_order, lmer.df = "satterthwaite")
  print(summary(contrast(emm_order, method = list("TL - SP" = c(-1, 1))), infer = TRUE))
}

# ── 10. Figure 16 ───────────────────────────────────────────────────────────

parameter_plot_data <- posterior_medians |>
  transmute(
    subject_id = as.character(subject_id),
    format,
    time_pressure,
    block_num,
    p8_1,
    p8_2,
    p8_3,
    N9
  ) |>
  pivot_longer(
    cols = c(p8_1, p8_2, p8_3, N9),
    names_to = "measure",
    values_to = "value"
  ) |>
  mutate(
    measure = recode(
      measure,
      p8_1 = "P(Integrated)",
      p8_2 = "P(Heuristic)",
      p8_3 = "P(Random)",
      N9 = "Sample size (N)"
    ),
    measure = factor(
      measure,
      levels = c(
        "P(Integrated)",
        "P(Heuristic)",
        "P(Random)",
        "Sample size (N)"
      )
    ),
    format = factor(
      format,
      levels = c("frequency", "probability"),
      labels = c("Frequency", "Probability")
    ),
    time_pressure = factor(
      time_pressure,
      levels = c("Self-paced", "Time-limited")
    )
  )

parameter_plot_summary <- parameter_plot_data |>
  group_by(format, measure, block_num, time_pressure, .drop = TRUE) |>
  summarise(
    n = sum(is.finite(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

p_parameters <- plot_exp2_blocks(parameter_plot_summary)
save_exp2_figure(p_parameters, "Figure_16")

# ── 11. Response times (exploratory) ────────────────────────────────────────

rt_data <- raw |>
  dplyr::mutate(subject_id = as.character(subject_id)) |>
  dplyr::filter(
    subject_id %in% valid_subjects,
    is.finite(analysis_rt_sec),
    analysis_rt_sec > 0
  ) |>
  dplyr::mutate(
    format = factor(format, levels = c("frequency", "probability")),
    time_pressure = factor(block_type, levels = c("Self-paced", "Time-limited")),
    block_num = factor(block_num, levels = c(1, 2)),
    log_rt = log(analysis_rt_sec)
  )

cat("\nRT sample:", dplyr::n_distinct(rt_data$subject_id), "participants;",
    nrow(rt_data), "trials\n")

rt_subject <- rt_data |>
  dplyr::group_by(subject_id, format, block_num, time_pressure) |>
  dplyr::summarise(
    mean_log_rt = mean(log_rt),
    mean_rt = mean(analysis_rt_sec),
    median_rt = median(analysis_rt_sec),
    n_trials = dplyr::n(),
    .groups = "drop"
  )

contrasts(rt_subject$format) <- contr.sum(2)
contrasts(rt_subject$block_num) <- contr.sum(2)
contrasts(rt_subject$time_pressure) <- contr.sum(2)

m_rt <- lmer(
  mean_log_rt ~ format * block_num * time_pressure +
    (1 | subject_id),
  data = rt_subject
)

rt_emm <- emmeans(
  m_rt,
  ~ block_num | format * time_pressure,
  lmer.df = "satterthwaite"
)

cat("\n=== Response times: Type III tests ===\n")
print(anova(m_rt, type = 3))

cat("\n=== Response times: first vs second block ===\n")
print(summary(pairs(rt_emm), infer = TRUE))

# ── 12. Fitted model structures ─────────────────────────────────────────────

## Records the specification that was actually estimated, i.e. after any
## fallback applied by fit_two_block(). Use this when writing up the
## random-effects structures.
report_fits <- function(models) {
  for (nm in names(models)) {
    cat("\n", nm, "\n", sep = "")
    cat("  formula:  ", paste(deparse(formula(models[[nm]])), collapse = " "), "\n", sep = "")
    cat("  singular: ", isSingular(models[[nm]]), "\n", sep = "")
  }
}

report_fits(list(
  `03 absolute error` = m_error,
  `04 conservatism` = m_conservatism,
  `05 non-exact responses` = m_nonexact,
  `06 Bayesian matching` = m_bayes,
  `06 heuristic matching` = m_heuristic_order,
  `07 mean judgments` = m_judgment,
  `11 response times` = m_rt
))

report_fits(setNames(
  parameter_models,
  paste("09 parameter:", MAIN_PARAMS[names(parameter_models)])
))

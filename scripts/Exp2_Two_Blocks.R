knitr::opts_knit$set(root.dir = here::here())
pacman::p_load(
  dplyr, tidyr, readr, lme4, lmerTest,
  emmeans, here, clubSandwich, purrr, data.table
)
source("R/metrics.R")
source("R/models.R")
source("R/posterior_helper.R")

PARAM_COLS <- c(paste0("p8_", 1:3), paste0("q_", 1:6), "N9", "v9")
MAIN_PARAMS <- c(p8_1 = "Integrated", p8_2 = "Heuristic", p8_3 = "Random", N9 = "N", v9 = "v")
Q_LABELS <- c(q_1 = "BO", q_2 = "HO", q_3 = "FO", q_4 = "JO", q_5 = "LS", q_6 = "50%")

raw <- read_csv(here::here("data", "Exp_2.csv"), show_col_types = FALSE)
order_metadata <- raw |>
  filter(block_num == 1) |>
  transmute(subject_id = as.character(subject_id), block_order = if_else(block_type == "Self-paced", "SP-first", "TL-first")) |>
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
    abs_error = abs(judgment - bayes),
    match_Bayes = round(judgment, 2) == round(bayes, 2),
    match_BO = round(judgment, 2) == round(BO, 2),
    match_HO = round(judgment, 2) == round(HO, 2),
    match_FO = round(judgment, 2) == round(FO, 2),
    match_JO = round(judgment, 2) == round(JO, 2),
    match_LS = round(judgment, 2) == round(LS, 2),
    match_50 = round(judgment, 2) == 0.5,
    exact_any = match_Bayes | match_BO | match_HO | match_FO | match_JO | match_LS | match_50
  )
cat("Two-block behavioural sample:", n_distinct(mf$subject_id), "participants;", nrow(mf), "trials\n")

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

m_judgment <- fit_two_block(
  judgment ~ format * time_pressure + block_order * time_pressure + (1 + TP | subject_id),
  judgment ~ format * time_pressure + block_order * time_pressure + (1 + TP || subject_id),
  judgment ~ format * time_pressure + block_order * time_pressure + (1 | subject_id),
  mf, "Judgment"
)
judgment_anova <- anova(m_judgment, type = 3)
judgment_emm <- emmeans(m_judgment, ~ time_pressure | format)
judgment_contrasts <- contrast(judgment_emm, method = list("TL - SP" = c(-1, 1)))

cat("\nJudgment model:\n")
print(judgment_anova)
print(judgment_contrasts)

m_error <- fit_two_block(
  abs_error ~ format * time_pressure + block_order * time_pressure + (1 + TP | subject_id),
  abs_error ~ format * time_pressure + block_order * time_pressure + (1 + TP || subject_id),
  abs_error ~ format * time_pressure + block_order * time_pressure + (1 | subject_id),
  mf, "Absolute error"
)
error_anova <- anova(m_error, type = 3)
error_emm <- emmeans(m_error, ~ time_pressure | format)
error_contrasts <- contrast(error_emm, method = list("TL - SP" = c(-1, 1)))
cat("\nAbsolute error model:\n")
print(error_anova)
print(error_contrasts)

individual_slopes <- mf |>
  group_by(subject_id, format, time_pressure) |>
  group_modify(~ {
    fit <- lm(judgment ~ bayes, data = .x)
    tibble(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]))
  }) |>
  ungroup()
individual_slope_summary <- individual_slopes |>
  group_by(format, time_pressure) |>
  summarise(n = n(), mean = mean(slope), sd = sd(slope), median = median(slope), .groups = "drop")
cat("\nIndividual posterior slopes:\n")
print(individual_slope_summary)

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

conservatism_anova <- anova(m_conservatism, type = 3)
conservatism_slopes <- emtrends(m_conservatism, ~ format * time_pressure, var = "bayes")
conservatism_contrasts <- contrast(conservatism_slopes, method = list("TL - SP" = c(-1, 1)), by = "format")
cat("\nPosterior sensitivity:\n")
print(conservatism_anova)
print(conservatism_slopes)
print(conservatism_contrasts)

conservatism_slopes_sat <- emmeans::emtrends(
  m_conservatism,
  ~ format * time_pressure,
  var = "bayes",
  lmer.df = "satterthwaite",
  lmerTest.limit = 10000
)

cat("\nPosterior slopes tested against 1:\n")
print(summary(
  conservatism_slopes_sat,
  null = 1,
  infer = c(TRUE, TRUE)
))

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
nonexact_slopes <- emtrends(m_nonexact, ~ format * time_pressure, var = "bayes")
nonexact_contrasts <- contrast(nonexact_slopes, method = list("TL - SP" = c(-1, 1)), by = "format")
cat("\nNon-exact responses:\n")
print(nonexact_anova)
print(nonexact_slopes)
print(nonexact_contrasts)

match_cols <- c("match_Bayes", "match_BO", "match_HO", "match_FO",
                "match_JO", "match_LS", "match_50")

match_trial <- mf |>
  mutate(
    n_matches = rowSums(across(all_of(match_cols))),
    across(all_of(match_cols), ~ if_else(n_matches > 0, as.numeric(.x) / n_matches, 0))
  )

exact_match_table <- match_trial |>
  group_by(format, time_pressure) |>
  summarise(across(all_of(match_cols), mean), .groups = "drop") |>
  mutate(Heuristic = match_BO + match_HO + match_FO + match_JO + match_LS + match_50) |>
  transmute(
    format, time_pressure,
    BO = match_BO, REP = match_HO, FC = match_FO,
    JO = match_JO, LS = match_LS, `50%` = match_50,
    Heuristic, Bayes = match_Bayes
  )
cat("\nFractional-credit exact-match table:\n")
print(exact_match_table |> mutate(across(where(is.numeric), ~ round(.x, 2))), n = Inf)

match_props <- match_trial |>
  dplyr::select(subject_id, format, time_pressure, block_order, TP, all_of(match_cols)) |>
  pivot_longer(all_of(match_cols), names_to = "match_type", values_to = "match") |>
  mutate(match_type = dplyr::recode(match_type,
                                    match_Bayes = "Bayes", match_BO = "BO", match_HO = "HO",
                                    match_FO = "FO", match_JO = "JO", match_LS = "LS",
                                    match_50 = "50%")) |>
  group_by(subject_id, format, time_pressure, block_order, TP, match_type) |>
  summarise(match_prop = mean(match), .groups = "drop")

match_types <- c("Bayes", "BO", "HO", "FO", "JO", "LS", "50%")

match_models <- lapply(match_types, function(mt) {
  data <- match_props |> filter(match_type == mt)
  fit_two_block(
    match_prop ~ format * time_pressure + block_order * time_pressure + (1 + TP | subject_id),
    match_prop ~ format * time_pressure + block_order * time_pressure + (1 + TP || subject_id),
    match_prop ~ format * time_pressure + block_order * time_pressure + (1 | subject_id),
    data, paste("Exact matches:", mt)
  )
})
names(match_models) <- match_types

for (mt in match_types) {
  cat("\n===== Exact matches:", mt, "=====\n")
  model <- match_models[[mt]]
  print(anova(model, type = 3))
  emm <- emmeans(model, ~ time_pressure | format)
  print(contrast(emm, method = list("TL - SP" = c(-1, 1))))
  if (isSingular(model)) {
    m_simple <- lm(match_prop ~ format * time_pressure + block_order * time_pressure,
                   data = match_props |> filter(match_type == mt))
    robust <- clubSandwich::coef_test(
      m_simple, vcov = "CR2",
      cluster = (match_props |> filter(match_type == mt))$subject_id,
      test = "Satterthwaite"
    )
    print(robust)
  }
}

heuristic_props <- match_trial |>
  mutate(
    heuristic_match =
      match_BO + match_HO + match_FO +
      match_JO + match_LS + match_50
  ) |>
  group_by(
    subject_id, format, time_pressure,
    block_order, TP
  ) |>
  summarise(
    match_prop = mean(heuristic_match),
    .groups = "drop"
  )

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

cat("\nOverall heuristic matching and block order:\n")
print(anova(m_heuristic_order, type = 3))
posterior_path <- here::here("data", "Exp2_ML_Posterior_Draws_r05.rds")
posterior_draws <- readRDS(posterior_path)
block_keys <- raw |>
  transmute(subject_id = as.character(subject_id), subject_block = paste(subject_id, block_type, sep = "__"),
            format = as.character(format), time_pressure = as.character(block_type), block_num) |>
  distinct() |>
  left_join(order_metadata, by = "subject_id")
posterior_medians <- posterior_draws |>
  semi_join(block_keys, by = "subject_block") |>
  group_by(subject_block) |>
  summarise(across(all_of(PARAM_COLS), ~ median(.x, na.rm = TRUE)), .groups = "drop") |>
  left_join(block_keys, by = "subject_block") |>
  mutate(format = factor(format, levels = c("frequency", "probability")),
         time_pressure = factor(time_pressure, levels = c("Self-paced", "Time-limited")),
         block_order = factor(block_order, levels = c("SP-first", "TL-first")),
         TP = as.integer(time_pressure == "Time-limited"))
complete_subjects <- posterior_medians |>
  dplyr::count(subject_id) |> filter(n == 2) |> pull(subject_id)
posterior_medians <- posterior_medians |> filter(subject_id %in% complete_subjects)
cat("Two-block ABC sample:", n_distinct(posterior_medians$subject_id), "participants\n")

parameter_long <- posterior_medians |>
  pivot_longer(all_of(PARAM_COLS), names_to = "parameter", values_to = "value")
describe_params <- function(data, group_col = "parameter") {
  data |>
    group_by(across(all_of(c(group_col, "format", "time_pressure")))) |>
    summarise(n = n(), mean = mean(value), sd = sd(value), median = median(value), .groups = "drop")
}
paired_test <- function(data, group_col = "parameter") {
  data |>
    group_by(across(all_of(c(group_col, "format")))) |>
    group_modify(~ {
      wide <- .x |> dplyr::select(subject_id, time_pressure, value) |>
        pivot_wider(names_from = time_pressure, values_from = value)
      tt <- t.test(wide$`Time-limited`, wide$`Self-paced`, paired = TRUE)
      tibble(n = nrow(wide), mean_SP = mean(wide$`Self-paced`), mean_TL = mean(wide$`Time-limited`),
             difference_TL_minus_SP = mean(wide$`Time-limited` - wide$`Self-paced`),
             t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
             CI_low = tt$conf.int[1], CI_high = tt$conf.int[2])
    }) |> ungroup()
}
main_long <- parameter_long |>
  filter(parameter %in% names(MAIN_PARAMS)) |>
  mutate(parameter = dplyr::recode(parameter, !!!MAIN_PARAMS))

q_long <- parameter_long |>
  filter(parameter %in% names(Q_LABELS)) |>
  mutate(heuristic = dplyr::recode(parameter, !!!Q_LABELS))

main_parameter_descriptives <- describe_params(main_long)
main_parameter_tests <- paired_test(main_long)

q_descriptives <- describe_params(q_long, "heuristic")
q_tests <- paired_test(q_long, "heuristic")

cat("\nTwo-block main parameter descriptives:\n")
print(main_parameter_descriptives, n = Inf)
cat("\nTwo-block paired parameter tests:\n")
print(main_parameter_tests, n = Inf)

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
  dplyr::summarise(
    mean = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
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
  tidyr::pivot_wider(
    names_from = condition,
    values_from = estimate
  ) |>
  dplyr::arrange(parameter)

cat("\nAll parameter means:\n")
print(full_parameter_table, n = Inf)

parameter_order_models <- lapply(names(MAIN_PARAMS), function(param) {
  data <- posterior_medians |>
    mutate(value = .data[[param]])
  model <- lmer(
    value ~ format * time_pressure + block_order * time_pressure +
      (1 | subject_id),
    data = data
  )
  list(
    model = model,
    anova = anova(model, type = 3)
  )
})

names(parameter_order_models) <- names(MAIN_PARAMS)

for (param in c("p8_1", "p8_2", "p8_3", "N9")) {
  cat("\nParameter:", MAIN_PARAMS[[param]], "\n")
  model <- parameter_order_models[[param]]$model
  emm <- emmeans(
    model,
    ~ time_pressure | format,
    lmer.df = "satterthwaite"
  )
  print(summary(
    contrast(
      emm,
      method = list("TL - SP" = c(-1, 1))
    ),
    infer = TRUE
  ))
}

parameter_models_sum <- lapply(
  names(MAIN_PARAMS),
  function(param) {
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
    model <- lmer(
      value ~ format * time_pressure +
        block_order * time_pressure +
        (1 | subject_id),
      data = data
    )
    anova(model, type = 3)
  }
)

names(parameter_models_sum) <- names(MAIN_PARAMS)

for (param in c("p8_1", "p8_2", "p8_3", "N9")) {
  cat("\nParameter:", MAIN_PARAMS[[param]], "\n")
  print(parameter_models_sum[[param]])
}

judgment_block_emm <- emmeans(
  m_judgment,
  ~ time_pressure | block_order
)

judgment_block_contrasts <- contrast(
  judgment_block_emm,
  method = list("TL - SP" = c(-1, 1))
)

cat("\nJudgment time-pressure contrasts by block order:\n")
print(summary(judgment_block_contrasts, infer = TRUE))


library(dplyr)
library(tidyr)
library(ggplot2)

valid_subjects <- unique(as.character(posterior_medians$subject_id))

order_metadata <- raw |>
  filter(
    block_num == 1,
    as.character(subject_id) %in% valid_subjects
  ) |>
  transmute(
    subject = as.character(subject_id),
    block_order = if_else(
      block_type == "Self-paced",
      "SP-first",
      "TL-first"
    )
  ) |>
  distinct()

posterior_plot_long <- posterior_medians |>
  transmute(
    subject = as.character(subject_id),
    format = as.character(format),
    time_pressure = as.character(time_pressure),
    block_num,
    block_order = as.character(block_order),
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
    block = factor(
      block_num,
      levels = c(1, 2),
      labels = c("First block", "Second block")
    )
  ) |>
  dplyr::select(
    subject, format, time_pressure,
    block_num, block_order, block,
    measure, value
  )

match_cols <- c(
  "match_Bayes",
  "match_BO",
  "match_HO",
  "match_FO",
  "match_JO",
  "match_LS",
  "match_50"
)

match_by_subject <- raw |>
  filter(as.character(subject_id) %in% valid_subjects) |>
  transmute(
    subject = as.character(subject_id),
    format = as.character(format),
    time_pressure = as.character(block_type),
    block_num,
    judgment = inputvalue / 100,
    bayes = correctAnswer / 100,
    BO = br,
    HO = hr,
    FO = 1 - far,
    JO = br * hr,
    LS = hr - far
  ) |>
  mutate(
    match_Bayes = round(judgment, 2) == round(bayes, 2),
    match_BO = round(judgment, 2) == round(BO, 2),
    match_HO = round(judgment, 2) == round(HO, 2),
    match_FO = round(judgment, 2) == round(FO, 2),
    match_JO = round(judgment, 2) == round(JO, 2),
    match_LS = round(judgment, 2) == round(LS, 2),
    match_50 = round(judgment, 2) == 0.5
  ) |>
  mutate(
    n_matches = rowSums(across(all_of(match_cols)))
  ) |>
  mutate(
    across(
      all_of(match_cols),
      ~ if_else(
        n_matches > 0,
        as.numeric(.x) / n_matches,
        0
      )
    )
  ) |>
  left_join(order_metadata, by = "subject") |>
  group_by(
    subject, format, time_pressure,
    block_num, block_order
  ) |>
  summarise(
    `Heuristic match` = mean(
      match_BO + match_HO + match_FO +
        match_JO + match_LS + match_50
    ),
    `Bayes match` = mean(match_Bayes),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = c(`Heuristic match`, `Bayes match`),
    names_to = "measure",
    values_to = "value"
  ) |>
  mutate(
    block = factor(
      block_num,
      levels = c(1, 2),
      labels = c("First block", "Second block")
    )
  ) |>
  dplyr::select(
    subject, format, time_pressure,
    block_num, block_order, block,
    measure, value
  )

measure_levels <- c(
  "P(Integrated)",
  "P(Heuristic)",
  "P(Random)",
  "Sample size (N)",
  "Heuristic match",
  "Bayes match"
)

block_plot_data <- bind_rows(
  posterior_plot_long,
  match_by_subject
) |>
  mutate(
    measure = factor(
      measure,
      levels = measure_levels
    ),
    format = factor(
      format,
      levels = c("frequency", "probability")
    ),
    time_pressure = factor(
      time_pressure,
      levels = c("Self-paced", "Time-limited")
    ),
    block_order = factor(
      block_order,
      levels = c("SP-first", "TL-first")
    ),
    block = factor(
      block,
      levels = c("First block", "Second block")
    )
  )

block_plot_summary <- block_plot_data |>
  group_by(
    format, measure, block_order,
    block, time_pressure
  ) |>
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

p_block_order <- ggplot(
  block_plot_summary,
  aes(
    x = block,
    y = mean,
    group = block_order
  )
) +
  geom_line(
    linewidth = 0.8,
    colour = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = mean - se,
      ymax = mean + se
    ),
    width = 0.08,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_point(
    aes(
      colour = time_pressure,
      shape = time_pressure
    ),
    size = 3.2
  ) +
  scale_colour_manual(
    values = c(
      "Self-paced" = "#E78AC3",
      "Time-limited" = "#4C78A8"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Self-paced" = 17,
      "Time-limited" = 16
    )
  ) +
  facet_grid(
    measure ~ format,
    scales = "free_y",
    switch = "y"
  ) +
  labs(
    x = NULL,
    y = NULL,
    colour = "Condition",
    shape = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    strip.text.y.left = element_text(angle = 0),
    
    axis.text.x = element_text(size = 10),
    
    panel.spacing.y = grid::unit(1.5, "lines"),
    panel.spacing.x = grid::unit(0.8, "lines"),
    
    panel.border = element_rect(
      colour = "grey80",
      fill = NA,
      linewidth = 0.4
    ),
    
    plot.margin = margin(10, 15, 10, 15)
  )

dir.create(
  here::here("fig"),
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  filename = here::here(
    "fig",
    "Exp2_Block_Order_r05.png"
  ),
  plot = p_block_order,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = here::here(
    "fig",
    "Exp2_Block_Order_r05.tiff"
  ),
  plot = p_block_order,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

pacman::p_load(
  dplyr, tidyr, readr, lme4, lmerTest,
  emmeans, here, ggplot2, clubSandwich
)
emm_options(lmerTest.limit = 5000)
emm_options(pbkrtest.limit = 5000)
# Parameter definitions
PARAM_COLS <- c(
  paste0("p8_", 1:3),
  paste0("q_", 1:6),
  "N9", "v9"
)

MAIN_PARAMS <- c(
  p8_1 = "Integrated",
  p8_2 = "Heuristic",
  p8_3 = "Random",
  N9 = "N",
  v9 = "v"
)

Q_LABELS <- c(
  q_1 = "BO",
  q_2 = "HO",
  q_3 = "FO",
  q_4 = "JO",
  q_5 = "LS",
  q_6 = "50%"
)

# Prepare first-block behavioural data
raw <- read_csv(
  here::here("data", "Exp2", "Exp_2.csv"),
  show_col_types = FALSE
)
mf <- raw |>
  filter(block_num == 1) |>
  transmute(
    subject_id = as.character(subject_id),
    format = factor(format, levels = c("frequency", "probability")),
    time_pressure = factor(block_type, levels = c("Self-paced", "Time-limited")),
    judgment = inputvalue / 100,
    bayes = correctAnswer / 100,
    BO = br,
    HO = hr,
    FO = 1 - far,
    JO = br * hr,
    LS = hr - far
  ) |>
  mutate(
    abs_error = abs(judgment - bayes),
    match_Bayes = round(judgment, 2) == round(bayes, 2),
    match_BO = round(judgment, 2) == round(BO, 2),
    match_HO = round(judgment, 2) == round(HO, 2),
    match_FO = round(judgment, 2) == round(FO, 2),
    match_JO = round(judgment, 2) == round(JO, 2),
    match_LS = round(judgment, 2) == round(LS, 2),
    match_50 = round(judgment, 2) == 0.5,
    exact_any =
      match_Bayes |
      match_BO |
      match_HO |
      match_FO |
      match_JO |
      match_LS |
      match_50
  )

# Judgment and absolute-error analyses
m_judgment <- lmer(
  judgment ~ format * time_pressure + (1 | subject_id),
  data = mf
)
judgment_anova <- anova(m_judgment, type = 3)
judgment_emm <- emmeans(m_judgment, ~ time_pressure | format)
judgment_contrasts <- contrast(
  judgment_emm,
  method = list("TL - SP" = c(-1, 1))
)

m_error <- lmer(
  abs_error ~ format * time_pressure + (1 | subject_id),
  data = mf
)

error_anova <- anova(m_error, type = 3)
error_emm <- emmeans(m_error, ~ time_pressure | format)
error_contrasts <- contrast(
  error_emm,
  method = list("TL - SP" = c(-1, 1))
)

# Posterior sensitivity and conservatism
individual_slopes <- mf |>
  group_by(subject_id, format, time_pressure) |>
  group_modify(
    ~ {
      fit <- lm(judgment ~ bayes, data = .x)
      tibble(
        intercept = coef(fit)[1],
        slope = coef(fit)[2]
      )
    }
  ) |>
  ungroup()

individual_slope_summary <- individual_slopes |>
  group_by(format, time_pressure) |>
  summarise(
    n = n(),
    mean = mean(slope),
    sd = sd(slope),
    median = median(slope),
    .groups = "drop"
  )

m_conservatism <- lmer(
  judgment ~ bayes * format * time_pressure +
    (1 + bayes | subject_id),
  data = mf
)

conservatism_anova <- anova(m_conservatism, type = 3)
conservatism_slopes <- emtrends(
  m_conservatism,
  ~ format * time_pressure,
  var = "bayes"
)
conservatism_contrasts <- contrast(
  conservatism_slopes,
  method = list("TL - SP" = c(-1, 1)),
  by = "format"
)

# Posterior sensitivity for non-exact responses
mf_nonexact <- mf |>
  filter(!exact_any)
m_nonexact <- lmer(
  judgment ~ bayes * format * time_pressure +
    (1 + bayes | subject_id),
  data = mf_nonexact
)

nonexact_anova <- anova(m_nonexact, type = 3)
nonexact_slopes <- emtrends(
  m_nonexact,
  ~ format * time_pressure,
  var = "bayes"
)
nonexact_contrasts <- contrast(
  nonexact_slopes,
  method = list("TL - SP" = c(-1, 1)),
  by = "format"
)

# Exact Bayesian and heuristic matching
match_cols <- c(
  "match_Bayes", "match_BO", "match_HO",
  "match_FO", "match_JO", "match_LS", "match_50"
)
match_props <- mf |>
  mutate(
    n_matches = rowSums(across(all_of(match_cols))),
    across(
      all_of(match_cols),
      ~ if_else(n_matches > 0, as.numeric(.x) / n_matches, 0)
    )
  ) |>
  dplyr::select(subject_id, format, time_pressure, all_of(match_cols)) |>
  pivot_longer(
    all_of(match_cols),
    names_to = "match_type",
    values_to = "match"
  ) |>
  mutate(
    match_type = dplyr::recode(
      match_type,
      match_Bayes = "Bayes",
      match_BO = "BO",
      match_HO = "HO",
      match_FO = "FO",
      match_JO = "JO",
      match_LS = "LS",
      match_50 = "50%"
    )
  ) |>
  group_by(subject_id, format, time_pressure, match_type) |>
  summarise(match_prop = mean(match), .groups = "drop")

match_types <- c("Bayes", "BO", "HO", "FO", "JO", "LS", "50%")

match_models <- lapply(match_types, function(mt) {
  lm(
    match_prop ~ format * time_pressure,
    data = match_props |> filter(match_type == mt)
  )
})

names(match_models) <- match_types
heuristic_props <- match_props |>
  filter(match_type != "Bayes") |>
  group_by(subject_id, format, time_pressure) |>
  summarise(
    match_prop = sum(match_prop),
    .groups = "drop"
  )

m_heuristic <- lm(
  match_prop ~ format * time_pressure,
  data = heuristic_props
)

heuristic_anova <- anova(m_heuristic)
heuristic_emm <- emmeans(
  m_heuristic,
  ~ time_pressure | format
)
heuristic_contrasts <- contrast(
  heuristic_emm,
  method = list("TL - SP" = c(-1, 1))
)

# First-block posterior parameter summaries
posterior_path <- here::here("data", "Exp2", "Exp2_ML_Posterior_Draws_r05.rds")
posterior_draws <- readRDS(posterior_path)

first_blocks <- raw |>
  dplyr::filter(block_num == 1) |>
  dplyr::transmute(
    subject = as.character(subject_id),
    time_pressure = as.character(block_type),
    subject_block = paste(subject_id, block_type, sep = "__"),
    format = as.character(format)
  ) |>
  dplyr::distinct()

first_draws <- posterior_draws |>
  dplyr::semi_join(first_blocks, by = "subject_block") |>
  dplyr::select(subject_block, dplyr::all_of(PARAM_COLS)) |>
  dplyr::left_join(first_blocks, by = "subject_block") |>
  dplyr::mutate(
    format = factor(format, levels = c("frequency", "probability")),
    time_pressure = factor(time_pressure, levels = c("Self-paced", "Time-limited"))
  )

posterior_medians <- first_draws |>
  dplyr::group_by(subject_block, subject, format, time_pressure) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(PARAM_COLS), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )

parameter_long <- posterior_medians |>
  tidyr::pivot_longer(dplyr::all_of(PARAM_COLS),
                      names_to = "parameter", values_to = "value")

# Between-group parameter comparisons
between_test <- function(data, group_col = "parameter") {
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_col, "format")))) |>
    dplyr::group_modify(~ {
      sp <- .x$value[.x$time_pressure == "Self-paced"]
      tl <- .x$value[.x$time_pressure == "Time-limited"]
      tt <- stats::t.test(tl, sp)
      tibble::tibble(
        n_SP = length(sp), n_TL = length(tl),
        mean_SP = mean(sp), mean_TL = mean(tl),
        difference_TL_minus_SP = mean(tl) - mean(sp),
        t = unname(tt$statistic), df = unname(tt$parameter),
        p = tt$p.value, CI_low = tt$conf.int[1], CI_high = tt$conf.int[2]
      )
    }) |>
    dplyr::ungroup()
}

describe_params <- function(data, group_col = "parameter") {
  data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_col, "format", "time_pressure")))) |>
    dplyr::summarise(n = dplyr::n(), mean = mean(value), sd = sd(value),
                     median = median(value), .groups = "drop")
}

main_long <- parameter_long |>
  dplyr::filter(parameter %in% names(MAIN_PARAMS)) |>
  dplyr::mutate(parameter = dplyr::recode(parameter, !!!MAIN_PARAMS))
main_parameter_descriptives <- describe_params(main_long)
main_parameter_tests <- between_test(main_long)

q_long <- parameter_long |>
  dplyr::filter(parameter %in% names(Q_LABELS)) |>
  dplyr::mutate(heuristic = dplyr::recode(parameter, !!!Q_LABELS))
q_descriptives <- describe_params(q_long, "heuristic")
q_tests <- between_test(q_long, "heuristic")

# Report behavioural and parameter results
options(width = 200, digits = 8)
cat("\n=== First-block sample ===\n")
cat("Behavioural participants:", dplyr::n_distinct(mf$subject_id), "; trials:", nrow(mf), "\n")
print(mf |> dplyr::distinct(subject_id, format, time_pressure) |> dplyr::count(format, time_pressure, name = "n_participants"), n = Inf)
cat("ABC participants:", dplyr::n_distinct(posterior_medians$subject), "\n")
print(posterior_medians |> dplyr::count(format, time_pressure, name = "n_participants"), n = Inf)
cat("\n=== Judgment accuracy ===\n")
print(error_anova)
print(summary(error_emm, infer = TRUE))
print(summary(error_contrasts, infer = TRUE))
cat("\n=== Posterior sensitivity ===\n")
print(conservatism_anova)
print(summary(conservatism_slopes, infer = TRUE))
print(summary(conservatism_contrasts, infer = TRUE))
cat("\n=== Non-exact responses ===\n")
cat("Participants:", dplyr::n_distinct(mf_nonexact$subject_id), "; trials:", nrow(mf_nonexact), "\n")
print(nonexact_anova)
print(summary(nonexact_slopes, infer = TRUE))
print(summary(nonexact_contrasts, infer = TRUE))
print(summary(nonexact_slopes, null = 0, infer = c(TRUE, TRUE)))
print(summary(nonexact_slopes, null = 1, infer = c(TRUE, TRUE)))
cat("\n=== Exact Bayesian matching ===\n")
print(anova(match_models[["Bayes"]]))
bayes_emm <- emmeans(match_models[["Bayes"]], ~ time_pressure | format)
print(summary(bayes_emm, infer = TRUE))
print(summary(contrast(bayes_emm, method = list("TL - SP" = c(-1, 1))), infer = TRUE))
cat("\n=== Overall heuristic matching ===\n")
print(heuristic_anova)
print(summary(heuristic_emm, infer = TRUE))
print(summary(heuristic_contrasts, infer = TRUE))
cat("\n=== Main parameter estimates ===\n")
print(main_parameter_descriptives, n = Inf, width = Inf)
print(main_parameter_tests, n = Inf, width = Inf)
cat("\n=== Heuristic anchor parameters ===\n")
print(q_descriptives, n = Inf, width = Inf)
print(q_tests, n = Inf, width = Inf)

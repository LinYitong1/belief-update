pacman::p_load(
  dplyr, tidyr, readr, lme4, lmerTest,
  emmeans, here
)

# Trial-level models: allow Satterthwaite d.f. beyond emmeans' 5000-row cap.
emmeans::emm_options(lmerTest.limit = 10000)

options(width = 200, digits = 8)

# ── 00. Constants ───────────────────────────────────────────────────────────

PARAM_COLS <- c(paste0("p8_", 1:3), paste0("q_", 1:6), "N9", "v9")
MAIN_PARAMS <- c(
  p8_1 = "Integrated", p8_2 = "Heuristic", p8_3 = "Random", N9 = "N", v9 = "v"
)
Q_LABELS <- c(q_1 = "BO", q_2 = "HO", q_3 = "FO", q_4 = "JO", q_5 = "LS", q_6 = "50%")

# ── 01. Data preparation ────────────────────────────────────────────────────

raw <- read_csv(here::here("data", "Exp2", "Exp_2.csv"), show_col_types = FALSE)

mf <- raw |>
  filter(block_num == 1) |>
  transmute(
    subject_id = as.character(subject_id),
    format = factor(format, levels = c("frequency", "probability")),
    time_pressure = factor(block_type, levels = c("Self-paced", "Time-limited")),
    judgment = inputvalue / 100,
    bayes = correctAnswer / 100,
    BO = br, HO = hr, FO = 1 - far, JO = br * hr, LS = hr - far
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
      match_Bayes | match_BO | match_HO | match_FO |
      match_JO | match_LS | match_50
  )

## Fractional credit: responses matching several predictions receive 1 / n.
match_cols <- c(
  "match_Bayes", "match_BO", "match_HO",
  "match_FO", "match_JO", "match_LS", "match_50"
)

match_trial <- mf |>
  mutate(
    n_matches = rowSums(across(all_of(match_cols))),
    across(all_of(match_cols), ~ if_else(n_matches > 0, as.numeric(.x) / n_matches, 0))
  )

## 01.1 Sample composition
cat("\n=== First-block sample ===\n")
cat("Behavioural participants:", n_distinct(mf$subject_id), "; trials:", nrow(mf), "\n")
print(
  mf |>
    distinct(subject_id, format, time_pressure) |>
    count(format, time_pressure, name = "n_participants"),
  n = Inf
)

# ── 02. Accuracy: absolute error ────────────────────────────────────────────

m_error <- lmer(
  abs_error ~ format * time_pressure + (1 | subject_id),
  data = mf
)

error_anova <- anova(m_error, type = 3)
error_emm <- emmeans(m_error, ~ time_pressure | format, lmer.df = "satterthwaite")
error_contrasts <- contrast(error_emm, method = list("TL - SP" = c(-1, 1)))

cat("\n=== Absolute error: Type III tests ===\n")
print(error_anova)

## Fixed effects: format x time-pressure interaction
cat("\n=== Absolute error: fixed effects ===\n")
print(summary(m_error)$coefficients)

## emmeans are on the 0-1 scale; multiply by 100 for percentages.
cat("\n=== Absolute error: condition means and TL - SP contrasts ===\n")
print(summary(error_emm, infer = TRUE))
print(summary(error_contrasts, infer = TRUE))

# ── 03. Conservatism: Bayesian posterior slopes ─────────────────────────────

m_conservatism <- lmer(
  judgment ~ bayes * format * time_pressure +
    (1 + bayes | subject_id),
  data = mf
)

conservatism_anova <- anova(m_conservatism, type = 3)
conservatism_slopes <- emtrends(
  m_conservatism,
  ~ format * time_pressure,
  var = "bayes",
  lmer.df = "satterthwaite"
)
conservatism_contrasts <- contrast(
  conservatism_slopes,
  method = list("TL - SP" = c(-1, 1)),
  by = "format"
)

cat("\n=== Conservatism: Type III tests ===\n")
print(conservatism_anova)

cat("\n=== Conservatism: posterior slopes and tests against 1 ===\n")
print(summary(conservatism_slopes, infer = TRUE))
print(summary(conservatism_slopes, null = 1, infer = c(TRUE, TRUE)))

cat("\n=== Conservatism: TL - SP contrasts ===\n")
print(summary(conservatism_contrasts, infer = TRUE))

# ── 04. Non-exact responses ─────────────────────────────────────────────────

## Restricted to trials matching no candidate prediction; distinguishes
## reduced sampling from a pure random-response account.
mf_nonexact <- mf |> filter(!exact_any)

m_nonexact <- lmer(
  judgment ~ bayes * format * time_pressure +
    (1 + bayes | subject_id),
  data = mf_nonexact
)

nonexact_anova <- anova(m_nonexact, type = 3)
nonexact_slopes <- emtrends(
  m_nonexact,
  ~ format * time_pressure,
  var = "bayes",
  lmer.df = "satterthwaite"
)
nonexact_contrasts <- contrast(
  nonexact_slopes,
  method = list("TL - SP" = c(-1, 1)),
  by = "format"
)

cat("\n=== Non-exact responses ===\n")
cat("Participants:", n_distinct(mf_nonexact$subject_id), "; trials:", nrow(mf_nonexact), "\n")
print(nonexact_anova)
print(summary(nonexact_slopes, infer = TRUE))
print(summary(nonexact_contrasts, infer = TRUE))

## Slopes above 0 rule out random responding; below 1 indicates conservatism.
print(summary(nonexact_slopes, null = 0, infer = c(TRUE, TRUE)))
print(summary(nonexact_slopes, null = 1, infer = c(TRUE, TRUE)))

# ── 05. Exact Bayesian and heuristic matching ───────────────────────────────

## Participant-level proportions; time pressure is between participants here,
## so these are fitted as linear models.
match_props <- match_trial |>
  dplyr::select(subject_id, format, time_pressure, all_of(match_cols)) |>
  pivot_longer(all_of(match_cols), names_to = "match_type", values_to = "match") |>
  mutate(match_type = recode(
    match_type,
    match_Bayes = "Bayes", match_BO = "BO", match_HO = "HO",
    match_FO = "FO", match_JO = "JO", match_LS = "LS",
    match_50 = "50%"
  )) |>
  group_by(subject_id, format, time_pressure, match_type) |>
  summarise(match_prop = mean(match), .groups = "drop")

m_bayes <- lm(
  match_prop ~ format * time_pressure,
  data = match_props |> filter(match_type == "Bayes")
)

heuristic_props <- match_props |>
  filter(match_type != "Bayes") |>
  group_by(subject_id, format, time_pressure) |>
  summarise(match_prop = sum(match_prop), .groups = "drop")

m_heuristic <- lm(match_prop ~ format * time_pressure, data = heuristic_props)

bayes_emm <- emmeans(m_bayes, ~ time_pressure | format)
heuristic_emm <- emmeans(m_heuristic, ~ time_pressure | format)

cat("\n=== Exact Bayesian matching ===\n")
print(anova(m_bayes))
print(summary(bayes_emm, infer = TRUE))
print(summary(contrast(bayes_emm, method = list("TL - SP" = c(-1, 1))), infer = TRUE))

cat("\n=== Overall heuristic matching ===\n")
print(anova(m_heuristic))
print(summary(heuristic_emm, infer = TRUE))
print(summary(contrast(heuristic_emm, method = list("TL - SP" = c(-1, 1))), infer = TRUE))

# ── 06. HABS_mixed parameter comparisons ────────────────────────────────────

posterior_draws <- readRDS(here::here("data", "Exp2", "Exp2_ML_Posterior_Draws_r05.rds"))

first_blocks <- raw |>
  filter(block_num == 1) |>
  transmute(
    subject = as.character(subject_id),
    time_pressure = as.character(block_type),
    subject_block = paste(subject_id, block_type, sep = "__"),
    format = as.character(format)
  ) |>
  distinct()

posterior_medians <- posterior_draws |>
  semi_join(first_blocks, by = "subject_block") |>
  dplyr::select(subject_block, all_of(PARAM_COLS)) |>
  left_join(first_blocks, by = "subject_block") |>
  mutate(
    format = factor(format, levels = c("frequency", "probability")),
    time_pressure = factor(time_pressure, levels = c("Self-paced", "Time-limited"))
  ) |>
  group_by(subject_block, subject, format, time_pressure) |>
  summarise(
    across(all_of(PARAM_COLS), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  )

## Independent-samples Welch t tests within each format.
between_test <- function(data, group_col = "parameter") {
  data |>
    group_by(across(all_of(c(group_col, "format")))) |>
    group_modify(~ {
      sp <- .x$value[.x$time_pressure == "Self-paced"]
      tl <- .x$value[.x$time_pressure == "Time-limited"]
      tt <- t.test(tl, sp)
      tibble(
        n_SP = length(sp), n_TL = length(tl),
        mean_SP = mean(sp), mean_TL = mean(tl),
        difference_TL_minus_SP = mean(tl) - mean(sp),
        t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
        CI_low = tt$conf.int[1], CI_high = tt$conf.int[2]
      )
    }) |>
    ungroup()
}

parameter_long <- posterior_medians |>
  pivot_longer(all_of(PARAM_COLS), names_to = "parameter", values_to = "value")

main_long <- parameter_long |>
  filter(parameter %in% names(MAIN_PARAMS)) |>
  mutate(parameter = recode(parameter, !!!MAIN_PARAMS))

q_long <- parameter_long |>
  filter(parameter %in% names(Q_LABELS)) |>
  mutate(heuristic = recode(parameter, !!!Q_LABELS))

cat("\n=== First-block ABC sample ===\n")
print(posterior_medians |> count(format, time_pressure, name = "n_participants"), n = Inf)

cat("\n=== Main parameter comparisons (TL - SP) ===\n")
print(between_test(main_long), n = Inf, width = Inf)

cat("\n=== Heuristic anchor parameters (TL - SP) ===\n")
print(between_test(q_long, "heuristic"), n = Inf, width = Inf)

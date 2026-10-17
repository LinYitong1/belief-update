# Dependencies
pacman::p_load(
  dplyr, magrittr, tidyr, stringr, data.table,
  ggplot2, ggpubr, gridExtra, cowplot, forcats,
  ranger, abcrf,
  foreach, doParallel,
  readxl,purrr,MCMCpack,
  afex, emmeans, abc,readr,parallel
)
pacman::p_load(data.table, dplyr, tidyr, magrittr, MCMCpack, purrr,readxl,parallel)

source("R/metrics.R")

# Load simulated predictions
dt <- as.data.table(readRDS("data/Exp1/prediction_dt_r05.rds"))
dt_s <- as.data.table(readRDS("data/Stengard/prediction_dt_s_r05.rds"))

# Stengard simulation summaries
A_Simulate_s <- compute_CLC_summary(dt_s)
M_Simulate_s <- compute_all_metrics(dt_s)
R_Simulate_s <- compute_SI_by(dt_s)

final_Simulate_s <- M_Simulate_s %>%
  left_join(R_Simulate_s, by = c("model", "Iteration"))

# Stengard human summaries
df_s <- read_csv("data/Stengard/s.csv") |>
  transmute(
    subject, format,trial,
    true_posterior,response,
    BR = br,
    HR = hr,
    FAR = far,
    subject_s = paste0(subject, "_", ifelse(format == "frequency", 1, 2))
  )

human_dt_s <- as.data.table(df_s)
M_Human_s <- compute_all_metrics(
  df = human_dt_s,
  column = "response",
  group_vars = c("subject_s","format")
)
R_Human_s <- compute_SI_by(
  dt = human_dt_s,
  group_vars = c("subject_s","format"),
  predictors = c("BR", "HR", "FAR"),
  column = "response"
)
final_Human_s <- M_Human_s %>%
  left_join(R_Human_s, by = c("subject_s","format"))

saveRDS(final_Human_s, "data/Stengard/Human_Summary_dts.rds")
saveRDS(final_Simulate_s, "data/Stengard/Simulate_Summary_dts_r05.rds")

# experiment 1 
A_Simulate <- compute_CLC_summary(dt)
M_Simulate <- compute_all_metrics(dt)
R_Simulate <- compute_SI_by(dt)
V_Simulate <- compute_variance_summary(dt)

final_Simulate <- M_Simulate %>%
  left_join(R_Simulate, by = c("model", "Iteration")) %>%
  left_join(V_Simulate, by = c("model", "Iteration"))

df <- read_csv("data/Exp1/df.csv") |>
  transmute(
    subject, format, n_trial, rt,
    BR = br,
    HR = hr,
    FAR = far,
    response = round(inputvalue / 100, 2),
    true_posterior = round(correctAnswer / 100, 2)
  )

human_dt <- as.data.table(df)
M_Human <- compute_all_metrics(
  df = human_dt,
  column = "response",
  group_vars = c("subject","format")
)
R_Human <- compute_SI_by(
  dt = human_dt,
  group_vars = c("subject","format"),
  predictors = c("BR", "HR", "FAR"),
  column = "response"
)
# Mean within-item response variance for each participant and format.
V_Human <- compute_variance_summary(
  dt = human_dt,
  column = "response",
  group_vars = c("subject", "BR", "HR", "FAR","format"),
  summary_vars = c("subject","format")
)

p_V_Human<-ggplot(V_Human, aes(x = mean_variance)) +
  geom_histogram(binwidth = 0.002, fill = "#4575b4", color = "black") +
  labs(
    x = "Mean Variance Across Repeated Items",
    y = "Number of Participants") +
  theme_minimal()

ggsave("fig/variance.png", p_V_Human, width = 6, height = 4.5, dpi = 300)

final_Human <- M_Human %>%
  left_join(R_Human, by = c("subject","format")) %>%
  left_join(V_Human, by = c("subject","format"))

saveRDS(final_Human, "data/Exp1/Human_Summary_dt.rds")
saveRDS(final_Simulate, "data/Exp1/Simulate_Summary_dt_r05.rds")

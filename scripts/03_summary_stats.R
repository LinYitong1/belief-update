# -----------------------------
# Analysis & Plots
# Description: Runs all summary statistics, model fitting (ABC-RF), and human model classification
# Created on: 2025-06-26
# Author: Yitong Lin
# -----------------------------

# ----- Load Required Packages -----
pacman::p_load(
  # Data manipulation
  dplyr, magrittr, tidyr, stringr, data.table,
  # Visualization
  ggplot2, ggpubr, gridExtra, cowplot, forcats,
  # Random forest
  ranger, abcrf,
  # Parallelism
  foreach, doParallel,
  # File I/O
  readxl,purrr,MCMCpack,
  # Statistics
  afex, emmeans, abc,readr,parallel
)
pacman::p_load(data.table, dplyr, tidyr, magrittr, MCMCpack, purrr,readxl,parallel)
# ----- Load Custom Metric Functions -----
source("R/metrics.R")

# ----- Load Simulated Model Prediction Data -----
dt <- as.data.table(readRDS("data/prediction_dt.rds"))
dt_s<-as.data.table(readRDS("data/prediction_dt_s.rds"))

# ===== Stengard =====
# ----- Compute Accuracy (CLC) -----
# Accuracy = % of predictions within 3% of the true posterior
A_Simulate_s <- compute_CLC_summary(dt_s)
# ----- Plot Accuracy Bar Chart -----
p1<-ggplot(A_Simulate_s, aes(x = reorder(model, Accuracy), y = Accuracy)) +
  geom_bar(stat = "identity", fill = "#4B6C8A", width = 0.6) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Model Accuracy Comparison",
    x = "Model",
    y = "Proportion Accurate (< 3%)"
  )
# ========== Figure 5: Model Comparison (accuracy) for Stengard ==========
ggsave("fig/Figure_5.png", p1, width = 10, height = 10, dpi = 600, units = "in")
# ----- Compute All Summary Statistics for Simulated Data -----
# EMP, Ad, PD metrics (by reference heuristics)
M_Simulate_s <- compute_all_metrics(dt_s)

# Linear regression (slope, intercept) by BR, HR, FAR
R_Simulate_s <- compute_SI_by(dt_s)

# ----- Merge All Metrics into One Data Table -----
final_Simulate_s <- M_Simulate_s %>%
  left_join(R_Simulate_s, by = c("model", "Iteration")) 
# ===== HUMAN DATA ANALYSIS =====

# ----- Load Raw Human Judgment Data -----
df_s <- read_csv("data/s.csv") |>
  transmute(
    subject, format,trial,
    true_posterior,response,
    BR = br,
    HR = hr,
    FAR = far,
    subject_s = paste0(subject, "_", ifelse(format == "frequency", 1, 2))
  )

human_dt_s <- as.data.table(df_s)

# ----- Compute Summary Statistics for Human Judgments -----

# EMP / Ad / PD metrics using response instead of predict
M_Human_s <- compute_all_metrics(
  df = human_dt_s,
  column = "response",
  group_vars = c("subject_s")
)

# Slope and intercept from regressions using response
R_Human_s <- compute_SI_by(
  dt = human_dt_s,
  group_vars = c("subject_s"),
  predictors = c("BR", "HR", "FAR"),
  column = "response"
)

# No Repeat-No Variance of response within each subject's trial set

# ----- Merge All Human Features -----
final_Human_s <- M_Human_s %>%
  left_join(R_Human_s, by = c("subject_s"))

saveRDS(final_Human_s, "data/Human_Summary_dts.rds")
saveRDS(final_Simulate_s, "data/Simulate_Summary_dts.rds")

# ===== Experiment =====
# ----- Compute Accuracy (CLC) -----
A_Simulate <- compute_CLC_summary(dt)
# ========== Figure 7: Model Comparison (Accuracy) for Experiment ==========
p2<-ggplot(A_Simulate, aes(x = reorder(model, Accuracy), y = Accuracy)) +
  geom_bar(stat = "identity", fill = "#4B6C8A", width = 0.6) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Model Accuracy Comparison",
    x = "Model",
    y = "Proportion Accurate (< 3%)"
  )
ggsave("fig/Figure_10.png", p2, width = 10, height = 10, dpi = 600, units = "in")
# ----- Compute All Summary Statistics for Simulated Data -----
# EMP, Ad, PD metrics (by reference heuristics)
M_Simulate <- compute_all_metrics(dt)

# Linear regression (slope, intercept) by BR, HR, FAR
R_Simulate <- compute_SI_by(dt)

# Variance of predictions across trials
V_Simulate <- compute_variance_summary(dt)

# ----- Merge All Metrics into One Data Table -----
final_Simulate <- M_Simulate %>%
  left_join(R_Simulate, by = c("model", "Iteration")) %>%
  left_join(V_Simulate, by = c("model", "Iteration"))
# ===== HUMAN DATA ANALYSIS =====

# ----- Load Raw Human Judgment Data -----
df <- readxl::read_excel("data/df.xlsx") |>
  transmute(
    subject, format, n_trial, rt,
    BR = br,
    HR = hr,
    FAR = far,
    response = round(inputvalue / 100, 2),
    true_posterior = round(correctAnswer / 100, 2)
  )

human_dt <- as.data.table(df)

# ----- Compute Summary Statistics for Human Judgments -----

# EMP / Ad / PD metrics using response instead of predict
M_Human <- compute_all_metrics(
  df = human_dt,
  column = "response",
  group_vars = c("subject")
)

# Slope and intercept from regressions using response
R_Human <- compute_SI_by(
  dt = human_dt,
  group_vars = c("subject"),
  predictors = c("BR", "HR", "FAR"),
  column = "response"
)

# Variance of response within each subject's trial set
V_Human <- compute_variance_summary(
  dt = human_dt,
  column = "response",
  group_vars = c("subject", "BR", "HR", "FAR"),
  summary_vars = c("subject")
)
p_V_Human<-ggplot(V_Human, aes(x = mean_variance)) +
  geom_histogram(binwidth = 0.002, fill = "#4575b4", color = "black") +
  labs(
       x = "Mean Variance Across Repeated Items",
       y = "Number of Participants") +
  theme_minimal()

ggsave("fig/variance.png", p_V_Human, width = 6, height = 4.5, dpi = 300)
# ----- Merge All Human Features -----
final_Human <- M_Human %>%
  left_join(R_Human, by = c("subject")) %>%
  left_join(V_Human, by = c("subject"))

saveRDS(final_Human, "data/Human_Summary_dt.rds")
saveRDS(final_Simulate, "data/Simulate_Summary_dt.rds")
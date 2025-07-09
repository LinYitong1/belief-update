library(broom.mixed)
library(kableExtra)
library(brms)
source("R/table_helper.R") 

# ----------- Table script usage -----------

# Paths
base_path <- "data"
out_path  <- "Table"
raw_path <- file.path(base_path, "df.xlsx")
S_path   <- file.path(base_path, "s.csv")
Si_path  <- file.path(base_path, "sirota.csv")

# Load data
df  <- load_clean_data(raw_path, type = "experiment")
s   <- load_clean_data(S_path,   type = "stengard")
si  <- load_clean_data(Si_path,  type = "sirota")

# Table 1: Slope summary for Stengard
table1 <- slope_t_summary(s, mu = 1, 
                          out_path = file.path(out_path, "Table1.tex"),
                          caption = "Slope summary for Stengard.")
cat("\n=== Table 1: Slope summary for own stengard ===\n")
print(table1)
l_s<-lmer(response ~ true_posterior * format + (1 | subject), data = s)
summary(l_s)

# Table 2: Exact match/heuristic summary
table2 <- exact_match_table(s, out_path = file.path(out_path, "Table2.tex"),
                            caption = "Exact match rates across conditions.")
cat("\n=== Table 2: Heuristic/exact match summary ===\n")
print(table2)

# Table 3 :Analysis of non exact match slope
s_nm <- filter_nonmatches(s)
fit <- brm(
  response ~ true_posterior * format + (1 | subject),
  data = s_nm,
  family = gaussian(),
  iter = 2000, warmup = 500, chains = 4, cores = 4,
  seed = 123
)

fit_tbl <- tidy(fit, effects = "fixed", conf.level = 0.95)

fit_tbl <- fit_tbl %>%
  dplyr::select(term, estimate, std.error, conf.low, conf.high) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 3))) %>%
  dplyr::rename(
    Term = term,
    Estimate = estimate,
    `SE` = std.error,
    `95% CrI (low)` = conf.low,
    `95% CrI (high)` = conf.high
  )

kbl(
  fit_tbl,
  format = "latex",
  booktabs = TRUE,
  caption = "Posterior estimates from the Bayesian linear mixed-effects model (brms)."
) %>%
  kable_classic_2(full_width = FALSE, font_size = 10) %>%
  save_kable("Table/Table3.tex")

# Table 4: Round
library(dplyr)
library(tidyr)

is_round_pct <- function(x, tol = 1e-9) {
  rem <- x %% 5
  (rem < tol) | (abs(rem - 5) < tol)
}

round_prop_tbl <- function(data, label, tol = 1e-9) {
  
  heur_cols <- c("BO", "REP", "JO", "FC", "LS", "Bayes")
  
  data %>% 
    mutate(
      BO    = BR    * 100,
      REP   = HR    * 100,
      JO    = BR*HR * 100,
      FC    = (1 - FAR) * 100,
      LS    = (HR - FAR) * 100,
      Bayes = BR*HR / (BR*HR + (1 - BR)*FAR) * 100
    ) %>% 
    summarise(across(all_of(heur_cols),
                     ~ mean(is_round_pct(.x, tol), na.rm = TRUE))) %>% 
    mutate(dataset = label, .before = 1)
}


tbl_round <- bind_rows(
  round_prop_tbl(s  %>% filter(subject == "1_1"), "Standard"),
  round_prop_tbl(df %>% filter(subject == 1),     "Experiment")
)

print(tbl_round)
results <- tbl_round %>% 
  mutate(across(-dataset, ~ percent(.x, accuracy = 0.1)))
make_round_table <- function(results, out_path,
                             caption = NULL) {
  
  default_cap <- "Proportion of round-valued heuristic predictions (\\%)"
  
  kbl(results,
      format   = "latex",
      booktabs = TRUE,
      caption  = ifelse(is.null(caption), default_cap, caption),
      align    = "lrrrrrr") %>%               
    kable_classic_2(full_width = FALSE,
                    font_size   = 10) %>% 
    save_kable(out_path)
}

make_round_table(results  ,
                 out_path = "Table/Table4.tex")

# Table 5 Slope summary for your own experiment
table5 <- slope_t_summary(df, mu = 1, 
                          out_path = file.path(out_path, "Table5.tex"),
                          caption = "Slope summary for own experiment.")
cat("\n=== Table 5: Slope summary for own experiment ===\n")
print(table5)
l_df<-lmer(response ~ true_posterior * format + (1 | subject), data = df)
summary(l_df)


# Table 6 :Analysis of non exact match slope
df_nm <- filter_nonmatches(df)
fit_df <- brm(
  response ~ true_posterior * format + (1 | subject),
  data = df_nm ,
  family = gaussian(),
  iter = 2000, warmup = 500, chains = 4, cores = 4,
  seed = 123
)

fit_tbl_df <- tidy(fit_df, effects = "fixed", conf.level = 0.95)

fit_tbl_df <- fit_tbl_df %>%
  dplyr::select(term, estimate, std.error, conf.low, conf.high) %>%
  dplyr::mutate(across(where(is.numeric), ~ round(., 3))) %>%
  dplyr::rename(
    Term = term,
    Estimate = estimate,
    `SE` = std.error,
    `95% CrI (low)` = conf.low,
    `95% CrI (high)` = conf.high
  )

kbl(
  fit_tbl_df,
  format = "latex",
  booktabs = TRUE,
  caption = "Posterior estimates from the Bayesian linear mixed-effects model (brms)."
) %>%
  kable_classic_2(full_width = FALSE, font_size = 10) %>%
  save_kable("Table/Table6.tex")

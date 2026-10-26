library(broom.mixed)
library(kableExtra)
library(brms)
library(binom)
source("R/table_helper.R") 

# ----------- Stats script usage -----------

# Paths
base_path <- "data"
out_path  <- "Table"
raw_path <- file.path(base_path, "Exp1", "df.csv")
S_path   <- file.path(base_path, "Stengard", "s.csv")
Si_path  <- file.path(base_path, "sirota.csv")

df  <- load_clean_data(raw_path, type = "experiment")
s   <- load_clean_data(S_path,   type = "stengard")
si  <- load_clean_data(Si_path,  type = "sirota")

# ----------- Reanalyzing Existing Data -----------
#Conservatism
slope_s <- slope_t_summary(s, mu = 1, 
                          out_path = file.path(out_path, "Table_slope_s.tex"),
                          caption = "Slope summary for Stengard.")
cat("\n=== Table: Slope summary for Stengard ===\n")
print(slope_s)

l_s<-lmer(response ~ true_posterior * format + (1 | subject), data = s)
summary(l_s)
#Chi-squared analysis: Heuristic distribution × Format
chi_s <- run_chi_test(
  data = s,
  dataset_name = "Stengard",
  seed = 111
)
print(chi_s$summary)

#Non-match trials
s_nm <- filter_nonmatches(s)
filtered_s <- s_nm  %>%
  group_by(subject) %>%
  filter(n() >= 3, n_distinct(true_posterior) > 1) %>%
  ungroup()

slope_s_non_0 <-slope_t_summary(filtered_s, mu = 0, 
                                caption = "Slope summary of non exact match for Stengard (0).")
print(slope_s_non_0)

slope_s_non <-slope_t_summary(filtered_s, mu = 1, 
                         caption = "Slope summary of non exact match for Stengard (1).")
print(slope_s_non)

l_s_n<-lmer(response ~ true_posterior * format + (1 | subject), data = filtered_s)
summary(l_s_n)

# ----------- Experiment -----------
#Conservatism
slope_df <- slope_t_summary(df, mu = 1, 
                           caption = "Slope summary for Experiment.")
cat("\n=== Table: Slope summary for Experiment ===\n")
print(slope_df)

l_df<-lmer(response ~ true_posterior * format + (1 | subject), data = df)
summary(l_df)
#Chi-squared analysis: Heuristic distribution × Format
chi_df <- run_chi_test(
  data = df,
  dataset_name = "Experiment 1",
  seed = 222
)
print(chi_df$summary)
#Non-match trials
df_nm <- filter_nonmatches(df)
filtered_df <- df_nm  %>%
  group_by(subject) %>%
  filter(n() >= 3, n_distinct(true_posterior) > 1) %>%
  ungroup()

slope_df_non_0 <-slope_t_summary(filtered_df, mu = 0, 
                               caption = "Slope summary of non exact match for Experiment.(0)")
print(slope_df_non_0)

slope_df_non <-slope_t_summary(filtered_df, mu = 1, 
                              caption = "Slope summary of non exact match for Experiment.(1)")
print(slope_df_non)

l_df_n<-lmer(response ~ true_posterior * format + (1 | subject), data = filtered_df)
summary(l_df_n)


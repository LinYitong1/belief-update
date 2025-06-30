
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
table1 <- slope_t_summary(df, mu = 1, out_path = file.path(out_path, "Table1.tex"),
                          caption = "Slope summary for Stengard.")
cat("\n=== Table 1: Slope summary for own experiment ===\n")
print(table1)

# Table 2: Exact match/heuristic summary
table2 <- exact_match_table(s, si, out_path = file.path(out_path, "Table2.tex"),
                            caption = "Exact match rates across conditions.")
cat("\n=== Table 2: Heuristic/exact match summary ===\n")
print(table2)


# Table 3: Slope summary for your own experiment
table3 <- slope_t_summary(df, mu = 1, out_path = file.path(out_path, "Table3.tex"),
                          caption = "Slope summary for own experiment.")
cat("\n=== Table 1: Slope summary for own experiment ===\n")
print(table3)

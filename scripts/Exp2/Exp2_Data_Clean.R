# Exp2 Data Clean ------------------------------------------------------------
# Consolidated pipeline: raw jsPsych sessions -> cleaned, anonymised data.
# Inputs : data/Exp2/exp_2_raw/*.csv; Prolific demographic export (optional).
# Outputs: data/Exp2/Exp_2.csv (19 columns, anonymised),
#          data/Exp2/Exp2_Demographics_by_subject.csv,
#          data/Exp2/subject_id_key_PRIVATE.csv (do not share).
# Run    : Rscript scripts/Exp2/Exp2_Data_Clean.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
})

# ── 00. Config ──────────────────────────────────────────────────────────────

source_dir   <- "data/Exp2/exp_2_raw"
output_file  <- "data/Exp2/Exp_2.csv"
key_file     <- "data/Exp2/subject_id_key_PRIVATE.csv"
demo_file    <- "~/Downloads/prolific_demographic_export_6a9a92c456d47512aa4522ea.csv"
demo_output  <- "data/Exp2/Exp2_Demographics_by_subject.csv"

raw_id_pattern <- "^[0-9a-f]{24}$"

to_num <- function(x) suppressWarnings(as.numeric(x))
to_lgl <- function(x) {
  x <- str_to_lower(str_trim(as.character(x)))
  case_when(
    x %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

NUM_COLS <- c(
  "rt", "time_elapsed", "speededDeadlineMs", "deadline_ms",
  "inputvalue", "br", "hr", "far",
  "error", "correctAnswer", "correctCount", "n_trial",
  "repeat_number", "block_num", "trial_in_block",
  "global_main_trial"
)
LGL_COLS <- c(
  "speededFirst", "speeded", "timed_out", "invalid_input", "correct",
  "practice", "counts_for_bonus"
)
# Raw columns read for cleaning (completion checks, ID resolution, RT
# reconstruction, block_type / is_check_trial derivation).
KEEP_COLS <- c(
  "trial_type", "response", "stimulus", "format", "subject_id",
  "blockType", NUM_COLS, LGL_COLS
)
# Columns written to Exp_2.csv; anything else is rebuilt by re-running.
OUTPUT_COLS <- c(
  "subject_id", "format", "block_type", "block_num", "n_trial",
  "trial_in_block", "file_id", "source_row",
  "br", "hr", "far", "inputvalue", "correctAnswer",
  "analysis_rt_ms", "analysis_rt_sec", "analysis_rt_source",
  "practice", "is_check_trial", "invalid_input"
)

show <- function(title, x) {
  cat("\n===== ", title, " =====\n", sep = "")
  print(x, n = if (inherits(x, "tbl_df")) Inf else NULL)
}

# ── 01. Trial-level cleaning ────────────────────────────────────────────────

resolve_sid <- function(dat, filename_id) {
  sid <- c(
    if ("subject_id" %in% names(dat)) dat$subject_id,
    if ("Participant id" %in% names(dat)) dat[["Participant id"]]
  )
  sid <- trimws(as.character(sid))
  sid <- sid[!is.na(sid) & nzchar(sid)]
  if (length(sid)) return(sid[1])
  
  # Fallback: Prolific ID in survey JSON (cc11)
  if (all(c("trial_type", "response") %in% names(dat))) {
    cc11 <- dat |>
      filter(str_detect(coalesce(trial_type, ""), "(?i)survey"), !is.na(response)) |>
      pull(response) |>
      map_chr(\(x) {
        out <- tryCatch(jsonlite::fromJSON(x), error = \(e) NULL)
        if (is.null(out) || is.null(out$cc11)) NA_character_ else as.character(out$cc11)
      }) |>
      str_extract("[a-f0-9]{24}") |>
      purrr::discard(is.na)
    if (length(cc11)) return(cc11[1])
  }
  as.character(filename_id)
}

process_csv <- function(file, file_id) {
  filename_id <- suppressWarnings(as.numeric(sub(
    "file_(\\d+)\\.csv", "\\1", basename(file)
  )))
  dat <- read_csv(
    file, col_types = cols(.default = col_character()),
    na = c("", "NA", "NaN", "nan", "null", "NULL"),
    show_col_types = FALSE, progress = FALSE
  )
  sid <- resolve_sid(dat, filename_id)
  dat[setdiff(KEEP_COLS, names(dat))] <- NA_character_
  
  dat |>
    dplyr::select(all_of(KEEP_COLS)) |>
    mutate(
      source_row = row_number(),
      across(all_of(NUM_COLS), to_num),
      across(all_of(LGL_COLS), to_lgl),
      subject_id = sid,
      file_id = as.character(file_id),
      rt_sec = rt / 1000,
      time_elapsed_sec = time_elapsed / 1000,
      deadline_sec = coalesce(deadline_ms, speededDeadlineMs) / 1000,
      block_type = case_when(
        speeded %in% TRUE ~ "Time-limited",
        speeded %in% FALSE ~ "Self-paced",
        TRUE ~ NA_character_
      ),
      subject_block = paste(subject_id, block_type, sep = "__"),
      is_check_trial = (counts_for_bonus %in% TRUE & br == 1 & hr == 1 & far == 0) |
        (br == 1 & hr == 1 & far == 0 & is.na(n_trial))
    )
}

add_reconstructed_rt <- function(df) {
  df |>
    arrange(subject_id, file_id, source_row) |>
    group_by(subject_id, file_id) |>
    mutate(
      raw_elapsed_gap_ms = time_elapsed - lag(time_elapsed),
      actual_rt_ms = raw_elapsed_gap_ms -
        if_else(!is.na(n_trial) & trial_in_block == 1, 500, 0),
      actual_rt_sec = actual_rt_ms / 1000
    ) |>
    ungroup()
}

add_analysis_flags <- function(df) {
  df |>
    mutate(
      deadline_ms_used = coalesce(deadline_ms, speededDeadlineMs),
      analysis_rt_ms = case_when(
        !is.na(actual_rt_ms) & actual_rt_ms > 0 ~ actual_rt_ms,
        block_type == "Self-paced" & !is.na(rt) & rt > 0 ~ rt,
        block_type == "Time-limited" & !is.na(rt) & rt > 0 &
          !is.na(deadline_ms_used) & rt <= deadline_ms_used ~ rt,
        TRUE ~ NA_real_
      ),
      analysis_rt_sec = analysis_rt_ms / 1000,
      analysis_rt_source = case_when(
        !is.na(actual_rt_ms) & actual_rt_ms > 0 ~ "reconstructed",
        block_type == "Self-paced" & !is.na(rt) & rt > 0 ~ "original_rt_fallback",
        block_type == "Time-limited" & !is.na(rt) & rt > 0 &
          !is.na(deadline_ms_used) & rt <= deadline_ms_used ~ "original_rt_fallback",
        TRUE ~ "missing"
      ),
      valid_numeric_response = !is.na(inputvalue) & between(inputvalue, 0, 100),
      submitted_enter = coalesce(
        str_to_lower(str_trim(as.character(response))) == "enter", FALSE
      ),
      within_response_window = case_when(
        block_type == "Self-paced" ~ TRUE,
        block_type == "Time-limited" & !is.na(analysis_rt_ms) &
          !is.na(deadline_ms_used) ~ analysis_rt_ms <= deadline_ms_used,
        TRUE ~ FALSE
      ),
      corrected_timed_out = case_when(
        block_type == "Time-limited" ~ !(submitted_enter & within_response_window),
        block_type == "Self-paced" ~ FALSE,
        TRUE ~ NA
      ),
      analysis_valid_response = case_when(
        block_type == "Self-paced" ~ valid_numeric_response,
        block_type == "Time-limited" ~
          valid_numeric_response & submitted_enter & within_response_window,
        TRUE ~ FALSE
      )
    )
}

files <- sort(list.files(source_dir, "\\.csv$", full.names = TRUE, ignore.case = TRUE))
if (!length(files)) stop("No CSV files in ", source_dir)
cat("Found", length(files), "files in", source_dir, "\n")

combined_p4_rt <- files |>
  imap(\(f, i) process_csv(f, i)) |>
  bind_rows() |>
  add_reconstructed_rt() |>
  add_analysis_flags()

# ── 02. Participant inclusion ───────────────────────────────────────────────

completed_p4 <- combined_p4_rt |>
  filter(str_detect(coalesce(stimulus, ""), fixed(
    "Thank you for participating in this experiment!", ignore_case = TRUE
  ))) |>
  distinct(subject_id) |>
  pull(subject_id)
passed_p4 <- combined_p4_rt |>
  filter(is_check_trial, inputvalue == 100) |>
  distinct(subject_id) |>
  pull(subject_id)
eligible_p4 <- intersect(completed_p4, passed_p4)

cat("\nTotal:", n_distinct(combined_p4_rt$subject_id),
    "| Completed:", length(completed_p4),
    "| Passed catch:", length(passed_p4),
    "| Completed + passed:", length(eligible_p4), "\n")

session_qc <- combined_p4_rt |>
  filter(subject_id %in% eligible_p4, !is.na(n_trial), !is_check_trial) |>
  group_by(subject_id, block_type) |>
  summarise(
    n_trials = n(),
    n_missing = sum(is.na(inputvalue)),
    n_deadline_text = sum(str_detect(coalesce(stimulus, ""),
                                     fixed("Time remaining", ignore_case = TRUE))),
    n_invalid_reconstructed_rt = sum(!is.na(actual_rt_ms) & actual_rt_ms <= 0),
    n_tiny_positive_rt = sum(!is.na(actual_rt_ms) & actual_rt_ms > 0 & actual_rt_ms < 100),
    .groups = "drop"
  )
qc_failed <- session_qc |>
  filter(block_type == "Self-paced", n_missing > 0 | n_deadline_text > 0)
corrupted_ids <- unique(qc_failed$subject_id)
final_p4 <- setdiff(eligible_p4, corrupted_ids)
show("Session QC failures", qc_failed)
cat("Excluded corrupted:", length(corrupted_ids),
    "| Final subjects:", length(final_p4), "\n")

main_trials <- combined_p4_rt |>
  filter(subject_id %in% final_p4, !is.na(n_trial), !is_check_trial)
mf_p4 <- main_trials |>
  filter(!is.na(block_type), analysis_valid_response,
         !is.na(br), !is.na(hr), !is.na(far))

# ── 03. Anonymisation ───────────────────────────────────────────────────────
# Raw 24-hex subject_id -> integer code; existing key numbering is preserved,
# new subjects are appended.

raw_ids_all <- sort(unique(mf_p4$subject_id[!is.na(mf_p4$subject_id)]))
stopifnot(all(grepl(raw_id_pattern, raw_ids_all)))

if (file.exists(key_file)) {
  key <- read_csv(key_file, show_col_types = FALSE) |>
    mutate(subject = as.integer(subject), subject_id = as.character(subject_id))
} else {
  key <- tibble(subject = integer(), subject_id = character())
}

new_ids <- setdiff(raw_ids_all, key$subject_id)
if (length(new_ids)) {
  next_no <- if (nrow(key)) max(key$subject) + 1L else 1L
  key <- bind_rows(
    key,
    tibble(subject = next_no + seq_along(new_ids) - 1L, subject_id = new_ids)
  )
  key <- key |> arrange(subject)
  cat("Added", length(new_ids), "new subject(s) to the key\n")
}
id_map <- setNames(as.character(key$subject), key$subject_id)
stopifnot(!anyDuplicated(key$subject_id), !anyDuplicated(key$subject))

write_csv(key, key_file)

mf_anon <- mf_p4 |>
  mutate(
    subject_id = unname(id_map[subject_id]),
    subject_block = paste(subject_id, block_type, sep = "__")
  ) |>
  dplyr::select(all_of(OUTPUT_COLS))

# ── 04. Save main data ──────────────────────────────────────────────────────

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(mf_anon, output_file)
cat("\nSaved:", output_file, "| rows:", nrow(mf_anon),
    "| cols:", ncol(mf_anon),
    "| participants:", n_distinct(mf_anon$subject_id), "\n")
show("Trial summary", mf_anon |>
       dplyr::count(subject_id, format, block_type, name = "n_trials") |>
       group_by(format, block_type) |>
       summarise(n_subjects = n(), median = median(n_trials),
                 min = min(n_trials), max = max(n_trials), .groups = "drop"))

# ── 05. Demographics ────────────────────────────────────────────────────────
# Match the final sample to the Prolific export via the private key; only the
# anonymous subject code is written to the output.

final_sample <- mf_anon |>
  distinct(subject_id, format) |>
  mutate(subject_id = as.character(subject_id)) |>
  left_join(
    key |> rename(prolific_id = subject_id) |>
      transmute(subject_id = as.character(subject), prolific_id),
    by = "subject_id"
  )

if (file.exists(demo_file)) {
  demo <- read_csv(path.expand(demo_file), show_col_types = FALSE) |>
    rename(prolific_id = `Participant id`) |>
    mutate(
      prolific_id = as.character(prolific_id),
      # Prolific uses DATA_EXPIRED / CONSENT_REVOKED placeholders for withheld fields
      across(c(Age, Sex, `Ethnicity simplified`), ~ na_if(.x, "CONSENT_REVOKED")),
      across(c(Age, Sex, `Ethnicity simplified`), ~ na_if(.x, "DATA_EXPIRED")),
      Age = suppressWarnings(as.numeric(Age)),
      Sex = str_to_title(str_trim(as.character(Sex)))
    ) |>
    distinct(prolific_id, .keep_all = TRUE)
  
  n_matched <- sum(final_sample$prolific_id %in% demo$prolific_id)
  cat("\nMatched to demographic export:", n_matched, "of", nrow(final_sample), "\n")
  if (n_matched < nrow(final_sample)) {
    warning("Unmatched participants: ",
            paste(setdiff(final_sample$prolific_id, demo$prolific_id), collapse = ", "))
  }
  
  demo_final <- final_sample |>
    left_join(
      demo |> dplyr::select(
        prolific_id, Age, Sex, `Ethnicity simplified`,
        `Country of residence`, `Student status`, `Employment status`
      ),
      by = "prolific_id"
    ) |>
    dplyr::select(-prolific_id)
  
  write_csv(demo_final, demo_output)
  cat("Saved:", demo_output, "| rows:", nrow(demo_final), "\n")
  
  show("Age by format", demo_final |>
         group_by(format) |>
         summarise(n = sum(!is.na(Age)), mean_age = mean(Age, na.rm = TRUE),
                   sd_age = sd(Age, na.rm = TRUE), min_age = min(Age, na.rm = TRUE),
                   max_age = max(Age, na.rm = TRUE), .groups = "drop"))
  cat("\nOverall age: M =", round(mean(demo_final$Age, na.rm = TRUE), 1),
      ", SD =", round(sd(demo_final$Age, na.rm = TRUE), 1),
      ", range =", min(demo_final$Age, na.rm = TRUE), "-",
      max(demo_final$Age, na.rm = TRUE), "\n")
  show("Sex by format", demo_final |> dplyr::count(format, Sex))
} else {
  cat("\nDemographic export not found (", demo_file, "); skipping.\n", sep = "")
}

# ── 06. Diagnostics (console only) ──────────────────────────────────────────

catch_by_subject <- combined_p4_rt |>
  filter(subject_id %in% completed_p4, is_check_trial) |>
  group_by(subject_id, format) |>
  summarise(catch_pass = any(inputvalue == 100, na.rm = TRUE), .groups = "drop")
catch_summary <- catch_by_subject |>
  dplyr::count(format, catch_pass, name = "n") |>
  dplyr::group_by(format) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::ungroup()
show("Catch-trial pass/fail by format", catch_summary)
catch_table <- table(catch_by_subject$format, catch_by_subject$catch_pass)
if (all(dim(catch_table) == c(2L, 2L))) {
  cat("\nFisher test for catch-trial pass by format:\n")
  print(fisher.test(catch_table))
} else {
  cat("\nFisher test skipped: 2x2 catch-trial table not available.\n")
}

reconstruction_qc <- main_trials |>
  summarise(
    n_main_trials = n(),
    n_negative_reconstructed_rt = sum(actual_rt_ms < 0, na.rm = TRUE),
    n_zero_reconstructed_rt = sum(actual_rt_ms == 0, na.rm = TRUE),
    n_nonpositive_reconstructed_rt = sum(actual_rt_ms <= 0, na.rm = TRUE),
    n_original_rt_fallback = sum(analysis_rt_source == "original_rt_fallback"),
    n_missing_analysis_rt = sum(is.na(analysis_rt_ms))
  )
show("RT reconstruction QC", reconstruction_qc)

tl_trials <- main_trials |> filter(block_type == "Time-limited")
show("Original vs corrected timeout", tl_trials |>
       dplyr::count(timed_out, corrected_timed_out, name = "n_trials"))
show("Timeout summary", tl_trials |>
       summarise(total_time_limited = n(),
                 n_original_timed_out = sum(timed_out %in% TRUE),
                 pct_original_timed_out = 100 * n_original_timed_out / total_time_limited,
                 n_corrected_timed_out = sum(corrected_timed_out %in% TRUE),
                 pct_corrected_timed_out = 100 * n_corrected_timed_out / total_time_limited))
show("Effect of timeout correction on numeric responses", main_trials |>
       filter(!is.na(block_type), valid_numeric_response) |>
       group_by(format, block_type) |>
       summarise(n_numeric_responses = n(),
                 n_original_timed_out = sum(timed_out %in% TRUE),
                 n_corrected_timed_out = sum(corrected_timed_out %in% TRUE),
                 n_original_valid_for_analysis = sum(!(timed_out %in% TRUE)),
                 n_corrected_valid_for_analysis = sum(analysis_valid_response),
                 change_in_valid_trials = n_corrected_valid_for_analysis - n_original_valid_for_analysis,
                 .groups = "drop"))

rt_validation <- main_trials |>
  filter(!is.na(rt), !is.na(actual_rt_ms), actual_rt_ms > 0,
         block_type == "Self-paced" |
           (block_type == "Time-limited" & !is.na(deadline_ms_used) &
              rt <= deadline_ms_used)) |>
  summarise(n = n(), correlation = cor(rt, actual_rt_ms, use = "complete.obs"),
            mean_difference_ms = mean(rt - actual_rt_ms, na.rm = TRUE),
            median_difference_ms = median(rt - actual_rt_ms, na.rm = TRUE))
show("RT reconstruction validation", rt_validation)

rt_p4 <- mf_anon |> filter(!is.na(analysis_rt_ms), analysis_rt_ms >= 100)
show("RT QC by format and condition", mf_anon |>
       group_by(format, block_type) |>
       summarise(n_judgment_trials = n(),
                 n_reconstructed = sum(analysis_rt_source == "reconstructed"),
                 n_fallback = sum(analysis_rt_source == "original_rt_fallback"),
                 n_missing_rt = sum(is.na(analysis_rt_ms)),
                 n_under_100ms = sum(!is.na(analysis_rt_ms) & analysis_rt_ms > 0 & analysis_rt_ms < 100),
                 n_rt_trials_retained = sum(!is.na(analysis_rt_ms) & analysis_rt_ms >= 100),
                 .groups = "drop"))
show("Response time summary", rt_p4 |>
       group_by(format, block_type) |>
       summarise(n_subjects = n_distinct(subject_id), n_trials = n(),
                 mean_rt = mean(analysis_rt_sec), sd_rt = sd(analysis_rt_sec),
                 median_rt = median(analysis_rt_sec),
                 q25 = quantile(analysis_rt_sec, 0.25),
                 q75 = quantile(analysis_rt_sec, 0.75), .groups = "drop"))

rt_subject <- rt_p4 |>
  group_by(subject_id, format, block_type) |>
  summarise(mean_rt = mean(analysis_rt_sec), median_rt = median(analysis_rt_sec),
            n_trials = n(), .groups = "drop")
show("Participant-level RT summary", rt_subject |>
       group_by(format, block_type) |>
       summarise(n_subjects = n(), group_mean_rt = mean(mean_rt),
                 group_sd_rt = sd(mean_rt), group_se_rt = sd(mean_rt)/sqrt(n()),
                 group_median_rt = median(mean_rt), .groups = "drop"))

rt_wide <- rt_subject |>
  dplyr::select(subject_id, format, block_type, mean_rt) |>
  pivot_wider(names_from = block_type, values_from = mean_rt) |>
  drop_na(`Self-paced`, `Time-limited`)
rt_test <- rt_wide |>
  group_by(format) |>
  group_modify(\(.x, .y) {
    test <- t.test(.x$`Time-limited`, .x$`Self-paced`, paired = TRUE)
    tibble(n = nrow(.x), mean_self_paced = mean(.x$`Self-paced`),
           mean_time_limited = mean(.x$`Time-limited`),
           mean_difference = mean(.x$`Time-limited` - .x$`Self-paced`),
           t = unname(test$statistic), df = unname(test$parameter), p = test$p.value)
  }) |>
  ungroup()
show("Paired RT manipulation check", rt_test)

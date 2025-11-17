# =========================================================
# Belief Updating Project — Data Clean
# Author: Yitong
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# ---------------------------------------------------------
# Section A: Copy/Rename raw files
# ---------------------------------------------------------
library(tidyverse)

original_folder <- "data/data_raw"
new_folder      <- "data/renamed_data"

dir.create(new_folder, showWarnings = FALSE, recursive = TRUE)

# Only CSVs (case-insensitive), with full paths
files <- list.files(original_folder, pattern = "(?i)\\.csv$", full.names = TRUE)
n <- length(files)
cat("Found", n, "CSV files.\n")

# Generate exactly n new names, not 1:200
new_file_names <- file.path(new_folder, sprintf("file_%03d.csv", seq_len(n)))

# COPY (keeps originals)
ok <- file.copy(from = files, to = new_file_names, overwrite = FALSE)

cat("Copied:", sum(ok), " | Failed:", sum(!ok), "\n")
if (any(!ok)) {
  cat("First few failures:\n")
  print(tibble(from = files[!ok], to = new_file_names[!ok]) %>% head())
}


# ---------------------------------------------------------
# Section B: Scan raw files for pay_1  pay_15 TRUE
# ---------------------------------------------------------

to_num <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(x))
}


to_num <- function(x) { if (is.numeric(x)) x else suppressWarnings(as.numeric(x)) }

# extract 24-hex ID from text like "... 583730c8b647c40001043cbf ..."
extract_hex24 <- function(txt) {
  if (is.na(txt) || !nzchar(txt)) return(NA_character_)
  m <- stringr::str_match(txt, "([a-f0-9]{24})")
  if (is.na(m[1,2])) NA_character_ else m[1,2]
}

process_csv <- function(file_path) {
  file_name <- basename(file_path)
  subject_from_filename <- as.numeric(sub("file_(\\d+)\\.csv", "\\1", file_name))
  
  # Read: include potential ID columns; ignore if some are absent
  dat <- tryCatch(
    readr::read_csv(
      file_path,
      col_types = readr::cols_only(
        trial_type       = readr::col_character(),
        response         = readr::col_character(),
        format           = readr::col_character(),
        rt               = readr::col_guess(),
        inputvalue       = readr::col_guess(),
        br               = readr::col_guess(),
        hr               = readr::col_guess(),
        far              = readr::col_guess(),
        brr              = readr::col_guess(),
        error            = readr::col_guess(),
        correctAnswer    = readr::col_guess(),
        n_trial          = readr::col_guess(),
        subject_id       = readr::col_character(),
        `Participant id` = readr::col_character()
      ),
      na = c("", "NA", "NaN", "nan", "null", "NULL"),
      show_col_types = FALSE
    ),
    error = function(e) readr::read_csv(
      file_path,
      na = c("", "NA", "NaN", "nan", "null", "NULL"),
      show_col_types = FALSE
    )
  )
  
  # Prefer subject_id, then Participant id
  sid_file <- c(
    if ("subject_id" %in% names(dat)) dat$subject_id,
    if ("Participant id" %in% names(dat)) dat$`Participant id`
  ) |> as.character() |> trimws()
  sid_file <- sid_file[!is.na(sid_file) & sid_file != ""]
  sid_file <- if (length(sid_file) == 0) NA_character_ else sid_file[1]
  
  id_source <- "subject_id_or_participant"
  
  # If still missing, try cc11 inside survey response JSON
  if (is.na(sid_file) && all(c("trial_type","response") %in% names(dat))) {
    cc11_texts <- dat %>%
      dplyr::filter(stringr::str_detect(trial_type %||% "", "(?i)survey"),
                    !is.na(response), nzchar(response)) %>%
      dplyr::pull(response)
    
    if (length(cc11_texts)) {
      cc11_vals <- purrr::map_chr(cc11_texts, function(s) {
        out <- tryCatch(jsonlite::fromJSON(s, simplifyVector = TRUE), error = function(e) NULL)
        if (is.null(out)) return(NA_character_)
        v <- out[["cc11"]]
        if (is.null(v)) NA_character_ else as.character(v)
      })
      sid_from_cc11 <- purrr::map_chr(cc11_vals, extract_hex24)
      sid_from_cc11 <- sid_from_cc11[!is.na(sid_from_cc11)]
      if (length(sid_from_cc11)) {
        sid_file <- sid_from_cc11[1]
        id_source <- "cc11_extracted"
      }
    }
  }
  
  # Final fallback to filename index
  if (is.na(sid_file)) {
    sid_file <- as.character(subject_from_filename)
    id_source <- "filename_fallback"
  }
  
  # Keep analysis columns safely
  keep <- intersect(
    c("format","rt","inputvalue","br","hr","far","brr","error","correctAnswer","n_trial"),
    names(dat)
  )
  dat <- dplyr::select(dat, dplyr::all_of(keep))
  
  # Coerce numeric columns
  num_cols <- intersect(c("rt","inputvalue","br","hr","far","brr","error","correctAnswer","n_trial"), names(dat))
  if (length(num_cols)) dat <- dplyr::mutate(dat, dplyr::across(dplyr::all_of(num_cols), to_num))
  
  # Propagate ID to all rows, then filter trials
  dplyr::mutate(dat,
                subject_id = sid_file,
                id_source  = id_source,
                subject    = subject_from_filename
  ) %>%
    dplyr::filter(!is.na(inputvalue) & !is.na(error))
}


folder_path <- "data/renamed_data"
file_names  <- list.files(path = folder_path, pattern = "(?i)\\.csv$", full.names = TRUE)
processed_data <- lapply(file_names, process_csv)
combined_data <- dplyr::bind_rows(processed_data, .id = "file_id") %>%
  dplyr::filter(inputvalue <= 100)


# ---------------------------------------------------------
# Section C: Count/filter/remove subjects based on check
# ---------------------------------------------------------
# NOTE: Your comment said "Count rows where brr=1 and inputvalue=100"
# but the code filtered inputvalue != 100. I did not change your logic line,
# just kept it as-is to respect your original code, but be aware of the mismatch.
count <- combined_data %>%
  filter(brr == 1 & inputvalue != 100) %>%
  nrow()
#23
# Subjects where brr=1 & inputvalue == 100
subjectst <- combined_data %>%
  filter(brr == 1 & inputvalue == 100) %>%
  pull(subject_id)

# keep those subjects from combined_data (fixed undefined variable 'subjects' -> 'subjectst')
dddf <- combined_data %>%
  filter(subject_id %in% subjectst)

mf <- dddf %>% filter(!is.na(n_trial))
write_csv(mf, "data/df.csv")
cat("Saved experiment_data.csv\n")

# ---------------------------------------------------------
# Section D:Demographics for participants who PASSED the check
#            Pass condition: brr == 1 & inputvalue == 100
#            Join with Prolific export that contains demographics.
# ---------------------------------------------------------
# 1) Identify passed subjects (from combined_data built above)
subjects_passed <- combined_data %>%
  filter(brr == 1, inputvalue == 100, !is.na(subject_id)) %>%
  distinct(subject_id) 

# 2) Read Prolific demographics CSV (the file that contains demographic fields)
#    Replace the path below with your actual Prolific export CSV path.
prolific_demographic_csv <- "data/demographic.csv"

if (file.exists(prolific_demographic_csv)) {
  demo_raw <- read_csv(prolific_demographic_csv, show_col_types = FALSE)
  
  # Prolific file uses "Participant id" as identifier — align to subject_id
  # Convert to character for safe joining
  demo <- demo_raw %>%
    rename(subject_id = `Participant id`) %>%
    mutate(subject_id = as.character(subject_id))
  
  # Our 'subjects_passed' are numeric from filenames; coerce to character for join
  subjects_passed <- subjects_passed %>%
    mutate(subject_id = as.character(subject_id))
  
  # Keep demographics only for passed subjects
  demo_passed <- demo %>%
    semi_join(subjects_passed, by = "subject_id")
  
  # Select commonly useful demographic columns (adjust if needed)
  keep_cols <- c(
    "subject_id",
    "Primary language",
    "Age",
    "Sex",
    "Ethnicity simplified",
    "Country of birth",
    "Country of residence",
    "Nationality",
    "Language",
    "Student status",
    "Employment status"
  )
  
  demo_passed_clean <- demo_passed %>% dplyr::select(any_of(keep_cols))
  cat("demographics success\n")
} else {
  cat("Failure\n")
}




demo_data_clean <- demo_passed_clean %>%
  mutate(
    Sex  = str_trim(Sex),
    Sex  = str_to_title(Sex),   # "female" -> "Female", "prefer not to say" -> "Prefer Not To Say"
    Age  = as.numeric(Age)
  )

summary_stats <- demo_data_clean %>%
  summarise(
    num_females = sum(Sex == "Female", na.rm = TRUE),
    num_males   = sum(Sex == "Male", na.rm = TRUE),
    num_pn      = sum(Sex == "Prefer Not To Say", na.rm = TRUE),
    mean_age    = mean(Age, na.rm = TRUE),
    age_span    = paste(min(Age, na.rm = TRUE), max(Age, na.rm = TRUE), sep = "–")
  )

# Print nicely
cat(paste0(
  summary_stats$num_females, " females, ",
  summary_stats$num_males, " males, ",
  summary_stats$num_pn, " prefer not to say; mean age ",
  round(summary_stats$mean_age, 1), " years, age span ",
  summary_stats$age_span, " years.\n"
))



#[1]107 females, 69 males, 1 prefer not to say; mean age 42.3 years, age span 20–76 years.
# ============================== END ==============================

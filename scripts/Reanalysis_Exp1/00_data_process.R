#Exp1 Data Process
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

library(tidyverse)

original_folder <- "data/Exp1/data_raw"
new_folder      <- "data/Exp1/renamed_data"

dir.create(new_folder, showWarnings = FALSE, recursive = TRUE)

files <- list.files(original_folder, pattern = "(?i)\\.csv$", full.names = TRUE)
n <- length(files)
cat("Found", n, "CSV files.\n")

new_file_names <- file.path(new_folder, sprintf("file_%03d.csv", seq_len(n)))

# Keep the original files and do not overwrite existing copies.
ok <- file.copy(from = files, to = new_file_names, overwrite = FALSE)

cat("Copied:", sum(ok), " | Failed:", sum(!ok), "\n")
if (any(!ok)) {
  cat("First few failures:\n")
  print(tibble(from = files[!ok], to = new_file_names[!ok]) %>% head())
}

# Read trial data and resolve participant identifiers
to_num <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(x))
}

to_num <- function(x) { if (is.numeric(x)) x else suppressWarnings(as.numeric(x)) }

# Extract a 24-character hexadecimal participant ID from text.
extract_hex24 <- function(txt) {
  if (is.na(txt) || !nzchar(txt)) return(NA_character_)
  m <- stringr::str_match(txt, "([a-f0-9]{24})")
  if (is.na(m[1,2])) NA_character_ else m[1,2]
}

process_csv <- function(file_path) {
  file_name <- basename(file_path)
  subject_from_filename <- as.numeric(sub("file_(\\d+)\\.csv", "\\1", file_name))

  # Read selected columns; fall back to reading all columns if this fails.
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

  keep <- intersect(
    c("format","rt","inputvalue","br","hr","far","brr","error","correctAnswer","n_trial"),
    names(dat)
  )
  dat <- dplyr::select(dat, dplyr::all_of(keep))

  num_cols <- intersect(c("rt","inputvalue","br","hr","far","brr","error","correctAnswer","n_trial"), names(dat))
  if (length(num_cols)) dat <- dplyr::mutate(dat, dplyr::across(dplyr::all_of(num_cols), to_num))

  # Attach the file-level ID and retain rows with a response and error value.
  dplyr::mutate(dat,
                subject_id = sid_file,
                id_source  = id_source,
                subject    = subject_from_filename
  ) %>%
    dplyr::filter(!is.na(inputvalue) & !is.na(error))
}

folder_path <- "data/Exp1/renamed_data"
file_names  <- list.files(path = folder_path, pattern = "(?i)\\.csv$", full.names = TRUE)
processed_data <- lapply(file_names, process_csv)
combined_data <- dplyr::bind_rows(processed_data, .id = "file_id") %>%
  dplyr::filter(inputvalue <= 100)

# Catch-trial counts and participant inclusion
count <- combined_data %>%
  filter(brr == 1 & inputvalue != 100) %>%
  nrow()

subjectst <- combined_data %>%
  filter(brr == 1 & inputvalue == 100) %>%
  pull(subject_id)

# Retain participants with a passing catch trial (previously corrected to use 'subjectst').
dddf <- combined_data %>%
  filter(subject_id %in% subjectst)

mf <- dddf %>% filter(!is.na(n_trial))
write_csv(mf, "data/Exp1/df.csv")
cat("Saved data/Exp1/df.csv\n")

subjects_passed <- combined_data %>%
  filter(brr == 1, inputvalue == 100, !is.na(subject_id)) %>%
  distinct(subject_id)

# The demographic export is expected at the path below.
prolific_demographic_csv <- "data/Exp1/demographic.csv"

if (file.exists(prolific_demographic_csv)) {
  demo_raw <- read_csv(prolific_demographic_csv, show_col_types = FALSE)

  # Align the Prolific participant identifier with subject_id.
  demo <- demo_raw %>%
    rename(subject_id = `Participant id`) %>%
    mutate(subject_id = as.character(subject_id))

  # Convert IDs to character for the demographic join.
  subjects_passed <- subjects_passed %>%
    mutate(subject_id = as.character(subject_id))

  demo_passed <- demo %>%
    semi_join(subjects_passed, by = "subject_id")

  # Retain available demographic fields.
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
    Sex  = str_to_title(Sex),
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

cat(paste0(
  summary_stats$num_females, " females, ",
  summary_stats$num_males, " males, ",
  summary_stats$num_pn, " prefer not to say; mean age ",
  round(summary_stats$mean_age, 1), " years, age span ",
  summary_stats$age_span, " years.\n"
))

#[1]107 females, 69 males, 1 prefer not to say; mean age 42.3 years, age span 20–76 years.

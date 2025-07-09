# =============================
# Belief Updating Table Analysis (Modular)
# Author: Yitong
# Date: 2025-06-09
# =============================

library(readr)
library(dplyr)
library(knitr)
library(kableExtra)
library(lme4)
library(readxl)

# ---- Data loading ----
load_clean_data <- function(path, type = c("experiment", "stengard", "sirota")) {
  type <- match.arg(type)
  if (type == "experiment") {
    readxl::read_excel(path) %>%
      transmute(
        subject, format, n_trial, rt,
        BR             = br,
        HR             = hr,
        FAR            = far,
        response       = round(inputvalue / 100, 2),
        true_posterior = round(correctAnswer / 100, 2)
      )
  } else if (type == "stengard") {
    readr::read_csv(path) %>%
      transmute(
        subject, format, trial,
        true_posterior, response,
        BR             = br,
        HR             = hr,
        FAR            = far,
        subject = paste0(subject, "_", ifelse(format == "frequency", 1, 2)),
      )
  } else if (type == "sirota") {
    readr::read_csv(path) %>%
      mutate(Man = factor(Man, levels = c("Probability", "Frequency"))) %>%
      transmute(
        format       = Man,
        true_posterior,
        response,
        BR, HR, FAR,
        subject_s    = Id
      )
  } else {
    stop("Unknown type")
  }
}

# ---- Regression slope/t-test summary ----
slope_t_summary <- function(df, mu = 1, out_path = NULL, caption = NULL) {
  # 1. Calculate slope per subject and format
  slopes <- df %>%
    group_by(subject, format) %>%
    do({
      fit <- lm(response ~ true_posterior, data = .)
      tibble(slope = coef(fit)[2])
    }) %>%
    ungroup()
  
  # 2. For each format, do t-test of slopes vs. mu
  results <- slopes %>%
    group_by(format) %>%
    summarise(
      mean_slope = mean(slope),
      sd_slope = sd(slope),
      n = n(),
      ci_low = mean_slope - qt(0.975, df = n - 1) * sd_slope / sqrt(n),
      ci_high = mean_slope + qt(0.975, df = n - 1) * sd_slope / sqrt(n),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      ttest = list(t.test(slopes$slope[slopes$format == format], mu = mu)),
      t = ttest$statistic,
      df = ttest$parameter,
      p = ttest$p.value,
      p_report = ifelse(p < 0.001, "< .001", sprintf("%.3f", p)),
      Mean_CI = sprintf("%.2f [%.2f, %.2f]", mean_slope, ci_low, ci_high)
    ) %>%
    ungroup() %>%
    dplyr::select(Format = format, Mean_CI, t, df, p = p_report)
  
  if (!is.null(out_path)) {
    kbl(results,
        format = "latex", booktabs = TRUE,
        caption = ifelse(is.null(caption), "Test of regression slopes vs. 1 (mean, 95\\% CI, t-test)", caption)) %>%
      kable_classic_2(full_width = FALSE, font_size = 10) %>%
      save_kable(out_path)
  }
  results
}


# ---- Heuristic/exact match summary ----
calculate_fractional_matches <- function(x) {
  pred_B  <- round(x$true_posterior, 2)
  pred_BO <- x$BR
  pred_HO <- x$HR
  pred_FO <- 1 - x$FAR
  pred_JO <- round(x$BR * x$HR, 2)
  pred_LS <- round(x$HR - x$FAR, 2)
  pred_H  <- rep(0.5, nrow(x))
  all_preds <- data.frame(
    Bayes  = pred_B,
    REP    = pred_HO,
    BO     = pred_BO,
    FC     = pred_FO,
    LS     = pred_LS,
    JO     = pred_JO,
    "50%"  = pred_H
  )
  scores <- matrix(0, nrow = nrow(x), ncol = ncol(all_preds))
  colnames(scores) <- colnames(all_preds)
  for (i in seq_len(nrow(x))) {
    matches <- which(all_preds[i, ] == x$response[i])
    if (length(matches) > 0) {
      frac <- 1 / length(matches)
      scores[i, matches] <- frac
    }
  }
  match_proportions <- colMeans(scores, na.rm = TRUE)
  round(match_proportions, 2)
}



exact_match_table <- function(stengard, sirota, out_path = NULL, caption = NULL) {
  stan <- stengard 
  stengard_matches <- rbind(
    calculate_fractional_matches(stan %>% filter(format == "probability")),
    calculate_fractional_matches(stan %>% filter(format == "frequency"))
  )
#  si_matches <- rbind(
#    calculate_fractional_matches(sirota %>% filter(format == "Probability")),
#    calculate_fractional_matches(sirota %>% filter(format == "Frequency"))
#  )
#  result <- rbind(stengard_matches, si_matches) %>% as.data.frame()
  result <- stengard_matches %>% as.data.frame()
  result <- result %>%
    mutate(Heuristic = BO + REP + FC + JO + LS + `X50.`)
  result <- result[, c("BO", "REP", "FC", "JO", "LS", "X50.", "Heuristic", "Bayes")]
  if (!is.null(out_path)) {
    kbl(result,
        format = "latex", booktabs = TRUE,
        caption = ifelse(is.null(caption), "Heuristic usage rates by participant.", caption)) %>%
      kable_classic_2(full_width = FALSE, font_size = 10) %>%
      save_kable(out_path)
  }
  result
}

# Function:Filter out responses that match any heuristic
filter_nonmatches <- function(df, tol = 0.03) {
  preds <- data.frame(
    Bayes = round(df$true_posterior, 2),
    REP   = df$HR,
    BO    = df$BR,
    FC    = 1 - df$FAR,
    JO    = round(df$BR * df$HR, 2),
    LS    = round(df$HR - df$FAR, 2),
    `50%` = 0.5
  )
  diffs <- abs(preds - df$response)
  df[rowSums(diffs <= tol) == 0, ]
}

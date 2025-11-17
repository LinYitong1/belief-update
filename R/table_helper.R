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

library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(lmerTest)
library(emmeans)
library(scales)

# ---- Data loading ----
load_clean_data <- function(path, type = c("experiment", "stengard", "sirota")) {
  type <- match.arg(type)
  if (type == "experiment") {
    readr::read_csv(path) %>%
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

# ---- Round summary ----
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
# Function:Overlap
# ---------- A. Model Prediction ----------
add_preds <- function(dt, digits = 2) {
  dt %>% mutate(
    B   = round(true_posterior, digits),
    BO  = round(BR,             digits),
    REP = round(HR,             digits),
    FC  = round(1 - FAR,        digits),
    JO  = round(BR * HR,        digits),
    LS  = round(HR - FAR,       digits),
    H   = round(0.5,            digits)
  )
}

get_overlap_mat <- function(dt, digits = 2,
                            cols = c("B","REP","BO","FC","LS","JO","H")) {
  preds <- add_preds(dt, digits) %>%
    dplyr::select(any_of(cols)) %>%
    mutate(across(everything(), as.numeric))  # guard against factors/char
  
  mod_names <- names(preds)
  n_mod <- ncol(preds)
  overlap <- matrix(0, n_mod, n_mod, dimnames = list(mod_names, mod_names))
  
  for (i in seq_len(n_mod)) {
    for (j in i:n_mod) {
      p <- mean(preds[[i]] == preds[[j]], na.rm = TRUE) * 100
      overlap[i, j] <- p
      overlap[j, i] <- p
    }
  }
  round(overlap, 0)
}

label_map <- c(B="Bayes", REP="REP", BO="BO", FC="FC", LS="LS", JO="JO", H="50%")

latex_overlap_table <- function(dt, digits=2,
                                cols=c("B","REP","BO","FC","LS","JO","H"),
                                caption="Model-prediction overlap (\\%)") {
  mat <- get_overlap_mat(dt, digits, cols)
  labs <- label_map[colnames(mat)]
  labs <- sub("%", "\\\\%", labs, fixed = TRUE)  
  colnames(mat) <- labs; rownames(mat) <- labs
  
  kbl(mat, format="latex", booktabs=TRUE, escape=FALSE,
      caption=caption, align="r") |>
    kable_styling(latex_options=c("hold_position"), full_width=FALSE)
}



# ---------- B. Core：Row level overlap & response match ----------
flag_overlaps <- function(dt, digits = 2, tol = NULL) {
  dt <- add_preds(dt, digits)
  
  pred_cols <- c("BO","REP","FC","JO","LS","H","B")  
  stopifnot(all(pred_cols %in% names(dt)))

  if (is.null(tol)) tol <- 5 * 10^-(digits + 1)
  
  preds_mat <- as.matrix(dt[, pred_cols])
  n <- nrow(preds_mat); k <- ncol(preds_mat)
  
  b_idx <- match("B", pred_cols)
  h_idx <- setdiff(seq_len(k), b_idx)
  
  bh_equal_mat <- abs(preds_mat[, h_idx, drop = FALSE] - preds_mat[, b_idx]) <= tol
  dt$stim_bh_overlap <- rowSums(bh_equal_mat) > 0
  
  r <- round(dt$response, digits)
  r_mat <- matrix(r, nrow = n, ncol = k)
  eq_to_resp <- abs(preds_mat - r_mat) <= tol
  
  bayes_hit  <- eq_to_resp[, b_idx]
  heur_hit   <- rowSums(eq_to_resp[, h_idx, drop = FALSE]) > 0
  dt$resp_bh_overlap <- bayes_hit & heur_hit

  bh_overlap_rates <- colMeans(bh_equal_mat) * 100
  names(bh_overlap_rates) <- pred_cols[h_idx]
  attr(dt, "bh_overlap_rates_pct") <- bh_overlap_rates
  
  attr(dt, "bh_union_pct") <- mean(dt$stim_bh_overlap) * 100
  
  dt
}


summarise_overlap <- function(dt_flagged){
  dt_flagged %>% 
    group_by(experiment, format) %>% 
    summarise(
      total_trials             = n(),
      # Stimulus Level
      stim_bh_overlap_n        = sum(stim_bh_overlap),
      stim_bh_overlap_rate     = stim_bh_overlap_n / total_trials,
      # Response Level
      resp_bh_overlap_n        = sum(resp_bh_overlap),
      resp_bh_overlap_rate     = resp_bh_overlap_n / total_trials,
      .groups = "drop"
    )
}

fmt_pct <- function(x, digits = 2) sprintf(paste0("%.", digits, "f%%"), 100 * x)

make_main_overlap_table <- function(dt_flagged,
                                    caption = "Bayes–Heuristic overlap (stimulus- and response-level, \\%)",
                                    digits_pct = 2) {
  sum_tbl <- dt_flagged %>%
    mutate(
      stim_bh_overlap = as.logical(stim_bh_overlap),
      resp_bh_overlap = as.logical(resp_bh_overlap)
    ) %>%
    group_by(experiment, format) %>%
    summarise(
      items_n   = n(),
      stim_hits = sum(stim_bh_overlap, na.rm = TRUE),
      resp_hits = sum(resp_bh_overlap, na.rm = TRUE),
      stim_rate = stim_hits / items_n,
      resp_rate = resp_hits / items_n,
      .groups   = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      stim_ci_low  = binom::binom.wilson(stim_hits, items_n)$lower,
      stim_ci_high = binom::binom.wilson(stim_hits, items_n)$upper,
      resp_ci_low  = binom::binom.wilson(resp_hits, items_n)$lower,
      resp_ci_high = binom::binom.wilson(resp_hits, items_n)$upper
    ) %>%
    ungroup() %>%
    mutate(
      `Items (n)`  = items_n,
      `Bayes–Heur. overlap (stim)` = fmt_pct(stim_rate, digits_pct),
      `95% CI (stim)`              = paste0(fmt_pct(stim_ci_low, digits_pct), "–", fmt_pct(stim_ci_high, digits_pct)),
      `Bayes–Heur. overlap (resp)` = fmt_pct(resp_rate, digits_pct),
      `95% CI (resp)`              = paste0(fmt_pct(resp_ci_low, digits_pct), "–", fmt_pct(resp_ci_high, digits_pct))
    ) %>%
    dplyr::select(experiment, format, `Items (n)`,
           `Bayes–Heur. overlap (stim)`, `95% CI (stim)`,
           `Bayes–Heur. overlap (resp)`, `95% CI (resp)`) %>%
    arrange(experiment, format)
  
  kbl(sum_tbl, format = "latex", booktabs = TRUE, escape = TRUE,
      caption = caption, align = c("l","l","r","r","l","r","l")) |>
    kable_styling(latex_options = c("hold_position"), full_width = FALSE)
}
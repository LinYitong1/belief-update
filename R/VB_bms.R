#' Bayesian model selection for group studies (SPM12 spm_BMS translation)
#'
#' This function is translated from spm_BMS.m of SPM12 (Stephan et al., 2009; Rigoux et al., 2014),
#' and edited to also return F0 and F1 alongside other outputs.
#'
#' @param m N x K matrix of log-model evidences (subjects x models).
#' @param n_samples Number of Monte Carlo samples for exceedance probabilities (default 1e6).
#' @return A list with elements:
#'   \item{alpha}{Posterior Dirichlet parameters.}
#'   \item{r}{Expected model frequencies.}
#'   \item{xp}{Exceedance probabilities.}
#'   \item{bor}{Bayes Omnibus Risk.}
#'   \item{F1}{Free energy under alternative (unequal frequencies).}
#'   \item{F0}{Free energy under null (equal frequencies).}
#'   \item{pxp}{Protected exceedance probabilities.}
#' @export
# ================== Core BMS Functions ==================

VB_bms_1 <- function(m, n_samples = 1e6) {
  max_val <- log(.Machine$double.xmax)
  Ni <- nrow(m); Nk <- ncol(m)
  tol <- 1e-4; alpha0 <- rep(1, Nk); alpha <- alpha0
  g <- matrix(0, Ni, Nk)
  diff <- Inf
  while (diff > tol) {
    for (i in seq_len(Ni)) {
      log_u <- m[i, ] + digamma(alpha) - digamma(sum(alpha))
      log_u <- log_u - mean(log_u)
      log_u <- sign(log_u) * pmin(abs(log_u), max_val)
      g[i, ] <- exp(log_u) / sum(exp(log_u))
    }
    prev <- alpha
    alpha <- alpha0 + colSums(g)
    diff <- sqrt(sum((alpha - prev)^2))
  }
  r <- alpha / sum(alpha)
  xp <- dirichlet_exceedance(alpha, n_samples)
  F1 <- FE(m, alpha, g, alpha0)
  F0 <- FE_null(m)
  bor <- 1 / (1 + exp(F1 - F0))
  pxp <- (1 - bor) * xp + bor / Nk
  list(alpha = alpha, r = r, xp = xp, bor = bor, F1 = F1, F0 = F0, pxp = pxp)
}

dirichlet_exceedance <- function(alpha, n_samples = 1e6) {
  K <- length(alpha)
  x <- matrix(rgamma(n_samples * K, shape = rep(alpha, each = n_samples)),
              ncol = K, byrow = TRUE)
  x <- x / rowSums(x)
  wins <- max.col(x)
  tabulate(wins, nbins = K) / n_samples
}

FE <- function(m, alpha, g, alpha0) {
  term1 <- sum(g * m)
  term2 <- -sum(g * log(g + .Machine$double.eps))
  logB0 <- sum(lgamma(alpha0)) - lgamma(sum(alpha0))
  logB1 <- sum(lgamma(alpha )) - lgamma(sum(alpha ))
  term1 + term2 + (logB1 - logB0)
}

FE_null <- function(m) {
  ev <- apply(m, 1, function(logu) {
    maxu <- max(logu)
    maxu + log(mean(exp(logu - maxu)))
  })
  sum(ev)
}

clean_model_names <- function(names) {
  names %>%
    stringr::str_remove("_NH$") %>%
    stringr::str_replace("^BS_R$", "A_BS") %>%
    stringr::str_replace("B$", "Bayes") %>%
    stringr::str_replace("M_L$", "Two Stage") %>%
    stringr::str_replace("^MIN", "Hybrid")
}

family_BMS <- function(m, family, n_samples = 1e6, aggregate = c("sum","mean")) {
  stopifnot(ncol(m) == length(family))
  aggregate <- match.arg(aggregate)
  fam_levels <- unique(family)
  lfam <- sapply(fam_levels, function(f) {
    cols <- which(family == f)
    lse  <- matrixStats::rowLogSumExps(m[, cols, drop = FALSE])
    if (aggregate == "sum") lse else lse - log(length(cols))
  })
  colnames(lfam) <- fam_levels
  fam_res <- VB_bms_1(lfam, n_samples)
  wfam <- fam_levels[which.max(fam_res$r)]
  cols_w <- which(family == wfam)
  m_sub  <- m[, cols_w, drop = FALSE]
  mod_res <- VB_bms_1(m_sub, n_samples)
  list(family_res = fam_res, model_res = mod_res, winning_family = wfam)
}

two_stage_bms <- function(m_p, model_names, fam_map, n_samples = 1e6) {
  stopifnot(ncol(m_p) == length(model_names))
  colnames(m_p) <- model_names
  family <- fam_map[model_names]
  if (any(is.na(family))) stop("Some model names not found in fam_map.")
  res <- family_BMS(m = m_p, family = family, n_samples = n_samples, aggregate = "mean")
  model_names <- colnames(m_p)
  family_vec  <- fam_map[model_names]
  family_names <- unique(family_vec)
  df_families <- tibble(family = family_names, r = as.numeric(res$family_res$r), xp = as.numeric(res$family_res$xp))
  df_models   <- tibble(model = colnames(m_p)[family == res$winning_family],
                        r = res$model_res$r, xp = res$model_res$xp)
  list(
    winning_family = res$winning_family,
    family_res     = res$family_res,
    model_res      = res$model_res,
    df_families    = df_families,
    df_models      = df_models
  )
}

fam_map <- c(
  "Bayes"            = "Bayes",
  "HO"               = "Single Heuristic",
  "BO"               = "Single Heuristic",
  "FO"               = "Single Heuristic",
  "JO"               = "Single Heuristic",
  "LS"               = "Single Heuristic",
  "H"                = "Single Heuristic",
  "MH"               = "Multiple Heuristic",
  "AH"               = "Multiple Heuristic",
  "LE"               = "Multiple Heuristic",
  "LA"               = "Linear Averaging",
  "BS"               = "Bayesian Sampler",
  "A_BS"             = "Bayesian Sampler",
  "R"                = "Random",
  "Hybrid_BS_BO"     = "Hybrid Models",
  "Hybrid_BS_HO"     = "Hybrid Models",
  "Hybrid_BS_FO"     = "Hybrid Models",
  "Hybrid_BS_JO"     = "Hybrid Models",
  "Hybrid_BS_LS"     = "Hybrid Models",
  "Hybrid_BS_H"      = "Hybrid Models",
  "Two Stage"        = "Hybrid Models"
)

# ================== Table and Plot Output Functions ==================

print_family_table <- function(res, title = "Family-level BMS") {
  df <- res$df_families %>%
    mutate(
      r_str      = sprintf("%.3f", r),
      xp_str     = sprintf("%.3f", xp),
      family_bold= if_else(family == res$winning_family,
                           paste0("**", family, "**"),
                           family)
    ) %>%
    dplyr::select(Family = family_bold, r = r_str, xp = xp_str)
  cat("###", title, "\n\n")
  df %>%
    knitr::kable(format = "markdown", align = c("l","r","r")) %>%
    kableExtra::kable_styling(full_width = FALSE) %>%
    print()
  cat("\n")
}

print_model_table <- function(res, title = "Model-level BMS") {
  df <- tibble(
    model = res$df_models$model,
    r     = res$model_res$r,
    xp    = res$model_res$xp
  ) %>%
    mutate(
      r     = sprintf("%.3f", r),
      xp    = sprintf("%.3f", xp),
      model = if_else(r == max(as.numeric(r)),
                      paste0("**", model, "**"),
                      model)
    )
  cat("###", title, "\n\n")
  df %>%
    kable(format = "markdown", align = c("l","r","r")) %>%
    kable_styling(full_width = FALSE) %>%
    print()
  cat("\n")
}

plot_two_stage_bms <- function(res, title = "Two-stage BMS: Family and Model", accent_color = "#0073B1") {
  df_fam <- res$df_families %>%
    transmute(
      name   = family,
      r      = as.numeric(r),
      level  = "Family",
      is_win = (family == res$winning_family)
    )
  df_mod <- res$df_models %>%
    transmute(
      name   = model,
      r      = as.numeric(r),
      level  = res$winning_family,
      is_win = (r == max(r))
    )
  df_all <- bind_rows(df_fam, df_mod) %>%
    mutate(
      level = factor(level, levels = c("Family", res$winning_family)),
      name  = forcats::fct_reorder(name, r)
    )
  ggplot(df_all, aes(x = r, y = name, fill = is_win)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = scales::percent(r, 0.1)),
              hjust = -0.05, size = 2, color = "black") +
    scale_fill_manual(values = c("FALSE" = "#BEDCED", "TRUE" = accent_color), guide = FALSE) +
    scale_x_continuous(labels = scales::percent_format(1), expand = expansion(c(0,0.1))) +
    facet_wrap(~level, scales = "free_y", ncol = 1, strip.position = "top") +
    labs(
      title = title,
      x     = "Group-level model probability",
      y     = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
      strip.text       = element_text(face = "bold", size = 10),
      axis.text.y      = element_text(size = 8),
      axis.title.x     = element_text(size =10),
      panel.spacing    = unit(1, "lines")
    )
}


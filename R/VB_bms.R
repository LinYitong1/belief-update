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

VB_bms_1 <- function(m, n_samples = 1e6, tol = 1e-6, max_iter = 1e4) {
  Ni <- nrow(m); Nk <- ncol(m)
  alpha0 <- rep(1, Nk); alpha <- alpha0
  g <- matrix(0, Ni, Nk)
  it <- 0L; diff <- Inf
  
  model_names <- colnames(m)
  if (is.null(model_names)) model_names <- paste0("M", seq_len(Nk))
  
  while (diff > tol && it < max_iter) {
    psi_alpha <- digamma(alpha)
    psi_alpha_sum <- digamma(sum(alpha))
    for (i in seq_len(Ni)) {
      log_u <- m[i, ] + psi_alpha - psi_alpha_sum
      log_u <- log_u - max(log_u)        
      u <- exp(log_u)
      g[i, ] <- u / sum(u)
    }
    prev <- alpha
    alpha <- alpha0 + colSums(g)
    diff <- sqrt(sum((alpha - prev)^2))
    it <- it + 1L
  }
  
  r  <- alpha / sum(alpha)
  
  names(alpha) <- model_names
  names(r)     <- model_names
  
  xp  <- dirichlet_exceedance(alpha, n_samples = n_samples)
  names(xp) <- model_names
  
  F1 <- FE(m, alpha, g, alpha0)
  F0 <- FE_null(m)
  bor <- 1 / (1 + exp(F1 - F0))
  pxp <- (1 - bor) * xp + bor / length(alpha)
  names(pxp) <- model_names
  
  list(alpha = alpha, r = r, xp = xp, bor = bor, F1 = F1, F0 = F0, pxp = pxp)
}

dirichlet_exceedance <- function(alpha, n_samples = 1e6) {
  K <- length(alpha)
  X <- matrix(
    rgamma(n_samples * K, shape = rep(alpha, each = n_samples), rate = 1),
    ncol = K, byrow = FALSE
  )
  X <- X / rowSums(X)
  tabulate(max.col(X), nbins = K) / n_samples
}

FE <- function(m, alpha, g, alpha0 = rep(1, ncol(m))) {
  # m: N x K log-evidence; g: N x K 责任; alpha: K; alpha0: K (先验Dirichlet)
  n <- nrow(m); K <- ncol(m)
  Elogr <- digamma(alpha) - digamma(sum(alpha))                 # E[log r_k]
  
  ELJ <- lgamma(sum(alpha0)) - sum(lgamma(alpha0)) + sum((alpha0 - 1) * Elogr)
  ELJ <- ELJ + sum(g * (matrix(Elogr, n, K, byrow = TRUE) + m)) # Σ_i,k g_ik(Elogr_k + m_ik)
  
  Sqf <- sum(lgamma(alpha)) - lgamma(sum(alpha)) - sum((alpha - 1) * Elogr)
  Sqm <- -sum(g * log(g + .Machine$double.eps))
  
  F1 <- ELJ + Sqf + Sqm
  return(F1)
}
FE_null <- function(m) {
  n <- nrow(m); K <- ncol(m)
  F0 <- 0
  for (i in 1:n) {
    tmp <- m[i, ] - max(m[i, ])
    g   <- exp(tmp) / sum(exp(tmp))                             # softmax
    F0  <- F0 + sum(g * (m[i, ] - log(K) - log(g + .Machine$double.eps)))
  }
  unname(F0)
}



clean_model_names <- function(names) {
  names %>%
    stringr::str_remove("_NH$") %>%
    stringr::str_replace("^BS_R$", "A_BS") %>%
    stringr::str_replace("B$", "Bayes") %>%
    stringr::str_replace("M_L$", "HA_MH") %>%
    stringr::str_replace("^MIN", "HA")
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
  
  fam_res <- VB_bms_1(lfam, n_samples = n_samples)
  
  for (nm in c("alpha","r","xp","pxp")) {
    names(fam_res[[nm]]) <- fam_levels
  }
  
  list(family_res = fam_res)
}

two_stage_bms <- function(m_p, model_names, fam_map, n_samples = 1e6,
                          aggregate = "mean", select_by = c("pxp","r")) {
  select_by <- match.arg(select_by)
  stopifnot(ncol(m_p) == length(model_names))
  colnames(m_p) <- model_names
  
  family <- fam_map[model_names]
  if (any(is.na(family))) stop("Some model names not found in fam_map.")
  
  fam <- family_BMS(m = m_p, family = family, n_samples = n_samples, aggregate = aggregate)
  fam_res <- fam$family_res
  
  pick <- if (select_by == "pxp") fam_res$pxp else fam_res$r
  if (is.null(names(pick))) stop("Internal error: family-level vector has no names.")
  wfam <- names(pick)[which.max(pick)]
  
  cols_w <- which(family == wfam)
  if (length(cols_w) == 0L) stop("Winning family has no models in fam_map/m_p. Check fam_map or names.")
  
  m_sub   <- m_p[, cols_w, drop = FALSE]
  mod_res <- VB_bms_1(m_sub, n_samples = n_samples)
  
  df_families <- tibble::tibble(
    family = names(fam_res$r),
    r      = as.numeric(fam_res$r),
    xp     = as.numeric(fam_res$xp),
    pxp    = as.numeric(fam_res$pxp)
  )
  
  df_models <- tibble::tibble(
    model = colnames(m_sub),
    r     = as.numeric(mod_res$r),
    xp    = as.numeric(mod_res$xp),
    pxp   = as.numeric(mod_res$pxp)
  )
  
  list(
    winning_family = wfam,
    family_res     = fam_res,
    model_res      = mod_res,
    df_families    = df_families,
    df_models      = df_models
  )
}

fam_map <- c(
  "Bayes"            = "Bayes",
  "HO"               = "Simple heuristics",
  "BO"               = "Simple heuristics",
  "FO"               = "Simple heuristics",
  "JO"               = "Simple heuristics",
  "LS"               = "Simple heuristics",
  "H"                = "Simple heuristics",
  "MH"               = "Complex heuristics",
  "AH"               = "Complex heuristics",
  "LE"               = "Complex heuristics",
  "LA"               = "Linear Additive",
  "BS"               = "Bayesian Sampler",
  "A_BS"             = "Bayesian Sampler",
 # "R"                = "Random",
  "HA_BS_BO"     = "Heuristic-Anchored (HA)",
  "HA_BS_HO"     = "Heuristic-Anchored (HA)",
  "HA_BS_FO"     = "Heuristic-Anchored (HA)",
  "HA_BS_JO"     = "Heuristic-Anchored (HA)",
  "HA_BS_LS"     = "Heuristic-Anchored (HA)",
  "HA_BS_H"      = "Heuristic-Anchored (HA)",
  "HA_MH"     = "Heuristic-Anchored (HA)"
)

# ================== Table and Plot Output Functions ==================

.format_pxp <- function(x) ifelse(x > 0.9995, ">.999",
                                  ifelse(x < 0.0005, "<.001", sprintf("%.3f", x)))

print_family_table <- function(res, title = "Family-level BMS") {
  df <- res$df_families %>%
    dplyr::mutate(
      r_str        = sprintf("%.3f", r),
      pxp_str      = .format_pxp(pxp),
      family_bold  = dplyr::if_else(
        as.character(family) == as.character(res$winning_family),
        paste0("**", family, "**"), as.character(family))
    ) %>%
    dplyr::select(Family = family_bold, r = r_str, pxp = pxp_str)
  
  cat("###", title, "\n\n")
  df %>%
    knitr::kable(format = "markdown", align = c("l","r","r")) %>%
    kableExtra::kable_styling(full_width = FALSE) %>%
    print()
  cat("\n")
}

print_model_table <- function(res, title = "Model-level BMS") {
  df <- tibble::tibble(
    model = res$df_models$model,
    r     = res$df_models$r,
    xp    = res$df_models$xp,
    pxp   = res$df_models$pxp
  ) %>%
    dplyr::mutate(
      r   = sprintf("%.3f", r),
      xp  = sprintf("%.3f", xp),
      pxp = .format_pxp(pxp),
      model = dplyr::if_else(
        as.numeric(r) == max(as.numeric(r)),
        paste0("**", model, "**"),
        model
      )
    )
  
  cat("###", title, "\n\n")
  df %>%
    knitr::kable(format = "markdown", align = c("l","r","r","r")) %>%
    kableExtra::kable_styling(full_width = FALSE) %>%
    print()
  cat("\n")
}

label_pm <- c(
  "Bayes" = "plain('Bayes')",
  "BS"    = "plain('BS')",
  "A_BS"  = "plain('BS-A')",
  
  "BO" = "plain('BO')",
  "HO" = "plain('HO')",
  "FO" = "plain('FO')",
  "JO" = "plain('JO')",
  "LS" = "plain('LS')",
  "H"  = "plain('H')",
  
  "MH" = "plain('MH')",
  "AH" = "plain('AH')",
  "LE" = "plain('LE')",
  "LA" = "plain('LA')",
  
  #"R"  = "plain('Random')",
  
  "HA_BS_BO" = "HA[plain('BO')]",
  "HA_BS_HO" = "HA[plain('Rep')]",
  "HA_BS_FO" = "HA[plain('FC')]",
  "HA_BS_JO" = "HA[plain('JO')]",
  "HA_BS_LS" = "HA[plain('LS')]",
  "HA_BS_H"  = "HA[plain('50%')]",

  "HA_MH" = "HA[plain('Mix')]"
)

plot_two_stage_bms <- function(res,
                               title = "Two-stage BMS: Family and Model",
                               accent_color = "#0073B1",
                               file_stub = NULL) {
  
  df_fam <- res$df_families %>%
    dplyr::transmute(
      name   = family,
      r      = as.numeric(r),
      level  = "Family",
      is_win = (family == res$winning_family)
    )
  
  df_mod <- res$df_models %>%
    dplyr::transmute(
      name   = model,
      r      = as.numeric(r),
      level  = res$winning_family,
      is_win = dplyr::near(r, max(r))
    )
  
  df_all <- dplyr::bind_rows(df_fam, df_mod) %>%
    dplyr::mutate(
      level = factor(level, levels = c("Family", res$winning_family)),
      name  = forcats::fct_reorder(name, r)
    )
  
  label_y <- function(x) {
    out <- vapply(x, function(k) {
      v <- unname(label_pm[k])
      if (is.na(v)) sprintf("plain('%s')", k) else v
    }, character(1))
    parse(text = out)
  }
  
  p <- ggplot(df_all, aes(x = r, y = name, fill = is_win)) +
    geom_col(width = 0.6) +
    geom_text(
      aes(label = scales::percent(r, accuracy = 0.1)),
      hjust = -0.15, size = 2.2, colour = "black"
    ) +
    scale_fill_manual(
      values = c("FALSE" = "#d5e5f2", "TRUE" = accent_color),
      guide  = "none"
    ) +
    scale_x_continuous(
      limits = c(0, 0.80),
      breaks = seq(0, 0.80, by = 0.20),
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0.005),
      name   = "Group-level model probability"
    ) +
    scale_y_discrete(labels = label_y) +
    facet_wrap(
      ~ level,
      scales = "free_y",
      ncol   = 1,
      strip.position = "top"
    ) +
    labs(
      title = title,
      y     = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size =10),
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text      = element_text(face = "bold", size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y     = element_text(size = 8),
      axis.text.x     = element_text(size = 8),
      axis.title.x    = element_text(size = 9),
      panel.spacing   = unit(0.8, "lines"),
      plot.margin     = margin(5.5, 5.5, 5.5, 5.5)
    )
  
  invisible(p)
} 



assign_family <- function(x) {
  dplyr::case_when(
    x %in% c("Bayes") ~ "Bayes",
    x %in% c("BS","A_BS") ~ "Bayesian Sampler",
    grepl("^HA_BS_", x) | x %in% c("HA_MH") ~ "Heuristic-Anchored (HA)",
    x %in% c("BO","HO","FO","JO","LS","H") ~ "Simple heuristics",
    x %in% c("MH","AH","LE") ~ "Complex heuristics",
    x %in% c("LA") ~ "Linear Additive",
    x %in% c("R") ~ "Random",
    TRUE ~ "Other"
  )
}

family_oob_metrics <- function(pred_model_oob, true_model, fam_map){
  if (is.factor(pred_model_oob)) pred_model_oob <- as.character(pred_model_oob)
  if (is.factor(true_model))     true_model     <- as.character(true_model)
  tf <- unname(fam_map[ true_model ])
  pf <- unname(fam_map[ pred_model_oob ])
  conf <- table(True = tf, Pred = pf)
  acc_micro <- sum(diag(conf))/sum(conf); err_micro <- 1 - acc_micro
  row_prop  <- prop.table(conf, 1)
  acc_macro <- mean(diag(row_prop), na.rm = TRUE); err_macro <- 1 - acc_macro
  list(confusion=conf, err_micro=err_micro, err_macro=err_macro,
       prior_chance = 1/length(unique(fam_map)))
}


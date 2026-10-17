# Reference predictions and response metrics
metric_reference <- function(x, reference) {
  switch(reference,
         "B"  = round(x$true_posterior, 2),
         "BO" = x$BR,
         "HO" = x$HR,
         "FO" = round(1 - x$FAR, 2),
         "JO" = round(x$BR * x$HR, 2),
         "LS" = round(x$HR - x$FAR, 2),
         "H"  = rep(0.5, nrow(x)),
         stop("Invalid reference"))
}

EMP <- function(x, reference, column = "predict") {
  mean(x[[column]] == metric_reference(x, reference))
}

Ad <- function(x, reference, column = "predict") {
  mean(abs(x[[column]] - metric_reference(x, reference)))
}

PD <- function(x, reference, column = "predict") {
  mean(sign(x[[column]] - metric_reference(x, reference)) == 1)
}

# EMP, Ad, and PD by model or participant
compute_all_metrics <- function(df, column = "predict", group_vars = c("model", "Iteration")) {
  refs <- c("B", "BO", "HO", "FO", "JO", "LS", "H")
  metrics <- df[, {
    wide <- rbindlist(lapply(refs, function(ref) {
      data.table(
        reference = ref,
        EMP = EMP(.SD, ref, column),
        Ad  = Ad(.SD, ref, column),
        PD  = PD(.SD, ref, column)
      )
    }))
    pivot_wider(wide, names_from = reference, values_from = c("EMP", "Ad", "PD"))
  }, by = group_vars]
  return(metrics)
}

# Accuracy: absolute error strictly below the threshold
compute_CLC_summary <- function(dt, column = "predict",
                                group_vars = c("model", "Iteration"),
                                summary_var = "model", threshold = 0.03) {
  clc_stat <- function(x) {
    mean(abs(x[[column]] - x$true_posterior) < threshold)
  }
  clc_dt <- dt[, .(CLC = clc_stat(.SD)), by = group_vars]
  Accuracy <- clc_dt[, .(Accuracy = mean(CLC, na.rm = TRUE)), by = summary_var]

  return(Accuracy)
}

# Slopes and transformed intercepts
SI <- function(x, intercept_cap = 100, column = "predict") {
  formula <- as.formula(paste(column, "~ true_posterior"))
  model <- lm(formula, data = x)
  slope <- coef(model)[2]
  intercept <- if (abs(1 - slope) < 1e-3) intercept_cap else coef(model)[1] / (1 - slope)
  data.frame(slope = slope, intercept = intercept)
}
# Slopes and transformed intercepts by predictor level
compute_SI_by <- function(dt, group_vars = c("model", "Iteration"),
                          predictors = c("BR", "HR", "FAR"),
                          column = "predict") {
  library(tidyr)
  si_tables <- lapply(predictors, function(varname) {
    dt[, SI(.SD, column = column), by = c(group_vars, varname)] %>%
      pivot_wider(
        names_from = !!sym(varname),
        values_from = c("slope", "intercept"),
        names_sep = paste0("_", varname)
      )
  })
  Reduce(function(x, y) left_join(x, y, by = group_vars), si_tables)
}

# Mean within-item variance
compute_variance_summary <- function(dt, column = "predict",
                                     group_vars = c("model", "Iteration", "BR", "HR", "FAR"),
                                     summary_vars = c("model", "Iteration")) {
  var_dt <- dt[, .(variance = var(get(column), na.rm = TRUE)), by = group_vars]
  var_dt[, .(mean_variance = mean(variance, na.rm = TRUE)), by = summary_vars]
}

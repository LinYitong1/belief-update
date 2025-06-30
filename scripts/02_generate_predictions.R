# -----------------------------
# Generate Predictions
# Description: Generate model predictions for all rows in parameter_dt
# Created on: 2025-06-26
# Author: Yitong Lin
# -----------------------------


library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(tidyr)

# Load data and models
source("R/models.R")
parameter_dt <- readRDS("data/parameter_dt.rds") 
parameter_dt_s <- readRDS("data/parameter_dt_s.rds") 

setDT(parameter_dt)

generate_predictions <- function(dt) {
  setDT(dt)
  pred <- dt[, .(
    B_NH  = B_NH(BR = BR, HR = HR, FAR = FAR, r = r1)$value,
    BO_NH = BO_NH(BR = BR, r = r2)$value,
    HO_NH = HO_NH(HR = HR, r = r3)$value,
    FO_NH = FO_NH(FAR = FAR, r = r4)$value,
    JO_NH = JO_NH(BR = BR, HR = HR, r = r5)$value,
    LS_NH = LS_NH(HR = HR, FAR = FAR, r = r6)$value,
    H_NH  = H_NH(r = r7)$value,
    LA_NH = LA_NH(BR = BR, HR = HR, FAR = FAR, wBR = wBR, wHR = wHR, wFAR = wFAR, r = r8)$value,
    AH_NH = select_best_heuristic_NH(BR = BR, HR = HR, FAR = FAR, r = r9),
    LE_NH = lexicographic_model_NH(BR = BR, HR = HR, FAR = FAR, delta_HR = delta_HR, delta_BR = delta_BR, delta_FAR = delta_FAR, r = r10)$value,
    BS_NH    = BS_NH(N = N1, v = v1, relative_frequency = relative_frequency1, r = r11)$value,
    BS_R_NH  = BS_R_NH(N = N8, v = v8, relative_frequency = relative_frequency8, u = u, r = r12)$value,
    MIN_BS_BO = MIN_BS_BO(N = N2, v = v2, relative_frequency = relative_frequency2, BR = BR, p_1 = p1_1, p_2 = p1_2, p_3 = p1_3)$value,
    MIN_BS_HO = MIN_BS_HO(N = N3, v = v3, relative_frequency = relative_frequency3, HR = HR, p_1 = p2_1, p_2 = p2_2, p_3 = p2_3)$value,
    MIN_BS_FO = MIN_BS_FO(N = N4, v = v4, relative_frequency = relative_frequency4, FAR = FAR, p_1 = p3_1, p_2 = p3_2, p_3 = p3_3)$value,
    MIN_BS_JO = MIN_BS_JO(N = N5, v = v5, relative_frequency = relative_frequency5, BR = BR, HR = HR, p_1 = p4_1, p_2 = p4_2, p_3 = p4_3)$value,
    MIN_BS_LS = MIN_BS_LS(N = N6, v = v6, relative_frequency = relative_frequency6, HR = HR, FAR = FAR, p_1 = p5_1, p_2 = p5_2, p_3 = p5_3)$value,
    MIN_BS_H  = MIN_BS_H(N = N7, v = v7, relative_frequency = relative_frequency7, p_1 = p6_1, p_2 = p6_2, p_3 = p6_3)$value,
    M_L = large_mixed_model(p_1 = p8_1, p_2 = p8_2, p_3 = p8_3, 
                            q_1 = q_1, q_2 = q_2, q_3 = q_3, q_4 = q_4, q_5 = q_5, q_6 = q_6,  
                            N = N9, v = v9, relative_frequency = relative_frequency9, 
                            BR = BR, HR = HR, FAR = FAR)$value,
    MH = mixed_heuristic_model(BR = BR, HR = HR, FAR = FAR, 
                               p_1 = p7_1, p_2 = p7_2, p_3 = p7_3, p_4 = p7_4, p_5 = p7_5, p_6 = p7_6, p_7 = p7_7)$value
  ), by = .I]
  pred[, I := NULL]
  return(cbind(dt, pred))
}

pv_full  <- generate_predictions(parameter_dt)
pv_full_s   <- generate_predictions(parameter_dt_s)

tidy_and_save <- function(pv_full, file) {
  out <- as.data.frame(pv_full) %>%
    pivot_longer(cols = B_NH:MH, names_to = "model", values_to = "predict") %>%
    mutate(
      predict = round(as.numeric(predict), 2),
      true_posterior = round(as.numeric(true_posterior), 2)
    )
  saveRDS(out, file)
}

tidy_and_save(pv_full,   "data/prediction_dt.rds")
tidy_and_save(pv_full_s, "data/prediction_dt_s.rds")

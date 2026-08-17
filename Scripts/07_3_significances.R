# =============================================================================
# 7.3 Significances
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# 1 - Read data
# -----------------------------------------------------------------------------

results_combined  <- read.csv("results/all_models_auc_tss.csv")
mae_bias_combined <- read.csv("results/mae_bias_true_suitability.csv")
raw <- read.csv("results/raw_predictions_long.csv")
raw_mae <- read.csv("results/raw_mae_predictions_long.csv")


# -----------------------------------------------------------------------------
# 2 - H1: Predictive performance decreases with increasing niche breadth, 
#         irrespective of modelling algorithm
# -----------------------------------------------------------------------------

library(PMCMRplus)

# TSS-Matrix bauen: Zeilen = Algorithmen (Blöcke), Spalten = Nischenbreite
# absteigend geordnet (breit -> schmal), da wir eine ABNAHME erwarten
tss_matrix <- results_combined %>%
  select(model_type, niche_breadth, TSS) %>%
  tidyr::pivot_wider(names_from = niche_breadth, values_from = TSS) %>%
  select(model_type, broad, high_mid, low_mid, narrow) %>%   # absteigende Reihenfolge!
  tibble::column_to_rownames("model_type") %>%
  as.matrix()

pageTest(tss_matrix, alternative = "greater")


auc_matrix <- results_combined %>%
  select(model_type, niche_breadth, AUC) %>%
  tidyr::pivot_wider(names_from = niche_breadth, values_from = AUC) %>%
  select(model_type, broad, high_mid, low_mid, narrow) %>%   # absteigende Reihenfolge!
  tibble::column_to_rownames("model_type") %>%
  as.matrix()

pageTest(auc_matrix, alternative = "greater")


mae_matrix <- mae_bias_combined %>%
  select(model_type, niche_breadth, MAE) %>%
  tidyr::pivot_wider(names_from = niche_breadth, values_from = MAE) %>%
  select(model_type, narrow, low_mid, high_mid, broad) %>%  # AUFSTEIGEND, da MAE steigen soll
  tibble::column_to_rownames("model_type") %>%
  as.matrix()

pageTest(mae_matrix, alternative = "greater")   # testet: MAE steigt mit Nischenbreite


# -----------------------------------------------------------------------------
# 3 - H2:  Maxent, RF and BRT outperform GLM across all niche breadths
# -----------------------------------------------------------------------------

library(pROC)

#------AUC------
# maxent
glm_narrow    <- raw %>% filter(niche_breadth == "narrow", model_type == "GLM")
maxent_narrow <- raw %>% filter(niche_breadth == "narrow", model_type == "Maxent")

paired <- glm_narrow %>%
  select(x, y, obs, pred_glm = pred) %>%
  inner_join(maxent_narrow %>% select(x, y, pred_maxent = pred), by = c("x", "y"))

roc.test(roc(paired$obs, paired$pred_glm, quiet = TRUE),
         roc(paired$obs, paired$pred_maxent, quiet = TRUE), method = "delong")

# brt
glm_narrow    <- raw %>% filter(niche_breadth == "narrow", model_type == "GLM")
brt_narrow <- raw %>% filter(niche_breadth == "narrow", model_type == "BRT")

paired <- glm_narrow %>%
  select(x, y, obs, pred_glm = pred) %>%
  inner_join(brt_narrow %>% select(x, y, pred_brt = pred), by = c("x", "y"))

roc.test(roc(paired$obs, paired$pred_glm, quiet = TRUE),
         roc(paired$obs, paired$pred_brt, quiet = TRUE), method = "delong")

# rf
glm_narrow    <- raw %>% filter(niche_breadth == "narrow", model_type == "GLM")
rf_narrow <- raw %>% filter(niche_breadth == "narrow", model_type == "RF")

paired <- glm_narrow %>%
  select(x, y, obs, pred_glm = pred) %>%
  inner_join(rf_narrow %>% select(x, y, pred_rf = pred), by = c("x", "y"))

roc.test(roc(paired$obs, paired$pred_glm, quiet = TRUE),
         roc(paired$obs, paired$pred_rf, quiet = TRUE), method = "delong")


#------MAE------
# maxent
glm_narrow    <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "GLM") %>%
  select(x, y, true_suit, pred_glm = pred)
maxent_narrow <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "Maxent") %>%
  select(x, y, pred_maxent = pred)

paired <- glm_narrow %>% inner_join(maxent_narrow, by = c("x", "y"))

err_glm    <- abs(paired$true_suit - paired$pred_glm)
err_maxent <- abs(paired$true_suit - paired$pred_maxent)

wilcox.test(err_maxent, err_glm, paired = TRUE)


# brt
glm_narrow    <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "GLM") %>%
  select(x, y, true_suit, pred_glm = pred)
brt_narrow <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "BRT") %>%
  select(x, y, pred_brt = pred)

paired <- glm_narrow %>% inner_join(brt_narrow, by = c("x", "y"))

err_glm    <- abs(paired$true_suit - paired$pred_glm)
err_brt <- abs(paired$true_suit - paired$pred_brt)

wilcox.test(err_brt, err_glm, paired = TRUE)


# rf
glm_narrow    <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "GLM") %>%
  select(x, y, true_suit, pred_glm = pred)
rf_narrow <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "RF") %>%
  select(x, y, pred_rf = pred)

paired <- glm_narrow %>% inner_join(rf_narrow, by = c("x", "y"))

err_glm    <- abs(paired$true_suit - paired$pred_glm)
err_rf <- abs(paired$true_suit - paired$pred_rf)

wilcox.test(err_rf, err_glm, paired = TRUE)


# -----------------------------------------------------------------------------
# 4 - H3: Maxent, RF and BRT differ in their sensitivity to niche breadth,
#         with Maxent reaching the highest performance
# -----------------------------------------------------------------------------

#------AUC Maxent------
# brt
maxent_narrow    <- raw %>% filter(niche_breadth == "narrow", model_type == "Maxent")
brt_narrow <- raw %>% filter(niche_breadth == "narrow", model_type == "BRT")

paired <- maxent_narrow %>%
  select(x, y, obs, pred_maxent = pred) %>%
  inner_join(brt_narrow %>% select(x, y, pred_brt = pred), by = c("x", "y"))

roc.test(roc(paired$obs, paired$pred_maxent, quiet = TRUE),
         roc(paired$obs, paired$pred_brt, quiet = TRUE), method = "delong")

# rf
maxent_narrow    <- raw %>% filter(niche_breadth == "narrow", model_type == "Maxent")
rf_narrow <- raw %>% filter(niche_breadth == "narrow", model_type == "RF")

paired <- maxent_narrow %>%
  select(x, y, obs, pred_maxent = pred) %>%
  inner_join(rf_narrow %>% select(x, y, pred_rf = pred), by = c("x", "y"))

roc.test(roc(paired$obs, paired$pred_maxent, quiet = TRUE),
         roc(paired$obs, paired$pred_rf, quiet = TRUE), method = "delong")


#------AUC RF------
# brt
rf_narrow    <- raw %>% filter(niche_breadth == "narrow", model_type == "RF")
brt_narrow <- raw %>% filter(niche_breadth == "narrow", model_type == "BRT")

paired <- rf_narrow %>%
  select(x, y, obs, pred_rf = pred) %>%
  inner_join(brt_narrow %>% select(x, y, pred_brt = pred), by = c("x", "y"))

roc.test(roc(paired$obs, paired$pred_rf, quiet = TRUE),
         roc(paired$obs, paired$pred_brt, quiet = TRUE), method = "delong")


#------MAE Maxent------
# brt
maxent_narrow    <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "Maxent") %>%
  select(x, y, true_suit, pred_maxent = pred)
brt_narrow <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "BRT") %>%
  select(x, y, pred_brt = pred)

paired <- maxent_narrow %>% inner_join(brt_narrow, by = c("x", "y"))

err_maxent   <- abs(paired$true_suit - paired$pred_maxent)
err_brt <- abs(paired$true_suit - paired$pred_brt)

wilcox.test(err_brt, err_maxent, paired = TRUE)


# rf
maxent_narrow    <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "Maxent") %>%
  select(x, y, true_suit, pred_maxent = pred)
rf_narrow <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "RF") %>%
  select(x, y, pred_rf = pred)

paired <- maxent_narrow %>% inner_join(rf_narrow, by = c("x", "y"))

err_maxent    <- abs(paired$true_suit - paired$pred_maxent)
err_rf <- abs(paired$true_suit - paired$pred_rf)

wilcox.test(err_rf, err_maxent, paired = TRUE)


#------MAE RF------
# brt
rf_narrow    <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "RF") %>%
  select(x, y, true_suit, pred_rf = pred)
brt_narrow <- raw_mae %>% filter(niche_breadth == "narrow", model_type == "BRT") %>%
  select(x, y, pred_brt = pred)

paired <- rf_narrow %>% inner_join(brt_narrow, by = c("x", "y"))

err_rf   <- abs(paired$true_suit - paired$pred_rf)
err_brt <- abs(paired$true_suit - paired$pred_brt)

wilcox.test(err_brt, err_rf, paired = TRUE)

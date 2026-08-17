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



library(pROC)

paired_bootstrap_tss_diff <- function(obs, pred_A, pred_B, n_boot = 1000) {
  n <- length(obs)
  diffs <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    if (length(unique(obs[idx])) < 2) next
    tss_A <- evalSDM(obs[idx], pred_A[idx], thresh.method = "ObsPrev")$TSS
    tss_B <- evalSDM(obs[idx], pred_B[idx], thresh.method = "ObsPrev")$TSS
    diffs[i] <- tss_A - tss_B
  }
  p_approx <- mean(diffs <= 0, na.rm = TRUE)
  list(diff_mean = mean(diffs, na.rm = TRUE),
       p = 2 * min(p_approx, 1 - p_approx, na.rm = TRUE))
}

h2_results <- list()

for (nm in niche_names) {
  glm_dat <- raw %>% filter(niche_breadth == nm, model_type == "GLM") %>%
    select(x, y, obs, pred_glm = pred)
  
  for (mtype in c("Maxent", "BRT", "RF")) {
    
    other_dat <- raw %>% filter(niche_breadth == nm, model_type == mtype) %>%
      select(x, y, pred_other = pred)
    
    paired <- glm_dat %>% inner_join(other_dat, by = c("x", "y"))
    
    # --- AUC: DeLong ---
    roc_glm   <- roc(paired$obs, paired$pred_glm,   quiet = TRUE)
    roc_other <- roc(paired$obs, paired$pred_other, quiet = TRUE)
    auc_test  <- roc.test(roc_glm, roc_other, method = "delong")
    
    # --- TSS: gepaarter Bootstrap ---
    tss_test <- paired_bootstrap_tss_diff(paired$obs, paired$pred_other, paired$pred_glm)
    
    h2_results[[paste(nm, mtype, sep = "_")]] <- data.frame(
      niche_breadth = nm,
      comparison    = paste(mtype, "vs GLM"),
      AUC_glm       = as.numeric(auc(roc_glm)),
      AUC_other     = as.numeric(auc(roc_other)),
      p_AUC         = auc_test$p.value,
      TSS_diff      = tss_test$diff_mean,
      p_TSS         = tss_test$p
    )
  }
}

h2_df <- dplyr::bind_rows(h2_results)
h2_df$p_AUC_adj <- p.adjust(h2_df$p_AUC, method = "holm")
h2_df$p_TSS_adj <- p.adjust(h2_df$p_TSS, method = "holm")

print(h2_df)


# -----------------------------------------------------------------------------
# 4 - H3: Maxent, RF and BRT differ in their sensitivity to niche breadth,
#         with Maxent reaching the highest performance
# -----------------------------------------------------------------------------


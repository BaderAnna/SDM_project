# =============================================================================
# 7.3 Significances
# =============================================================================

library(tidyverse)
library(PMCMRplus)
library(pROC)

# -----------------------------------------------------------------------------
# 1 - Read data
# -----------------------------------------------------------------------------

results_combined  <- read.csv("results/all_models_auc_tss.csv")
mae_bias_combined <- read.csv("results/mae_bias_true_suitability.csv")
raw               <- read.csv("results/raw_predictions_long.csv")
raw_mae           <- read.csv("results/raw_mae_predictions_long.csv")

# Nischenbreite als geordneten Faktor definieren
niche_order <- c("narrow", "low_mid", "high_mid", "broad")

results_combined <- results_combined %>%
  mutate(niche_breadth = factor(niche_breadth, levels = niche_order))

mae_bias_combined <- mae_bias_combined %>%
  mutate(niche_breadth = factor(niche_breadth, levels = niche_order))

# Sammelliste für alle Ergebnisse
all_results <- list()

# -----------------------------------------------------------------------------
# 2 - H1: Predictive performance decreases with increasing niche breadth,
#         irrespective of modelling algorithm
# -----------------------------------------------------------------------------

# --- TSS ---
tss_matrix <- results_combined %>%
  select(model_type, niche_breadth, TSS) %>%
  pivot_wider(names_from = niche_breadth, values_from = TSS) %>%
  select(model_type, broad, high_mid, low_mid, narrow) %>%  # absteigend
  column_to_rownames("model_type") %>%
  as.matrix()

tss_page <- pageTest(tss_matrix, alternative = "greater")

all_results[["H1_TSS_Page"]] <- tibble(
  Hypothesis  = "H1",
  Metric      = "TSS",
  Test        = "Page Test",
  Statistic   = round(tss_page$statistic, 3),
  p_value     = round(tss_page$p.value, 4),
  Significant = ifelse(tss_page$p.value < 0.05, "yes", "no")
)

# --- AUC ---
auc_matrix <- results_combined %>%
  select(model_type, niche_breadth, AUC) %>%
  pivot_wider(names_from = niche_breadth, values_from = AUC) %>%
  select(model_type, broad, high_mid, low_mid, narrow) %>%  # absteigend
  column_to_rownames("model_type") %>%
  as.matrix()

auc_page <- pageTest(auc_matrix, alternative = "greater")

all_results[["H1_AUC_Page"]] <- tibble(
  Hypothesis  = "H1",
  Metric      = "AUC",
  Test        = "Page Test",
  Statistic   = round(auc_page$statistic, 3),
  p_value     = round(auc_page$p.value, 4),
  Significant = ifelse(auc_page$p.value < 0.05, "yes", "no")
)

# --- MAE ---
mae_matrix <- mae_bias_combined %>%
  select(model_type, niche_breadth, MAE) %>%
  pivot_wider(names_from = niche_breadth, values_from = MAE) %>%
  select(model_type, narrow, low_mid, high_mid, broad) %>%  # aufsteigend
  column_to_rownames("model_type") %>%
  as.matrix()

mae_page <- pageTest(mae_matrix, alternative = "greater")

all_results[["H1_MAE_Page"]] <- tibble(
  Hypothesis  = "H1",
  Metric      = "MAE",
  Test        = "Page Test",
  Statistic   = round(mae_page$statistic, 3),
  p_value     = round(mae_page$p.value, 4),
  Significant = ifelse(mae_page$p.value < 0.05, "yes", "no")
)

# -----------------------------------------------------------------------------
# 3 - H2: Maxent, RF and BRT outperform GLM across all niche breadths
# -----------------------------------------------------------------------------

# Hilfsfunktion: AUC-Vergleich über alle Nischenbreiten
run_roc_test <- function(data, nb, m1, m2, pred1_name, pred2_name) {
  d1 <- data %>% filter(niche_breadth == nb, model_type == m1)
  d2 <- data %>% filter(niche_breadth == nb, model_type == m2)
  
  joined <- d1 %>%
    select(x, y, obs, pred1 = pred) %>%
    inner_join(d2 %>% select(x, y, pred2 = pred), by = c("x", "y"))
  
  test <- roc.test(
    roc(joined$obs, joined$pred1, quiet = TRUE),
    roc(joined$obs, joined$pred2, quiet = TRUE),
    method = "delong"
  )
  
  tibble(
    niche_breadth = nb,
    model1        = m1,
    model2        = m2,
    AUC_model1    = round(test$roc1$auc, 3),
    AUC_model2    = round(test$roc2$auc, 3),
    Statistic     = round(test$statistic, 3),
    p_value       = round(test$p.value, 4),
    Significant   = ifelse(test$p.value < 0.05, "yes", "no")
  )
}

# Hilfsfunktion: MAE-Vergleich über alle Nischenbreiten
run_mae_test <- function(data, nb, m1, m2, pred1_name, pred2_name) {
  d1 <- data %>% filter(niche_breadth == nb, model_type == m1) %>%
    select(x, y, true_suit, pred1 = pred)
  d2 <- data %>% filter(niche_breadth == nb, model_type == m2) %>%
    select(x, y, pred2 = pred)
  
  joined <- d1 %>% inner_join(d2, by = c("x", "y"))
  
  err1 <- abs(joined$true_suit - joined$pred1)
  err2 <- abs(joined$true_suit - joined$pred2)
  
  test <- wilcox.test(err1, err2, paired = TRUE)
  
  tibble(
    niche_breadth = nb,
    model1        = m1,
    model2        = m2,
    MAE_model1    = round(mean(err1), 4),
    MAE_model2    = round(mean(err2), 4),
    Statistic     = round(test$statistic, 3),
    p_value       = round(test$p.value, 4),
    Significant   = ifelse(test$p.value < 0.05, "yes", "no")
  )
}

# --- H2: AUC - GLM vs. Maxent, BRT, RF über alle Nischenbreiten ---
h2_comparisons <- list(
  c("GLM", "Maxent"),
  c("GLM", "BRT"),
  c("GLM", "RF")
)

h2_auc_results <- map_dfr(niche_order, function(nb) {
  map_dfr(h2_comparisons, function(pair) {
    run_roc_test(raw, nb, pair[1], pair[2])
  })
}) %>%
  mutate(Hypothesis = "H2", Metric = "AUC", Test = "DeLong") %>%
  select(Hypothesis, Metric, Test, niche_breadth, model1, model2,
         AUC_model1, AUC_model2, Statistic, p_value, Significant)

all_results[["H2_AUC"]] <- h2_auc_results

# --- H2: MAE - GLM vs. Maxent, BRT, RF über alle Nischenbreiten ---
h2_mae_results <- map_dfr(niche_order, function(nb) {
  map_dfr(h2_comparisons, function(pair) {
    run_mae_test(raw_mae, nb, pair[1], pair[2])
  })
}) %>%
  mutate(Hypothesis = "H2", Metric = "MAE", Test = "Wilcoxon") %>%
  select(Hypothesis, Metric, Test, niche_breadth, model1, model2,
         MAE_model1, MAE_model2, Statistic, p_value, Significant)

all_results[["H2_MAE"]] <- h2_mae_results

# -----------------------------------------------------------------------------
# 4 - H3: Models differ in sensitivity to niche breadth;
#         Maxent reaches highest performance
# -----------------------------------------------------------------------------

# ~~~ Teil 1: Unterschiedliche Sensitivität = Interaktion Modell × Nischenbreite ~~~
#
# Logik: Wenn Modelle sich unterschiedlich stark mit der Nischenbreite verändern,
# zeigt sich das als signifikante Interaktion in einem gemischten Modell.
# Da wir keine Replikation haben (1 Wert pro Modell × Nischenbreite),
# verwenden wir den Scheirer-Ray-Hare Test (nichtparametrische 2-Faktor ANOVA)

library(rcompanion)

# --- H3 Teil 1: AUC - Interaktion Modell × Nischenbreite ---
auc_long <- results_combined %>%
  mutate(
    model_type    = factor(model_type),
    niche_breadth = factor(niche_breadth, levels = niche_order)
  ) %>%
  filter(model_type != "GLM")  # nur die 3 Modelle aus H3

srh_auc <- scheirerRayHare(AUC ~ model_type + niche_breadth,
                           data = auc_long)

# Interaktionszeile extrahieren
srh_auc_interaction <- srh_auc %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(grepl(":", Term)) %>%
  transmute(
    Hypothesis  = "H3",
    Metric      = "AUC",
    Test        = "Scheirer-Ray-Hare (Interaktion)",
    Term        = Term,
    Df          = Df,
    Statistic   = round(H, 3),
    p_value     = round(p.value, 4),
    Significant = ifelse(p.value < 0.05, "yes", "no")
  )

all_results[["H3_Interaction_AUC"]] <- srh_auc_interaction

# --- H3 Teil 1: TSS - Interaktion Modell × Nischenbreite ---
tss_long <- results_combined %>%
  mutate(
    model_type    = factor(model_type),
    niche_breadth = factor(niche_breadth, levels = niche_order)
  ) %>%
  filter(model_type != "GLM")

srh_tss <- scheirerRayHare(TSS ~ model_type + niche_breadth,
                           data = tss_long)

srh_tss_interaction <- srh_tss %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(grepl(":", Term)) %>%
  transmute(
    Hypothesis  = "H3",
    Metric      = "TSS",
    Test        = "Scheirer-Ray-Hare (Interaktion)",
    Term        = Term,
    Df          = Df,
    Statistic   = round(H, 3),
    p_value     = round(p.value, 4),
    Significant = ifelse(p.value < 0.05, "yes", "no")
  )

all_results[["H3_Interaction_TSS"]] <- srh_tss_interaction

# --- H3 Teil 1: MAE - Interaktion Modell × Nischenbreite ---
mae_long <- mae_bias_combined %>%
  mutate(
    model_type    = factor(model_type),
    niche_breadth = factor(niche_breadth, levels = niche_order)
  ) %>%
  filter(model_type != "GLM")

srh_mae <- scheirerRayHare(MAE ~ model_type + niche_breadth,
                           data = mae_long)

srh_mae_interaction <- srh_mae %>%
  as.data.frame() %>%
  rownames_to_column("Term") %>%
  filter(grepl(":", Term)) %>%
  transmute(
    Hypothesis  = "H3",
    Metric      = "MAE",
    Test        = "Scheirer-Ray-Hare (Interaktion)",
    Term        = Term,
    Df          = Df,
    Statistic   = round(H, 3),
    p_value     = round(p.value, 4),
    Significant = ifelse(p.value < 0.05, "yes", "no")
  )

all_results[["H3_Interaction_MAE"]] <- srh_mae_interaction

# ~~~ Teil 2: MaxEnt hat höchste Performance ~~~
# Vergleich MaxEnt vs. BRT und MaxEnt vs. RF über alle Nischenbreiten

h3_comparisons <- list(
  c("Maxent", "BRT"),
  c("Maxent", "RF"),
  c("RF", "BRT")   # vollständiger paarweiser Vergleich
)

# --- H3 Teil 2: AUC ---
h3_auc_results <- map_dfr(niche_order, function(nb) {
  map_dfr(h3_comparisons, function(pair) {
    run_roc_test(raw, nb, pair[1], pair[2])
  })
}) %>%
  mutate(Hypothesis = "H3", Metric = "AUC", Test = "DeLong") %>%
  select(Hypothesis, Metric, Test, niche_breadth, model1, model2,
         AUC_model1, AUC_model2, Statistic, p_value, Significant)

all_results[["H3_AUC_MaxEnt_highest"]] <- h3_auc_results

# --- H3 Teil 2: MAE ---
h3_mae_results <- map_dfr(niche_order, function(nb) {
  map_dfr(h3_comparisons, function(pair) {
    run_mae_test(raw_mae, nb, pair[1], pair[2])
  })
}) %>%
  mutate(Hypothesis = "H3", Metric = "MAE", Test = "Wilcoxon") %>%
  select(Hypothesis, Metric, Test, niche_breadth, model1, model2,
         MAE_model1, MAE_model2, Statistic, p_value, Significant)

all_results[["H3_MAE_MaxEnt_highest"]] <- h3_mae_results

# -----------------------------------------------------------------------------
# 5 - Ergebnistabelle zusammenführen und speichern
# -----------------------------------------------------------------------------

# H1-Ergebnisse in kompatibles Format bringen
h1_table <- bind_rows(
  all_results[["H1_TSS_Page"]],
  all_results[["H1_AUC_Page"]],
  all_results[["H1_MAE_Page"]]
) %>%
  mutate(
    niche_breadth = "all",
    model1        = "all models",
    model2        = NA_character_,
    AUC_model1    = NA_real_,
    AUC_model2    = NA_real_,
    MAE_model1    = NA_real_,
    MAE_model2    = NA_real_,
    Term          = NA_character_,
    Df            = NA_real_
  )

# H3 Interaktions-Ergebnisse in kompatibles Format bringen
h3_int_table <- bind_rows(
  all_results[["H3_Interaction_AUC"]],
  all_results[["H3_Interaction_TSS"]],
  all_results[["H3_Interaction_MAE"]]
) %>%
  mutate(
    niche_breadth = "all",
    model1        = "Maxent/RF/BRT",
    model2        = NA_character_,
    AUC_model1    = NA_real_,
    AUC_model2    = NA_real_,
    MAE_model1    = NA_real_,
    MAE_model2    = NA_real_
  )

# Alle Tabellen zusammenführen
final_table <- bind_rows(
  h1_table,
  all_results[["H2_AUC"]],
  all_results[["H2_MAE"]],
  h3_int_table,
  all_results[["H3_AUC_MaxEnt_highest"]],
  all_results[["H3_MAE_MaxEnt_highest"]]
) %>%
  select(Hypothesis, Metric, Test, niche_breadth, model1, model2,
         Statistic, p_value, Significant) %>%
  arrange(Hypothesis, Metric, niche_breadth)

# Speichern
write.csv(final_table, "results/significance_results.csv", row.names = FALSE)

# Übersichtliche Ausgabe in der Konsole
cat("\n========================================\n")
cat("    SIGNIFICANCE TEST RESULTS SUMMARY\n")
cat("========================================\n\n")

cat("--- H1: Performance decreases with niche breadth ---\n")
print(bind_rows(all_results[["H1_TSS_Page"]],
                all_results[["H1_AUC_Page"]],
                all_results[["H1_MAE_Page"]]))

cat("\n--- H2: Maxent/BRT/RF outperform GLM (AUC) ---\n")
print(all_results[["H2_AUC"]])

cat("\n--- H2: Maxent/BRT/RF outperform GLM (MAE) ---\n")
print(all_results[["H2_MAE"]])

cat("\n--- H3 Teil 1: Interaction Model x Niche Breadth ---\n")
print(bind_rows(
  all_results[["H3_Interaction_AUC"]],
  all_results[["H3_Interaction_TSS"]],
  all_results[["H3_Interaction_MAE"]]
))

cat("\n--- H3 Teil 2: MaxEnt highest performance (AUC) ---\n")
print(all_results[["H3_AUC_MaxEnt_highest"]])

cat("\n--- H3 Teil 2: MaxEnt lowest MAE ---\n")
print(all_results[["H3_MAE_MaxEnt_highest"]])
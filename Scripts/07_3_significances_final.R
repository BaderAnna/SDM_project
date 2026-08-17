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

niche_order <- c("narrow", "low_mid", "high_mid", "broad")

results_combined <- results_combined %>%
  mutate(niche_breadth = factor(niche_breadth, levels = niche_order))

mae_bias_combined <- mae_bias_combined %>%
  mutate(niche_breadth = factor(niche_breadth, levels = niche_order))

all_results <- list()

# -----------------------------------------------------------------------------
# 2 - Hilfsfunktionen
# -----------------------------------------------------------------------------

format_p <- function(p) {
  ifelse(p < 0.001, "< 0.001", as.character(round(p, 3)))
}

run_roc_test <- function(data, nb, m1, m2) {
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
    Statistic     = round(test$statistic, 3),
    p_value       = format_p(test$p.value), 
    Significant   = ifelse(test$p.value < 0.05, "yes", "no")
  )
}

run_mae_test <- function(data, nb, m1, m2) {
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
    Statistic     = round(test$statistic, 3),
    p_value       = format_p(test$p.value), 
    Significant   = ifelse(test$p.value < 0.05, "yes", "no")
  )
}

# -----------------------------------------------------------------------------
# 3 - H1: Predictive performance decreases with increasing niche breadth,
#         irrespective of modelling algorithm
# -----------------------------------------------------------------------------

tss_matrix <- results_combined %>%
  select(model_type, niche_breadth, TSS) %>%
  pivot_wider(names_from = niche_breadth, values_from = TSS) %>%
  select(model_type, broad, high_mid, low_mid, narrow) %>%
  column_to_rownames("model_type") %>%
  as.matrix()

tss_page <- pageTest(tss_matrix, alternative = "greater")

all_results[["H1_TSS"]] <- tibble(
  Hypothesis  = "H1",
  Metric      = "TSS",
  Test        = "Page Test",
  Comparison  = "all models",
  niche_breadth = "all",
  Statistic   = round(tss_page$statistic, 3),
  p_value     = format_p(tss_page$p.value),
  Significant = ifelse(tss_page$p.value < 0.05, "yes", "no")
)

auc_matrix <- results_combined %>%
  select(model_type, niche_breadth, AUC) %>%
  pivot_wider(names_from = niche_breadth, values_from = AUC) %>%
  select(model_type, broad, high_mid, low_mid, narrow) %>%
  column_to_rownames("model_type") %>%
  as.matrix()

auc_page <- pageTest(auc_matrix, alternative = "greater")

all_results[["H1_AUC"]] <- tibble(
  Hypothesis    = "H1",
  Metric        = "AUC",
  Test          = "Page Test",
  Comparison    = "all models",
  niche_breadth = "all",
  Statistic     = round(auc_page$statistic, 3),
  p_value       = format_p(auc_page$p.value),
  Significant   = ifelse(auc_page$p.value < 0.05, "yes", "no")
)

mae_matrix <- mae_bias_combined %>%
  select(model_type, niche_breadth, MAE) %>%
  pivot_wider(names_from = niche_breadth, values_from = MAE) %>%
  select(model_type, narrow, low_mid, high_mid, broad) %>%
  column_to_rownames("model_type") %>%
  as.matrix()

mae_page <- pageTest(mae_matrix, alternative = "greater")

all_results[["H1_MAE"]] <- tibble(
  Hypothesis    = "H1",
  Metric        = "MAE",
  Test          = "Page Test",
  Comparison    = "all models",
  niche_breadth = "all",
  Statistic     = round(mae_page$statistic, 3),
  p_value       = format_p(mae_page$p.value),
  Significant   = ifelse(mae_page$p.value < 0.05, "yes", "no")
)

# -----------------------------------------------------------------------------
# 4 - H2: Maxent, RF and BRT outperform GLM across all niche breadths
# -----------------------------------------------------------------------------

h2_comparisons <- list(
  c("GLM", "Maxent"),
  c("GLM", "BRT"),
  c("GLM", "RF")
)

# Sortierung: erst nach Vergleichspaar, dann narrow -> broad
h2_auc_results <- map_dfr(h2_comparisons, function(pair) {
  map_dfr(niche_order, function(nb) {
    run_roc_test(raw, nb, pair[1], pair[2])
  })
}) %>%
  mutate(
    Hypothesis  = "H2",
    Metric      = "AUC",
    Test        = "DeLong",
    Comparison  = paste(model1, "vs", model2)
  ) %>%
  select(Hypothesis, Metric, Test, Comparison, niche_breadth, Statistic, p_value, Significant)

all_results[["H2_AUC"]] <- h2_auc_results

h2_mae_results <- map_dfr(h2_comparisons, function(pair) {
  map_dfr(niche_order, function(nb) {
    run_mae_test(raw_mae, nb, pair[1], pair[2])
  })
}) %>%
  mutate(
    Hypothesis  = "H2",
    Metric      = "MAE",
    Test        = "Wilcoxon",
    Comparison  = paste(model1, "vs", model2)
  ) %>%
  select(Hypothesis, Metric, Test, Comparison, niche_breadth, Statistic, p_value, Significant)

all_results[["H2_MAE"]] <- h2_mae_results

# -----------------------------------------------------------------------------
# 5 - H3: Models differ in sensitivity to niche breadth;
#         Maxent reaches highest performance
#
# H3 Teil 1 (sensitivity): Verändern sich die Modelle unterschiedlich stark?
#   → Ablesbar daran, ob Signifikanz/Differenz sich über Nischenbreiten verändert
#   → Interpretation: p-Werte narrow vs. broad vergleichen
#
# H3 Teil 2 (highest performance): Hat Maxent generell die höchsten Werte?
#   → Ablesbar daran, ob Maxent in den Vergleichen konsistent gewinnt
#
# Beide Fragen werden durch dieselben paarweisen Tests beantwortet.
# Sortierung: erst Vergleichspaar, dann narrow -> broad
# -----------------------------------------------------------------------------

h3_comparisons <- list(
  c("Maxent", "BRT"),
  c("Maxent", "RF"),
  c("RF",     "BRT")
)

# --- H3: AUC ---
h3_auc_results <- map_dfr(h3_comparisons, function(pair) {
  map_dfr(niche_order, function(nb) {
    run_roc_test(raw, nb, pair[1], pair[2])
  })
}) %>%
  mutate(
    Hypothesis = "H3",
    Metric     = "AUC",
    Test       = "DeLong",
    Comparison = paste(model1, "vs", model2)
  ) %>%
  select(Hypothesis, Metric, Test, Comparison, niche_breadth, Statistic, p_value, Significant)

all_results[["H3_AUC"]] <- h3_auc_results

# --- H3: MAE ---
h3_mae_results <- map_dfr(h3_comparisons, function(pair) {
  map_dfr(niche_order, function(nb) {
    run_mae_test(raw_mae, nb, pair[1], pair[2])
  })
}) %>%
  mutate(
    Hypothesis = "H3",
    Metric     = "MAE",
    Test       = "Wilcoxon",
    Comparison = paste(model1, "vs", model2)
  ) %>%
  select(Hypothesis, Metric, Test, Comparison, niche_breadth, Statistic, p_value, Significant)

all_results[["H3_MAE"]] <- h3_mae_results

# -----------------------------------------------------------------------------
# 6 - Ergebnistabelle zusammenführen und speichern
# -----------------------------------------------------------------------------

final_table <- bind_rows(
  # H1: Page Tests
  all_results[["H1_AUC"]],
  all_results[["H1_TSS"]],
  all_results[["H1_MAE"]],
  # H2: GLM vs. andere
  all_results[["H2_AUC"]],
  all_results[["H2_MAE"]],
  # H3: Maxent/BRT/RF paarweise
  all_results[["H3_AUC"]],
  all_results[["H3_MAE"]]
)

write.csv(final_table, "results/significance_results.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 7 - Konsolenausgabe
# -----------------------------------------------------------------------------

cat("\n========================================\n")
cat("    SIGNIFICANCE TEST RESULTS SUMMARY\n")
cat("========================================\n")

cat("\n--- H1: Performance ~ niche breadth (Page Test) ---\n")
all_results[["H1_AUC"]] %>% print()
all_results[["H1_TSS"]] %>% print()
all_results[["H1_MAE"]] %>% print()

cat("\n--- H2: GLM vs. Maxent/BRT/RF (AUC - DeLong) ---\n")
all_results[["H2_AUC"]] %>% print(n = Inf)

cat("\n--- H2: GLM vs. Maxent/BRT/RF (MAE - Wilcoxon) ---\n")
all_results[["H2_MAE"]] %>% print(n = Inf)

cat("\n--- H3: Pairwise comparisons (AUC - DeLong) ---\n")
cat("    Teil 1: Sensitivity = changes in significance across niche breadths\n")
cat("    Teil 2: Maxent performance = Maxent wins consistently\n")
all_results[["H3_AUC"]] %>% print(n = Inf)

cat("\n--- H3: Pairwise comparisons (MAE - Wilcoxon) ---\n")
cat("    Teil 1: Sensitivity = changes in significance across niche breadths\n")
cat("    Teil 2: Maxent performance = Maxent wins consistently\n")
all_results[["H3_MAE"]] %>% print(n = Inf)

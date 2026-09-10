# =============================================================================
# MAE und Bias 
# =============================================================================

library(Metrics)
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# 1 - Hilfsfunktion: Vorhersagen für einen Modelltyp + Nischenbreite holen
#     (identisch zur predict_model()-Logik aus dem Evaluierungsskript)
# -----------------------------------------------------------------------------

predict_model <- function(model_obj, newdata, model_type) {
  
  switch(model_type,
         
         "GLM" = predict(model_obj$model,
                         newdata = newdata,
                         type    = "response"),
         
         "GAM" = predict(model_obj$model,
                         newdata = newdata,
                         type    = "response"),
         
         "BRT" = {
           gbm_fit <- model_obj$model
           gbm::predict.gbm(gbm_fit,
                            newdata = newdata,
                            n.trees = gbm_fit$n.trees,
                            type    = "response")
         },
         
         "RF"  = {
           rf_fit <- model_obj$model
           as.numeric(
             predict(rf_fit, newdata = newdata, type = "prob")[, "1"]
           )
         },
         
         # Fallback
         predict(model_obj, newdata = newdata, type = "response")
  )
}

# -----------------------------------------------------------------------------
# 2 - Kernfunktion: MAE + Bias für eine Modell-Nischen-Kombination
# -----------------------------------------------------------------------------
# Metrics::mae()  = mean(|obs - pred|)   → mittlerer absoluter Fehler
# Metrics::bias() = mean(actual - predicted) = mean(obs - pred)
#                   positive = model underestimates (obs > pred)
#                   negative = model overestimates (obs < pred)
# -----------------------------------------------------------------------------

calc_mae_bias <- function(model_obj, test_df, env,
                          model_type    = "GLM",
                          n_boot        = 1000,
                          boot_seed     = 42) {
  
  # Nur Testzeilen
  test_df  <- test_df[test_df$split == "test", ]
  
  # Umweltwerte extrahieren
  env_vals <- terra::extract(env, test_df[, c("x", "y")])[, -1, drop = FALSE]
  newdat   <- cbind(test_df, env_vals)
  
  # Vorhersage
  pred <- tryCatch(
    predict_model(model_obj, newdat, model_type),
    error = function(e) {
      warning(sprintf("[%s] predict() fehlgeschlagen: %s", model_type, e$message))
      rep(NA_real_, nrow(newdat))
    }
  )
  
  # NA-Zeilen entfernen
  ok   <- !is.na(pred)
  obs  <- newdat$presence[ok]
  pred <- pred[ok]
  
  if (sum(ok) == 0) {
    warning(sprintf("[%s] Keine gültigen Vorhersagen.", model_type))
    return(data.frame(MAE = NA, Bias = NA,
                      MAE_lo = NA, MAE_hi = NA,
                      Bias_lo = NA, Bias_hi = NA,
                      n_test = 0L))
  }
  
  # --- Punktschätzer ---------------------------------------------------------
  mae_val  <- Metrics::mae(actual    = obs, predicted = pred)
  bias_val <- Metrics::bias(actual   = obs, predicted = pred)
  
  # --- Bootstrap-KI ----------------------------------------------------------
  set.seed(boot_seed)
  n        <- length(obs)
  mae_boot <- numeric(n_boot)
  bias_boot <- numeric(n_boot)
  
  for (i in seq_len(n_boot)) {
    idx          <- sample(seq_len(n), size = n, replace = TRUE)
    mae_boot[i]  <- Metrics::mae( actual = obs[idx], predicted = pred[idx])
    bias_boot[i] <- Metrics::bias(actual = obs[idx], predicted = pred[idx])
  }
  
  mae_ci  <- quantile(mae_boot,  probs = c(0.025, 0.975), na.rm = TRUE)
  bias_ci <- quantile(bias_boot, probs = c(0.025, 0.975), na.rm = TRUE)
  
  data.frame(
    MAE     = mae_val,
    Bias    = bias_val,
    MAE_lo  = mae_ci["2.5%"],
    MAE_hi  = mae_ci["97.5%"],
    Bias_lo = bias_ci["2.5%"],
    Bias_hi = bias_ci["97.5%"],
    n_test  = sum(ok),
    row.names = NULL
  )
}

# -----------------------------------------------------------------------------
# 3 - Hauptschleife: alle Modelltypen × alle Nischenbreiten
# -----------------------------------------------------------------------------

niche_names <- c("narrow", "low_mid", "high_mid", "broad")

# Modelle und Testdaten (wie im Evaluierungsskript geladen)
# models_all und test_data werden als bereits vorhanden vorausgesetzt;
# falls dieses Skript standalone läuft, hier erneut laden:
#
# env       <- terra::rast("Data/raster/env_raster_masked.tif")
# models_all <- list(GLM = models_glm, GAM = models_gam, BRT = models_brt)
# test_data  <- list(narrow = ..., low_mid = ..., ...)

metrics_list <- list()

for (mtype in names(models_all)) {
  
  cat("\n==============================\n")
  cat("Berechne MAE/Bias für:", mtype, "\n")
  cat("==============================\n")
  
  res_list <- lapply(niche_names, function(nm) {
    
    cat("  Nischenbreite:", nm, "\n")
    
    calc_mae_bias(
      model_obj  = models_all[[mtype]][[nm]],
      test_df    = test_data[[nm]],
      env        = env,
      model_type = mtype,
      n_boot     = 1000
    )
  })
  
  res_df <- do.call(rbind, res_list)
  res_df$niche_breadth <- niche_names
  res_df$model_type    <- mtype
  rownames(res_df)     <- NULL
  
  metrics_list[[mtype]] <- res_df
}

# Alle Ergebnisse zusammenführen
metrics_combined <- do.call(rbind, metrics_list) %>%
  mutate(
    niche_breadth = factor(niche_breadth, levels = niche_names),
    model_type    = factor(model_type,    levels = names(models_all))
  )

print(metrics_combined)

# -----------------------------------------------------------------------------
# 4 - Speichern
# -----------------------------------------------------------------------------

write.csv(metrics_combined, "results/all_models_mae_bias.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 5 - Plots
# -----------------------------------------------------------------------------

n_models <- length(unique(metrics_combined$model_type))
farben   <- RColorBrewer::brewer.pal(max(3, n_models), "Set1")[seq_len(n_models)]

# --- 5a: MAE mit Bootstrap-KI -----------------------------------------------
p_mae <- ggplot(metrics_combined,
                aes(x     = niche_breadth,
                    y     = MAE,
                    color = model_type,
                    group = model_type)) +
  
  geom_ribbon(aes(ymin = MAE_lo,
                  ymax = MAE_hi,
                  fill = model_type),
              alpha = 0.15,
              color = NA) +
  
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  
  # MAE = 0 wäre perfekt
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             color      = "grey50") +
  
  scale_color_manual(values = farben) +
  scale_fill_manual(values  = farben) +
  
  labs(
    title    = "MAE: Modelltypen × Nischenbreite",
    subtitle = "Punkte = MAE; Bänder = 95%-Bootstrap-KI",
    x        = "Nischenbreite",
    y        = "MAE",
    color    = "Modelltyp",
    fill     = "Modelltyp"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

print(p_mae)

# --- 5b: Bias mit Bootstrap-KI ----------------------------------------------
p_bias <- ggplot(metrics_combined,
                 aes(x     = niche_breadth,
                     y     = Bias,
                     color = model_type,
                     group = model_type)) +
  
  geom_ribbon(aes(ymin = Bias_lo,
                  ymax = Bias_hi,
                  fill = model_type),
              alpha = 0.15,
              color = NA) +
  
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  
  # Bias = 0 bedeutet keine systematische Verzerrung
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             color      = "grey50") +
  
  # Beschriftung der Nulllinie
  annotate("text",
           x     = 0.6,           # links vom ersten Punkt
           y     = 0.003,         # knapp über der Nulllinie
           label = "kein Bias",
           color = "grey40",
           size  = 3.5,
           hjust = 0) +
  
  scale_color_manual(values = farben) +
  scale_fill_manual(values  = farben) +
  
  labs(
    title    = "Bias: Modelltypen × Nischenbreite",
    subtitle = "positiv = Überschätzung | negativ = Unterschätzung",
    x        = "Nischenbreite",
    y        = "Bias  (mean predicted − mean observed)",
    color    = "Modelltyp",
    fill     = "Modelltyp"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

print(p_bias)

# --- 5c: MAE + Bias kombiniert (Facets) ------------------------------------
metrics_long <- metrics_combined %>%
  select(model_type, niche_breadth, MAE, Bias, MAE_lo, MAE_hi, Bias_lo, Bias_hi) %>%
  pivot_longer(
    cols      = c(MAE, Bias),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    # KI-Grenzen je nach Metrik zuordnen
    ci_lo = if_else(metric == "MAE", MAE_lo, Bias_lo),
    ci_hi = if_else(metric == "MAE", MAE_hi, Bias_hi)
  )

p_combined <- ggplot(metrics_long,
                     aes(x     = niche_breadth,
                         y     = value,
                         color = model_type,
                         group = model_type)) +
  
  geom_ribbon(aes(ymin = ci_lo,
                  ymax = ci_hi,
                  fill = model_type),
              alpha = 0.15,
              color = NA) +
  
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             color      = "grey50") +
  
  facet_wrap(~ metric,
             scales = "free_y",   # y-Achsen unabhängig skalieren
             nrow   = 1) +
  
  scale_color_manual(values = farben) +
  scale_fill_manual(values  = farben) +
  
  labs(
    title    = "MAE und Bias: Modellvergleich",
    subtitle = "Bänder = 95%-Bootstrap-KI",
    x        = "Nischenbreite",
    y        = "Wert",
    color    = "Modelltyp",
    fill     = "Modelltyp"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

print(p_combined)

# --- Plots speichern --------------------------------------------------------
ggsave("results/plot_mae.png",      p_mae,      width = 8,  height = 5, dpi = 300)
ggsave("results/plot_bias.png",     p_bias,     width = 8,  height = 5, dpi = 300)
ggsave("results/plot_mae_bias.png", p_combined, width = 12, height = 5, dpi = 300)
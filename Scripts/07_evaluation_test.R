# =============================================================================
# TSS und AUC - Erweiterter Vergleich mehrerer Modelltypen
# =============================================================================

library(mecofun)
library(PresenceAbsence)
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(maxnet)
library(ranger)

# -----------------------------------------------------------------------------
# 1 - Umweltraster laden
# -----------------------------------------------------------------------------
env <- terra::rast("Data/raster/env_raster_masked.tif")
print(names(env))

# -----------------------------------------------------------------------------
# 2 - Alle Modelle und Testdaten laden
#     STRUKTUR: models_all[[modelltyp]][[nischenbreite]]
#
#     Neue Modelle einfach als weiteren Block hinzufügen - der Rest läuft
#     automatisch durch.
# -----------------------------------------------------------------------------

niche_names <- c("narrow", "low_mid", "high_mid", "broad")

# --- GLM ---------------------------------------------------------------------
models_glm <- list(
  narrow   = readRDS("Data/models/knndm/glm_narrow_train.RDS"),
  low_mid  = readRDS("Data/models/knndm/glm_low_mid_train.RDS"),
  high_mid = readRDS("Data/models/knndm/glm_high_mid_train.RDS"),
  broad    = readRDS("Data/models/knndm/glm_broad_train.RDS")
)

# --- Maxent ---------------------------------------------------------------------
models_maxent <- list(
  narrow   = readRDS("Data/models/knndm/maxnet_narrow_train.RDS"),
  low_mid  = readRDS("Data/models/knndm/maxnet_low_mid_train.RDS"),
  high_mid = readRDS("Data/models/knndm/maxnet_high_mid_train.RDS"),
  broad    = readRDS("Data/models/knndm/maxnet_broad_train.RDS")
)

# --- BRT ---------------------------------------------------------------------
models_brt <- list(
  narrow   = readRDS("Data/models/knndm/brt_narrow_train.RDS"),
  low_mid  = readRDS("Data/models/knndm/brt_low_mid_train.RDS"),
  high_mid = readRDS("Data/models/knndm/brt_high_mid_train.RDS"),
  broad    = readRDS("Data/models/knndm/brt_broad_train.RDS")
)

# --- RF ----------------------------------------------------------------------
models_rf <- list(
  narrow   = readRDS("Data/models/knndm/rf_narrow_train.RDS"),
  low_mid  = readRDS("Data/models/knndm/rf_low_mid_train.RDS"),
  high_mid = readRDS("Data/models/knndm/rf_high_mid_train.RDS"),
  broad    = readRDS("Data/models/knndm/rf_broad_train.RDS")
)

# --- Alle Modelltypen zusammenfassen -----------------------------------------
# Hier einfach weitere Modelltypen eintragen (Namen werden später für den Plot
# verwendet)
models_all <- list(
  GLM = models_glm,
  Maxent = models_maxent,
  BRT = models_brt,
  RF  = models_rf
)

# --- Testdaten (einmal laden, für alle Modelle gleich) -----------------------
test_data <- list(
  narrow   = readRDS("Data/species/split_knndm/test_narrow.RDS"),
  low_mid  = readRDS("Data/species/split_knndm/test_low_mid.RDS"),
  high_mid = readRDS("Data/species/split_knndm/test_high_mid.RDS"),
  broad    = readRDS("Data/species/split_knndm/test_broad.RDS")
)

# -----------------------------------------------------------------------------
# 3 - Generische Predict-Funktion
#     Jeder Modelltyp hat u.U. eine andere predict()-Syntax.
#     Hier zentral hinterlegen, damit eval_one_niche() sauber bleibt.
# -----------------------------------------------------------------------------

predict_model <- function(model_obj, newdata, model_type) {
  
  switch(model_type,
         
         # GLM: einzelnes Modell unter $model
         "GLM" = predict(model_obj$model,
                         newdata = newdata,
                         type    = "response"),
         
         # Maxent: einzelnes Modell, aber Daten müssen skaliert werden!
         "Maxent" = {
           # ✅ Skalierungsparameter aus dem Modellobjekt
           scale_center <- model_obj$scale_center
           scale_scale  <- model_obj$scale_scale
           env_vars     <- names(scale_center)
           
           # ✅ Testdaten skalieren
           newdata_scaled <- newdata
           newdata_scaled[, env_vars] <- scale(
             newdata[, env_vars],
             center = scale_center,
             scale  = scale_scale
           )
           
           # ✅ Vorhersage mit skalierten Daten
           as.vector(predict(model_obj$model,
                             newdata = newdata_scaled,
                             type    = "cloglog"))
         },
         
         # BRT (gbm): Ensemble aus $models
         "BRT" = {
           gbm_models <- model_obj$models
           n_trees    <- gbm_models[[1]]$n.trees
           preds <- sapply(gbm_models, function(gbm_fit) {
             gbm::predict.gbm(gbm_fit,
                              newdata = newdata,
                              n.trees = n_trees,
                              type    = "response")
           })
           rowMeans(preds)
         },
         
         # RF (ranger): Ensemble aus $models
         # ✅ Fix: Spalte "1" (Präsenz) verwenden
         "RF" = {
           ranger_models <- model_obj$models
           preds <- sapply(ranger_models, function(rf_fit) {
             pred_list <- predict(rf_fit, data = newdata, type = "response")
             # ✅ Spalte "1" = Präsenz-Wahrscheinlichkeit
             pred_list$predictions[, "1"]
           })
           rowMeans(preds)
         },
         
         # Fallback
         predict(model_obj, newdata = newdata, type = "response")
  )
}

# -----------------------------------------------------------------------------
# 4 - Evaluierungsfunktion (jetzt modelltyp-agnostisch)
# -----------------------------------------------------------------------------

eval_one_niche <- function(model_obj, test_df, env,
                           model_type    = "GLM",
                           thresh.method = "ObsPrev",
                           debug         = FALSE) {
  
  # Nur Testzeilen
  test_df <- test_df[test_df$split == "test", ]
  
  # Umweltwerte extrahieren
  env_vals <- terra::extract(env, test_df[, c("x", "y")])[, -1, drop = FALSE]
  newdat   <- cbind(test_df, env_vals)
  
  # Vorhersage über die generische Funktion
  pred <- tryCatch(
    predict_model(model_obj, newdat, model_type),
    error = function(e) {
      warning(sprintf("[%s] predict() fehlgeschlagen: %s", model_type, e$message))
      rep(NA_real_, nrow(newdat))
    }
  )
  
  # ✅ Sicherstellen, dass pred gleich lang wie newdat ist
  # Falls predict() intern NA-Zeilen entfernt (z.B. maxnet), auffüllen
  if (length(pred) != nrow(newdat)) {
    warning(sprintf(
      "[%s] pred (%d) kürzer als newdat (%d) – fülle fehlende Zeilen mit NA auf.",
      model_type, length(pred), nrow(newdat)
    ))
    # Finde NA-Zeilen in newdat (Umweltvariablen)
    env_cols <- names(env_vals)
    na_rows  <- apply(newdat[, env_cols, drop = FALSE], 1, anyNA)
    pred_full <- rep(NA_real_, nrow(newdat))
    pred_full[!na_rows] <- pred
    pred <- pred_full
  }
  
  # NA-Zeilen entfernen
  ok <- !is.na(pred)
  if (any(!ok)) {
    warning(sprintf("[%s] %d von %d Testpunkten ausgeschlossen (NA).",
                    model_type, sum(!ok), length(ok)))
  }
  
  # Absicherung: mindestens beide Klassen vorhanden
  if (length(unique(newdat$presence[ok])) < 2) {
    warning(sprintf("[%s] Nur eine Klasse im Testdatensatz - Evaluation übersprungen.",
                    model_type))
    return(data.frame(AUC = NA, TSS = NA, thresh = NA,
                      true_prevalence = NA, n_test = sum(ok)))
  }
  
  # Reine Vektoren sicherstellen
  observation <- as.vector(newdat$presence[ok])
  predictions <- as.vector(pred[ok])
  
  # Debug
  if (debug) {
    cat("nrow(newdat):", nrow(newdat), "\n")
    cat("length(pred):", length(pred), "\n")
    cat("sum(ok):", sum(ok), "\n")
    cat("length(observation):", length(observation), "\n")
    cat("length(predictions):", length(predictions), "\n")
  }
  
  eval_res <- evalSDM(observation = observation,
                      predictions = predictions,
                      thresh.method = thresh.method)
  
  eval_res$true_prevalence <- mean(newdat$presence[ok])
  eval_res$n_test          <- sum(ok)
  
  eval_res
}

# -----------------------------------------------------------------------------
# 5 - Bootstrap-KI (unverändert, aber jetzt mit model_type-Argument)
# -----------------------------------------------------------------------------

bootstrap_tss <- function(obs, pred, thresh.method = "ObsPrev", n_boot = 1000) {
  n        <- length(obs)
  tss_boot <- numeric(n_boot)
  
  for (i in seq_len(n_boot)) {
    idx   <- sample(seq_len(n), size = n, replace = TRUE)
    obs_b <- obs[idx]
    pred_b <- pred[idx]
    if (length(unique(obs_b)) < 2) next
    res <- tryCatch(
      evalSDM(observation = obs_b, predictions = pred_b,
              thresh.method = thresh.method),
      error = function(e) NULL
    )
    tss_boot[i] <- if (!is.null(res)) res$TSS else NA
  }
  quantile(tss_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
}


# -----------------------------------------------------------------------------
# 6 - Hauptschleife: alle Modelltypen × alle Nischenbreiten
# -----------------------------------------------------------------------------

all_results <- list()
all_ci      <- list()

for (mtype in names(models_all)) {
  
  cat("\n==============================\n")
  cat("Evaluiere Modelltyp:", mtype, "\n")
  cat("==============================\n")
  
  local({
    mtype_local  <- mtype
    models_mtype <- models_all[[mtype_local]]
    
    res_list <- Map(
      function(mod, tdat) eval_one_niche(mod, tdat, env,
                                         model_type    = mtype_local,
                                         thresh.method = "ObsPrev",
                                         debug         = FALSE),
      models_mtype, test_data
    )
    
    res_df <- do.call(rbind, res_list)
    res_df$niche_breadth <- niche_names
    res_df$model_type    <- mtype_local
    rownames(res_df)     <- NULL
    all_results[[mtype_local]] <<- res_df
    
    ci_list <- lapply(niche_names, function(nm) {
      test_df  <- test_data[[nm]][test_data[[nm]]$split == "test", ]
      env_vals <- terra::extract(env, test_df[, c("x", "y")])[, -1, drop = FALSE]
      newdat   <- cbind(test_df, env_vals)
      pred     <- tryCatch(
        predict_model(models_mtype[[nm]], newdat, mtype_local),
        error = function(e) rep(NA_real_, nrow(newdat))
      )
      # ✅ Auch hier: pred auffüllen falls nötig
      if (length(pred) != nrow(newdat)) {
        env_cols  <- names(env_vals)
        na_rows   <- apply(newdat[, env_cols, drop = FALSE], 1, anyNA)
        pred_full <- rep(NA_real_, nrow(newdat))
        pred_full[!na_rows] <- pred
        pred <- pred_full
      }
      ok <- !is.na(pred)
      bootstrap_tss(as.vector(newdat$presence[ok]),
                    as.vector(pred[ok]))
    })
    names(ci_list)        <- niche_names
    all_ci[[mtype_local]] <<- ci_list
  })
}

# -----------------------------------------------------------------------------
# 7 - Ergebnisse zusammenführen
# -----------------------------------------------------------------------------

results_combined <- do.call(rbind, all_results)
print(results_combined)

# CI in ein übersichtliches Data Frame umwandeln
ci_df <- do.call(rbind, lapply(names(all_ci), function(mtype) {
  do.call(rbind, lapply(names(all_ci[[mtype]]), function(nm) {
    ci <- all_ci[[mtype]][[nm]]
    data.frame(
      model_type   = mtype,
      niche_breadth = nm,
      TSS_lo       = ci["2.5%"],
      TSS_med      = ci["50%"],
      TSS_hi       = ci["97.5%"],
      row.names    = NULL
    )
  }))
}))

# Für den Plot: results_combined mit CI zusammenführen
plot_df <- results_combined %>%
  left_join(ci_df, by = c("model_type", "niche_breadth")) %>%
  mutate(
    niche_breadth = factor(niche_breadth, levels = niche_names),
    model_type    = factor(model_type,    levels = names(models_all))
  )

# -----------------------------------------------------------------------------
# 8 - Speichern
# -----------------------------------------------------------------------------

write.csv(results_combined, "results/all_models_auc_tss.csv", row.names = FALSE)
write.csv(ci_df,            "results/all_models_tss_ci.csv",  row.names = FALSE)

# -----------------------------------------------------------------------------
# 9 - Vergleichsplots
# -----------------------------------------------------------------------------




# Farbpalette (automatisch skaliert mit Anzahl Modelltypen)
n_models <- length(unique(plot_df$model_type))
farben   <- RColorBrewer::brewer.pal(max(3, n_models), "Set1")[seq_len(n_models)]

# --- 9a: TSS-Vergleich mit Bootstrap-KI ------------------------------------
p_tss <- ggplot(plot_df,
                aes(x     = niche_breadth,
                    y     = TSS,
                    color = model_type,
                    group = model_type)) +
  
  # Konfidenzband
  geom_ribbon(aes(ymin  = TSS_lo,
                  ymax  = TSS_hi,
                  fill  = model_type),
              alpha        = 0.15,
              color        = NA) +
  
  # Linie + Punkte
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  
  # Referenzlinie: TSS = 0 (kein Skill)
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             color      = "grey50") +
  
  scale_color_manual(values = farben) +
  scale_fill_manual(values  = farben) +
  
  labs(
    title    = "TSS-Vergleich: Modelltypen × Nischenbreite",
    subtitle = "Punkte = TSS; Bänder = 95%-Bootstrap-KI",
    x        = "Nischenbreite",
    y        = "TSS",
    color    = "Modelltyp",
    fill     = "Modelltyp"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

print(p_tss)

# --- 9b: AUC-Vergleich (ohne KI, da AUC-Bootstrap hier nicht berechnet) ----
p_auc <- ggplot(plot_df,
                aes(x     = niche_breadth,
                    y     = AUC,
                    color = model_type,
                    group = model_type)) +
  
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  
  # Referenzlinie: AUC = 0.5 (Zufallsmodell)
  geom_hline(yintercept = 0.5,
             linetype   = "dashed",
             color      = "grey50") +
  
  scale_color_manual(values = farben) +
  
  labs(
    title = "AUC-Vergleich: Modelltypen × Nischenbreite",
    x     = "Nischenbreite",
    y     = "AUC",
    color = "Modelltyp"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

print(p_auc)

# --- 9c: TSS + AUC kombiniert (Facets) --------------------------------------
plot_long <- plot_df %>%
  select(model_type, niche_breadth, TSS, AUC) %>%
  pivot_longer(cols      = c(TSS, AUC),
               names_to  = "metric",
               values_to = "value")

p_combined <- ggplot(plot_long,
                     aes(x     = niche_breadth,
                         y     = value,
                         color = model_type,
                         group = model_type)) +
  
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  
  facet_wrap(~ metric,
             scales = "free_y",
             nrow   = 1) +
  
  scale_color_manual(values = farben) +
  
  labs(
    title = "Modellvergleich: TSS und AUC",
    x     = "Nischenbreite",
    y     = "Wert",
    color = "Modelltyp"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

print(p_combined)

# --- Plots speichern --------------------------------------------------------
ggsave("results/plot_tss_vergleich.png",  p_tss,      width = 8, height = 5, dpi = 300)
ggsave("results/plot_auc_vergleich.png",  p_auc,      width = 8, height = 5, dpi = 300)
ggsave("results/plot_combined.png",       p_combined, width = 12, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 10 - Presence-Übersicht (einmalig, da Testdaten für alle Modelle gleich)
# -----------------------------------------------------------------------------

presence_overview <- data.frame(
  niche_breadth    = niche_names,
  n_presence_test  = sapply(test_data, function(d)
    sum(d$presence[d$split == "test"] == 1)),
  n_absence_test   = sapply(test_data, function(d)
    sum(d$presence[d$split == "test"] == 0))
)
presence_overview$prevalence_test <-
  presence_overview$n_presence_test /
  (presence_overview$n_presence_test + presence_overview$n_absence_test)

print(presence_overview)


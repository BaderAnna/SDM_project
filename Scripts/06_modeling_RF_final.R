# =============================================================================
# RANDOM FOREST MODELLIERUNG
# WorldClim Variablen nach VIF-Reduktion
# Spatial CV: CAST::knndm, 1 Durchlauf
# 10 Runs pro Art (unterschiedliche pseudo-absences), dann mitteln
# Für narrow, low_mid, high_mid, broad
# =============================================================================

library(terra)
library(sf)
library(ranger)   # schneller als randomForest
library(CAST)
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

# =============================================================================
# 1 - Raster + Variablen laden
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars     <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - Hilfsfunktionen
# =============================================================================

# --- Daten vorbereiten für einen RF-Run ---
# Koordinaten bereits EPSG:3035 -> keine Transformation!
prepare_rf_data <- function(run_data, env_raster) {
  
  coords <- as.matrix(run_data[, c("x", "y")])
  
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  dat <- cbind(run_data, env_vals)
  dat <- na.omit(dat)
  
  # presence als Faktor für Klassifikation
  dat$presence <- as.factor(dat$presence)
  
  return(dat)
}

# --- knndm Fold-Erstellung ---
run_knndm_folds <- function(dat_sf_geom, predictor_raster,
                            k = 5, seed = 42, n_predpoints = 2000) {
  set.seed(seed)
  
  pred_pts <- terra::spatSample(
    predictor_raster, size = n_predpoints,
    method = "random", as.points = TRUE, na.rm = TRUE
  ) |> sf::st_as_sf()
  
  pred_pts_geom <- sf::st_as_sf(
    data.frame(geometry = sf::st_geometry(pred_pts)),
    crs = sf::st_crs(pred_pts)
  )
  
  target_crs    <- sf::st_crs(predictor_raster)
  dat_sf_geom   <- sf::st_transform(dat_sf_geom,   target_crs)
  pred_pts_geom <- sf::st_transform(pred_pts_geom, target_crs)
  
  CAST::knndm(
    tpoints    = dat_sf_geom,
    predpoints = pred_pts_geom,
    k          = k
  )
}

# =============================================================================
# 3 - RF mit knndm Spatial-CV
#     10 Runs mit unterschiedlichen pseudo-absences -> mitteln
# =============================================================================

run_rf_knndm <- function(sampling, env_raster,
                         n_runs = 10, k = 5,
                         n_trees = 500) {
  
  env_cols   <- names(env_raster)
  rf_formula <- as.formula(paste("presence ~",
                                 paste(env_cols, collapse = " + ")))
  
  auc_scores <- numeric(n_runs)
  tss_scores <- numeric(n_runs)
  pred_list  <- list()
  model_list <- list()
  
  for (run in seq_len(n_runs)) {
    
    message(sprintf("  -> Run %d/%d", run, n_runs))
    
    # --- Daten für diesen Run vorbereiten ---
    run_data <- prepare_rf_data(sampling$runs[[run]], env_raster)
    
    if (nrow(run_data) == 0) {
      message("  -> Run ", run, ": dat leer -> überspringe")
      next
    }
    
    message(sprintf("     %d Zeilen | %d Präsenzen | %d Absences",
                    nrow(run_data),
                    sum(run_data$presence == 1),
                    sum(run_data$presence == 0)))
    
    # sf-Objekt für knndm - NUR Geometrie
    dat_sf      <- sf::st_as_sf(run_data, coords = c("x", "y"), crs = 3035)
    dat_sf_geom <- sf::st_as_sf(
      data.frame(geometry = sf::st_geometry(dat_sf)),
      crs = sf::st_crs(dat_sf)
    )
    
    # --- knndm Folds ---
    knn <- tryCatch(
      run_knndm_folds(dat_sf_geom, env_raster, k = k, seed = run),
      error = function(e) {
        message("  knndm Fehler in Run ", run, ": ", e$message)
        NULL
      }
    )
    
    if (!is.null(knn)) {
      
      fold_ids    <- knn$clusters
      fold_preds  <- rep(NA_real_, nrow(run_data))
      fold_truths <- rep(NA_real_, nrow(run_data))
      
      for (fold in seq_len(k)) {
        
        test_idx  <- which(fold_ids == fold)
        train_idx <- which(fold_ids != fold)
        if (length(test_idx) == 0 || length(train_idx) == 0) next
        
        train <- run_data[train_idx, ]
        test  <- run_data[test_idx,  ]
        if (length(unique(train$presence)) < 2) next
        
        # RF fitten mit ranger
        rf_model <- tryCatch(
          ranger::ranger(
            rf_formula,
            data        = train,
            num.trees   = n_trees,
            probability = TRUE,   # für Wahrscheinlichkeiten
            seed        = run * fold
          ),
          error = function(e) {
            message("  RF Fehler in Fold ", fold, ": ", e$message)
            NULL
          }
        )
        if (is.null(rf_model)) next
        
        # Vorhersage: Wahrscheinlichkeit für Klasse "1"
        fold_preds[test_idx] <- predict(
          rf_model, data = test
        )$predictions[, "1"]
        
        fold_truths[test_idx] <- as.numeric(
          as.character(test$presence)
        )
      }
      
      # Metriken berechnen
      valid        <- !is.na(fold_preds) & !is.na(fold_truths)
      preds_valid  <- fold_preds[valid]
      truths_valid <- fold_truths[valid]
      
      if (length(preds_valid) >= 10 && length(unique(truths_valid)) == 2) {
        
        auc <- as.numeric(
          pROC::auc(pROC::roc(truths_valid, preds_valid, quiet = TRUE))
        )
        
        pred_bin <- ifelse(preds_valid >= 0.5, 1, 0)
        tp   <- sum(pred_bin == 1 & truths_valid == 1)
        fp   <- sum(pred_bin == 1 & truths_valid == 0)
        tn   <- sum(pred_bin == 0 & truths_valid == 0)
        fn   <- sum(pred_bin == 0 & truths_valid == 1)
        sens <- tp / (tp + fn)
        spec <- tn / (tn + fp)
        tss  <- sens + spec - 1
        
        auc_scores[run] <- auc
        tss_scores[run] <- tss
        
        message(sprintf("     AUC = %.3f | TSS = %.3f", auc, tss))
      }
    }
    
    # --- Finales RF-Modell für diesen Run auf ALLEN Daten ---
    rf_final <- tryCatch(
      ranger::ranger(
        rf_formula,
        data        = run_data,
        num.trees   = n_trees,
        probability = TRUE,
        seed        = run
      ),
      error = function(e) {
        message("  RF Final-Fehler in Run ", run, ": ", e$message)
        NULL
      }
    )
    if (is.null(rf_final)) next
    
    # Vorhersage auf Raster
    pred_run <- terra::predict(
      env_raster,
      rf_final,
      fun = function(model, data) {
        predict(model, data = data)$predictions[, "1"]
      },
      na.rm = TRUE
    )
    
    model_list[[run]] <- rf_final
    pred_list[[run]]  <- pred_run
  }
  
  # --- Über alle Runs mitteln ---
  message("  -> Mittele über ", length(pred_list), " Runs...")
  final_pred <- terra::app(terra::rast(pred_list), mean)
  
  # Nur gültige Scores
  valid_auc <- auc_scores[auc_scores != 0]
  valid_tss <- tss_scores[tss_scores != 0]
  
  message(sprintf("  -> AUC mean = %.3f (sd = %.3f) | TSS mean = %.3f (sd = %.3f)",
                  mean(valid_auc, na.rm = TRUE),
                  sd(valid_auc,   na.rm = TRUE),
                  mean(valid_tss, na.rm = TRUE),
                  sd(valid_tss,   na.rm = TRUE)))
  
  list(
    models     = model_list,
    prediction = final_pred,
    auc_mean   = mean(valid_auc, na.rm = TRUE),
    auc_sd     = sd(valid_auc,   na.rm = TRUE),
    tss_mean   = mean(valid_tss, na.rm = TRUE),
    tss_sd     = sd(valid_tss,   na.rm = TRUE),
    auc_all    = auc_scores,
    tss_all    = tss_scores
  )
}

# =============================================================================
# 4 - Sampling-Daten laden + RF ausführen
# =============================================================================

sampling_narrow   <- readRDS("Data/species/sampling_narrow.RDS")
sampling_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
sampling_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
sampling_broad    <- readRDS("Data/species/sampling_broad.RDS")

message("===== Narrow =====")
rf_narrow   <- run_rf_knndm(sampling_narrow,   env_raster_masked, n_runs = 10)

message("===== Low-Mid =====")
rf_low_mid  <- run_rf_knndm(sampling_low_mid,  env_raster_masked, n_runs = 10)

message("===== High-Mid =====")
rf_high_mid <- run_rf_knndm(sampling_high_mid, env_raster_masked, n_runs = 10)

message("===== Broad =====")
rf_broad    <- run_rf_knndm(sampling_broad,    env_raster_masked, n_runs = 10)

# =============================================================================
# 5 - Ergebnisse zusammenfassen
# =============================================================================

results_rf <- data.frame(
  species  = c("narrow", "low_mid", "high_mid", "broad"),
  sigma    = c(0.2, 0.4, 0.6, 0.8),
  auc_mean = c(rf_narrow$auc_mean,
               rf_low_mid$auc_mean,
               rf_high_mid$auc_mean,
               rf_broad$auc_mean),
  auc_sd   = c(rf_narrow$auc_sd,
               rf_low_mid$auc_sd,
               rf_high_mid$auc_sd,
               rf_broad$auc_sd),
  tss_mean = c(rf_narrow$tss_mean,
               rf_low_mid$tss_mean,
               rf_high_mid$tss_mean,
               rf_broad$tss_mean),
  tss_sd   = c(rf_narrow$tss_sd,
               rf_low_mid$tss_sd,
               rf_high_mid$tss_sd,
               rf_broad$tss_sd)
)

print(results_rf)

# =============================================================================
# 6 - Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
plot(rf_narrow$prediction,
     main = sprintf("RF Narrow\nAUC=%.3f | TSS=%.3f",
                    rf_narrow$auc_mean, rf_narrow$tss_mean))
plot(rf_low_mid$prediction,
     main = sprintf("RF Low-Mid\nAUC=%.3f | TSS=%.3f",
                    rf_low_mid$auc_mean, rf_low_mid$tss_mean))
plot(rf_high_mid$prediction,
     main = sprintf("RF High-Mid\nAUC=%.3f | TSS=%.3f",
                    rf_high_mid$auc_mean, rf_high_mid$tss_mean))
plot(rf_broad$prediction,
     main = sprintf("RF Broad\nAUC=%.3f | TSS=%.3f",
                    rf_broad$auc_mean, rf_broad$tss_mean))
par(mfrow = c(1, 1))

# =============================================================================
# 7 - Speichern
# =============================================================================

dir.create("Data/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(rf_narrow,   "Data/models/rf_narrow.RDS")
saveRDS(rf_low_mid,  "Data/models/rf_low_mid.RDS")
saveRDS(rf_high_mid, "Data/models/rf_high_mid.RDS")
saveRDS(rf_broad,    "Data/models/rf_broad.RDS")
saveRDS(results_rf,  "Data/models/results_rf.RDS")

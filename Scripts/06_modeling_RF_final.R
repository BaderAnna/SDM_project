# =============================================================================
# RANDOM FOREST MODELLIERUNG
# WorldClim Variablen nach VIF-Reduktion
# 10 Runs pro Art (unterschiedliche pseudo-absences), dann mitteln
# Reines Modelltraining (Evaluierung erfolgt in separatem Skript)
# =============================================================================

library(terra)
library(sf)
library(ranger)   # schneller als randomForest

# =============================================================================
# 1 - Raster + Variablen laden
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars     <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - Daten vorbereiten für einen RF-Run
#     Koordinaten bereits EPSG:3035 -> keine Transformation!
# =============================================================================

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

# =============================================================================
# 3 - RF trainieren
#     10 Runs mit unterschiedlichen pseudo-absences -> Vorhersagen mitteln
#     (kein CV/Evaluierung - nur zur Absicherung gegen Pseudo-Absence-Zufall)
# =============================================================================

run_rf <- function(sampling, env_raster, n_runs = 10, n_trees = 500) {
  
  env_cols   <- names(env_raster)
  rf_formula <- as.formula(paste("presence ~",
                                 paste(env_cols, collapse = " + ")))
  
  pred_list  <- list()
  model_list <- list()
  
  for (run in seq_len(n_runs)) {
    
    message(sprintf("  -> Run %d/%d", run, n_runs))
    
    run_data <- prepare_rf_data(sampling$runs[[run]], env_raster)
    
    if (nrow(run_data) == 0) {
      message("  -> Run ", run, ": dat leer -> überspringe")
      next
    }
    
    message(sprintf("     %d Zeilen | %d Präsenzen | %d Absences",
                    nrow(run_data),
                    sum(run_data$presence == 1),
                    sum(run_data$presence == 0)))
    
    # Finales RF-Modell für diesen Run auf ALLEN Daten
    rf_final <- tryCatch(
      ranger::ranger(
        rf_formula,
        data        = run_data,
        num.trees   = n_trees,
        probability = TRUE,
        seed        = run
      ),
      error = function(e) {
        message("  RF Fehler in Run ", run, ": ", e$message)
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
  
  list(
    models     = model_list,
    prediction = final_pred
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
rf_narrow   <- run_rf(sampling_narrow,   env_raster_masked, n_runs = 10)

message("===== Low-Mid =====")
rf_low_mid  <- run_rf(sampling_low_mid,  env_raster_masked, n_runs = 10)

message("===== High-Mid =====")
rf_high_mid <- run_rf(sampling_high_mid, env_raster_masked, n_runs = 10)

message("===== Broad =====")
rf_broad    <- run_rf(sampling_broad,    env_raster_masked, n_runs = 10)

# =============================================================================
# 5 - Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
plot(rf_narrow$prediction,   main = "RF Narrow")
plot(rf_low_mid$prediction,  main = "RF Low-Mid")
plot(rf_high_mid$prediction, main = "RF High-Mid")
plot(rf_broad$prediction,    main = "RF Broad")
par(mfrow = c(1, 1))

# =============================================================================
# 6 - Speichern
# =============================================================================

dir.create("Data/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(rf_narrow,   "Data/models/rf_narrow.RDS")
saveRDS(rf_low_mid,  "Data/models/rf_low_mid.RDS")
saveRDS(rf_high_mid, "Data/models/rf_high_mid.RDS")
saveRDS(rf_broad,    "Data/models/rf_broad.RDS")

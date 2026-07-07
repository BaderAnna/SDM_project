# =============================================================================
# BRT MODELLIERUNG
# WorldClim Variablen nach VIF-Reduktion
# 10 Runs pro Art (unterschiedliche pseudo-absences), dann mitteln
# Reines Modelltraining (Evaluierung erfolgt in separatem Skript)
# =============================================================================

library(terra)
library(sf)
library(gbm)

# =============================================================================
# 1 - Raster + Variablen laden
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars     <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - Formel für BRT bauen
# =============================================================================

build_brt_formula <- function(vars) {
  as.formula(paste("presence ~", paste(vars, collapse = " + ")))
}

brt_formula <- build_brt_formula(selected_vars)
print(brt_formula)

# =============================================================================
# 3 - Daten vorbereiten für einen BRT-Run
#     Koordinaten bereits EPSG:3035 -> keine Transformation!
# =============================================================================

prepare_brt_data <- function(run_data, env_raster) {
  
  coords <- as.matrix(run_data[, c("x", "y")])
  
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  dat <- cbind(run_data, env_vals)
  dat <- na.omit(dat)
  
  return(dat)
}

# =============================================================================
# 4 - BRT trainieren
#     10 Runs mit unterschiedlichen pseudo-absences -> Vorhersagen mitteln
#     (kein CV/Evaluierung - nur zur Absicherung gegen Pseudo-Absence-Zufall)
# =============================================================================

run_brt <- function(sampling, env_raster, brt_formula,
                    n_runs = 10,
                    n_trees = 2000, interaction_depth = 3,
                    shrinkage = 0.01, bag_fraction = 0.75) {
  
  pred_list  <- list()
  model_list <- list()
  
  for (run in seq_len(n_runs)) {
    
    message(sprintf("  -> Run %d/%d", run, n_runs))
    
    # runs enthält 10 verschiedene pseudo-absence Datensätze
    run_data <- prepare_brt_data(sampling$runs[[run]], env_raster)
    
    if (nrow(run_data) == 0) {
      message("  -> Run ", run, ": dat leer -> überspringe")
      next
    }
    
    message(sprintf("     %d Zeilen | %d Präsenzen | %d Absences",
                    nrow(run_data),
                    sum(run_data$presence == 1),
                    sum(run_data$presence == 0)))
    
    # BRT fitten (cv.folds hier nur zur internen Bestimmung der optimalen
    # Baumanzahl via gbm.perf, NICHT zur externen Modellevaluierung)
    brt_final <- tryCatch(
      gbm::gbm(
        brt_formula,
        data              = run_data,
        distribution      = "bernoulli",
        n.trees           = n_trees,
        interaction.depth = interaction_depth,
        shrinkage         = shrinkage,
        bag.fraction      = bag_fraction,
        cv.folds          = 5,
        verbose           = FALSE
      ),
      error = function(e) {
        message("  BRT Fehler in Run ", run, ": ", e$message)
        NULL
      }
    )
    if (is.null(brt_final)) next
    
    best_trees <- gbm::gbm.perf(brt_final, method = "cv", plot.it = FALSE)
    
    # Vorhersage auf Raster
    pred_run <- terra::predict(
      env_raster,
      brt_final,
      n.trees = best_trees,
      type    = "response",
      na.rm   = TRUE
    )
    
    model_list[[run]] <- brt_final
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
# 5 - Sampling-Daten laden + BRT ausführen
# =============================================================================

sampling_narrow   <- readRDS("Data/species/sampling_narrow.RDS")
sampling_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
sampling_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
sampling_broad    <- readRDS("Data/species/sampling_broad.RDS")

message("===== Narrow =====")
brt_narrow   <- run_brt(sampling_narrow,   env_raster_masked, brt_formula, n_runs = 10)

message("===== Low-Mid =====")
brt_low_mid  <- run_brt(sampling_low_mid,  env_raster_masked, brt_formula, n_runs = 10)

message("===== High-Mid =====")
brt_high_mid <- run_brt(sampling_high_mid, env_raster_masked, brt_formula, n_runs = 10)

message("===== Broad =====")
brt_broad    <- run_brt(sampling_broad,    env_raster_masked, brt_formula, n_runs = 10)

# =============================================================================
# 6 - Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
plot(brt_narrow$prediction,   main = "BRT Narrow")
plot(brt_low_mid$prediction,  main = "BRT Low-Mid")
plot(brt_high_mid$prediction, main = "BRT High-Mid")
plot(brt_broad$prediction,    main = "BRT Broad")
par(mfrow = c(1, 1))

# =============================================================================
# 7 - Speichern
# =============================================================================

dir.create("Data/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(brt_narrow,   "Data/models/brt_narrow.RDS")
saveRDS(brt_low_mid,  "Data/models/brt_low_mid.RDS")
saveRDS(brt_high_mid, "Data/models/brt_high_mid.RDS")
saveRDS(brt_broad,    "Data/models/brt_broad.RDS")
# =============================================================================
# BRT MODELLIERUNG (TRAINING)
# WorldClim Variablen nach VIF-Reduktion
# Nutzt die im Skript "02_VirtualSpecies_presence_CV.R" (split_data_knndm())
# bereits erzeugten Trainingsdaten (train_brt_<sp>.RDS = Liste von 10
# data.frames mit unterschiedlichen Pseudo-Absence-Sets pro Art, jeweils
# mit Spalte "split" für train/test) -> 10 Runs, dann mitteln.
# Evaluierung erfolgt in separatem Skript auf Basis von test_<sp>.RDS.
# =============================================================================

library(terra)
library(sf)
library(gbm)

# =============================================================================
# 1 - Raster + Variablen laden
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars      <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - BRT-Formel dynamisch bauen
# =============================================================================

build_brt_formula <- function(vars) {
  as.formula(paste("presence ~", paste(vars, collapse = " + ")))
}

brt_formula <- build_brt_formula(selected_vars)
print(brt_formula)

# =============================================================================
# 3 - Trainingsdaten vorbereiten (pro Run)
#     (Fold-Split + Pseudo-Absence-Sets kommen bereits aus train_brt_<sp>.RDS
#      -> hier nur noch Umweltwerte an den Koordinaten extrahieren)
# =============================================================================

prepare_brt_train_data <- function(run_df, env_raster) {
  
  # run_df enthält train + test (Spalte "split") -> nur Trainingsdaten nutzen
  if ("split" %in% names(run_df)) {
    run_df <- run_df[run_df$split == "train", , drop = FALSE]
  }
  
  if (nrow(run_df) == 0) stop("run_df ist leer!")
  if (any(is.na(run_df$x)) || any(is.na(run_df$y))) {
    stop("x oder y enthalten NA-Werte!")
  }
  
  coords_3035  <- as.matrix(run_df[, c("x", "y")])
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords_3035))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  dat <- cbind(run_df, env_vals)
  dat <- na.omit(dat)
  
  message("prepare_brt_train_data: ", nrow(dat), " Zeilen | ",
          sum(dat$presence == 1), " Presence | ",
          sum(dat$presence == 0), " Pseudo-Absence")
  
  return(dat)
}

# =============================================================================
# 4 - BRT auf Trainingsdaten fitten
#     10 Runs mit unterschiedlichen pseudo-absences (train_df ist die Liste
#     der 10 Runs) -> Vorhersagen mitteln
#     (kein CV/Evaluierung - nur zur Absicherung gegen Pseudo-Absence-Zufall)
# =============================================================================

run_brt_train <- function(train_df, env_raster, brt_formula, seed = 42,
                          n_runs = 10, n_trees = 2000, interaction_depth = 3,
                          shrinkage = 0.01, bag_fraction = 0.75) {
  
  n_runs <- min(n_runs, length(train_df))
  
  set.seed(seed)
  
  pred_list  <- list()
  model_list <- list()
  
  for (run in seq_len(n_runs)) {
    
    message(sprintf("  -> Run %d/%d", run, n_runs))
    
    dat <- prepare_brt_train_data(train_df[[run]], env_raster)
    
    if (nrow(dat) == 0) {
      message("  -> Run ", run, ": dat leer -> überspringe")
      next
    }
    
    # Prüfe, ob alle Variablen der Formel in dat existieren
    missing_vars <- setdiff(all.vars(brt_formula), names(dat))
    if (length(missing_vars) > 0) {
      stop("Fehlende Variablen in dat: ", paste(missing_vars, collapse = ", "))
    }
    
    # BRT fitten (cv.folds hier nur zur internen Bestimmung der optimalen
    # Baumanzahl via gbm.perf, NICHT zur externen Modellevaluierung)
    model <- withCallingHandlers(
      tryCatch(
        gbm::gbm(
          brt_formula,
          data              = dat,
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
      ),
      warning = function(w) {
        message("Warnung beim Modelltraining (Run ", run, "): ", w$message)
      }
    )
    if (is.null(model)) next
    
    best_trees <- gbm::gbm.perf(model, method = "cv", plot.it = FALSE)
    
    prediction <- terra::predict(
      env_raster,
      model,
      n.trees = best_trees,
      type    = "response",
      na.rm   = TRUE
    )
    
    model_list[[run]] <- model
    pred_list[[run]]  <- prediction
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
# 5 - Für alle Arten ausführen
# =============================================================================

species_list     <- c("narrow", "low_mid", "high_mid", "broad")
brt_models_train  <- list()

for (sp in species_list) {
  
  message("===== ", sp, " (Training) =====")
  
  train_brt <- readRDS(paste0("Data/species/split_knndm/train_brt_", sp, ".RDS"))
  
  brt_models_train[[sp]] <- run_brt_train(
    train_df    = train_brt,
    env_raster  = env_raster_masked,
    brt_formula = brt_formula,
    n_runs      = 10
  )
}

# =============================================================================
# 6 - Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
for (sp in species_list) {
  plot(brt_models_train[[sp]]$prediction, main = paste("BRT", sp))
}
par(mfrow = c(1, 1))

# =============================================================================
# 7 - Speichern der Trainingsmodelle
# =============================================================================

dir.create("Data/models/knndm", recursive = TRUE, showWarnings = FALSE)

for (sp in species_list) {
  saveRDS(
    brt_models_train[[sp]],
    paste0("Data/models/knndm/brt_", sp, "_train.RDS")
  )
  message("✓ Trainingsmodell gespeichert: brt_", sp, "_train.RDS")
}





'# =============================================================================
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
saveRDS(brt_broad,    "Data/models/brt_broad.RDS")'

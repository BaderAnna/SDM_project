# =============================================================================
# DATEN-SPLITTING MIT kNNDM (80/20)
# kNNDM nur auf Presence Points
# =============================================================================

library(terra)
library(sf)
library(CAST)

# =============================================================================
# 1. Funktion: kNNDM-basiertes Split
# =============================================================================

split_data_knndm <- function(sampling, species_PA,
                             train_ratio = 0.8, k = 5, seed = 42,
                             n_true_abs  = 1000) {
  
  # -------------------------------------------------------------------
  # 1. Presence Points
  # -------------------------------------------------------------------
  pres <- sampling$po$sample.points[
    sampling$po$sample.points$Observed == TRUE, ]
  pres <- data.frame(x = pres$x, y = pres$y, presence = 1)
  
  message("Presence Points: ", nrow(pres))
  
  # -------------------------------------------------------------------
  # 2. Räumliche Geometrie (EPSG:3035)
  # -------------------------------------------------------------------
  pres_sf <- sf::st_as_sf(pres, coords = c("x", "y"), crs = 3035)
  
  # -------------------------------------------------------------------
  # 3. Vorhersagepunkte für kNNDM (aus Raster)
  # -------------------------------------------------------------------
  predictor_raster <- terra::rast("Data/raster/pc_raster_masked.tif")
  predictor_raster <- terra::project(predictor_raster, "EPSG:3035")
  
  set.seed(seed)
  pred_pts <- terra::spatSample(
    predictor_raster,
    size      = 10000,
    method    = "random",
    as.points = TRUE,
    na.rm     = TRUE
  ) |> sf::st_as_sf()
  
  pred_pts <- sf::st_transform(pred_pts, crs = 3035)
  
  # -------------------------------------------------------------------
  # 4. kNNDM NUR auf Presence Points
  # -------------------------------------------------------------------
  knn <- tryCatch(
    CAST::knndm(
      tpoints    = pres_sf,
      predpoints = pred_pts,
      k          = k
    ),
    error = function(e) {
      message("kNNDM fehlgeschlagen: ", e$message)
      stop("kNNDM konnte nicht durchgeführt werden.")
    }
  )
  
  # Fold-Nummern
  fold_ids    <- knn$clusters
  pres$fold   <- fold_ids
  n_folds     <- max(fold_ids)
  
  message("Fold-Verteilung: ", paste(table(fold_ids), collapse = ", "))
  
  # -------------------------------------------------------------------
  # 5. Presence Points splitten (80/20)
  # -------------------------------------------------------------------
  set.seed(seed)
  train_folds <- sample(1:n_folds, size = floor(train_ratio * n_folds))
  test_folds  <- setdiff(1:n_folds, train_folds)
  
  train_pres  <- pres[pres$fold %in% train_folds, c("x", "y", "presence")]
  test_pres   <- pres[pres$fold %in% test_folds,  c("x", "y", "presence")]
  
  message("Train Presence: ", nrow(train_pres))
  message("Test  Presence: ", nrow(test_pres))
  
  # Hilfsstring zum Koordinaten-Matching (Train-Presence identifizieren)
  train_pres_key <- paste(train_pres$x, train_pres$y)
  
  # -------------------------------------------------------------------
  # 6a. GLM: Train-Presence + bg_glm
  # -------------------------------------------------------------------
  bg_glm <- as.data.frame(sampling$bg_glm)
  bg_glm <- data.frame(x = bg_glm$x, y = bg_glm$y, presence = 0)
  
  train_glm       <- rbind(train_pres, bg_glm)
  train_glm$split <- "train"
  
  message("GLM Train: ", nrow(train_glm), " (",
          sum(train_glm$presence == 1), " Presence / ",
          sum(train_glm$presence == 0), " Pseudo-Absence)")
  
  # -------------------------------------------------------------------
  # 6b. MaxEnt: Train-Presence + bg_maxent
  # -------------------------------------------------------------------
  bg_maxent <- as.data.frame(sampling$bg_maxent)
  bg_maxent <- data.frame(x = bg_maxent$x, y = bg_maxent$y, presence = 0)
  
  train_maxent       <- rbind(train_pres, bg_maxent)
  train_maxent$split <- "train"
  
  message("MaxEnt Train: ", nrow(train_maxent), " (",
          sum(train_maxent$presence == 1), " Presence / ",
          sum(train_maxent$presence == 0), " Background)")
  
  # -------------------------------------------------------------------
  # 6c. BRT/RF: 10 Runs, jeweils Train-Presence + Run-eigene Pseudo-Absences
  # -------------------------------------------------------------------
  train_brt_runs <- lapply(sampling$runs, function(run) {
    
    # Presence-Zeilen des Runs, die auch im Train-Split liegen
    run_pres <- run[run$presence == 1, ]
    run_pres_key <- paste(run_pres$x, run_pres$y)
    run_pres_train <- run_pres[run_pres_key %in% train_pres_key, ]
    
    # Pseudo-Absences des Runs (alle, unverändert)
    run_abs <- run[run$presence == 0, ]
    
    run_train        <- rbind(run_pres_train, run_abs)
    run_train$split  <- "train"
    run_train
  })
  
  message("BRT/RF Train: ", length(train_brt_runs), " Runs, je ~",
          round(mean(sapply(train_brt_runs, nrow))), " Punkte im Schnitt")
  
  # -------------------------------------------------------------------
  # 7. True Absences → gemeinsamer Testdatensatz für alle Modelle
  # -------------------------------------------------------------------
  pa_raster <- species_PA$pa.raster
  
  # terra-Objekte werden beim saveRDS() "gepackt" -> vor Nutzung entpacken
  if (inherits(pa_raster, "PackedSpatRaster")) {
    pa_raster <- terra::unwrap(pa_raster)
  }
  
  set.seed(seed)
  true_abs_pts <- terra::spatSample(
    pa_raster,
    size      = n_true_abs * 2,  # Mehr samplen, dann filtern
    method    = "random",
    na.rm     = TRUE,
    as.points = TRUE
  ) |> as.data.frame(geom = "XY")
  
  # Nur echte Abwesenheiten (pa == 0)
  true_abs <- true_abs_pts[true_abs_pts[, 1] == 0, c("x", "y")]
  true_abs <- true_abs[1:min(n_true_abs, nrow(true_abs)), ]
  true_abs$presence <- 0
  
  message("True Absences (Evaluation): ", nrow(true_abs))
  
  # -------------------------------------------------------------------
  # 8. Gemeinsamer Testdatensatz (für alle Modelle gleich)
  # -------------------------------------------------------------------
  test_data        <- rbind(test_pres, true_abs)
  test_data$split  <- "test"
  
  message("─────────────────────────────────────")
  message("TEST (gemeinsam): ", nrow(test_data), " Punkte | ",
          sum(test_data$presence == 1), " Presence | ",
          sum(test_data$presence == 0), " True Absence")
  message("─────────────────────────────────────")
  
  list(
    train_glm      = train_glm,
    train_maxent   = train_maxent,
    train_brt_runs = train_brt_runs,   # Liste mit 10 Trainingssets
    test           = test_data,        # gemeinsamer Testdatensatz
    knn            = knn,
    fold_ids       = fold_ids
  )
}

library(sf)

# =============================================================================
# 1.2 Hilfsfunktion: Datensätze als GeoPackage für QGIS speichern
# =============================================================================

save_as_gpkg <- function(split, sp, out_dir = "Data/species/split_knndm/gpkg") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gpkg_path <- file.path(out_dir, paste0("split_", sp, ".gpkg"))
  
  # -------------------------------------------------------------------
  # GLM Training
  # -------------------------------------------------------------------
  glm_sf <- sf::st_as_sf(split$train_glm, coords = c("x", "y"), crs = 3035)
  glm_sf$type <- ifelse(glm_sf$presence == 1, "presence", "pseudo_absence")
  sf::st_write(glm_sf, gpkg_path, layer = "train_glm",
               delete_layer = TRUE, quiet = TRUE)
  
  # -------------------------------------------------------------------
  # MaxEnt Training
  # -------------------------------------------------------------------
  maxent_sf <- sf::st_as_sf(split$train_maxent, coords = c("x", "y"), crs = 3035)
  maxent_sf$type <- ifelse(maxent_sf$presence == 1, "presence", "background")
  sf::st_write(maxent_sf, gpkg_path, layer = "train_maxent",
               delete_layer = TRUE, quiet = TRUE)
  
  # -------------------------------------------------------------------
  # BRT/RF Training (alle 10 Runs in einem Layer, mit Run-ID als Spalte)
  # -------------------------------------------------------------------
  brt_combined <- do.call(rbind, lapply(seq_along(split$train_brt_runs), function(i) {
    run_df <- split$train_brt_runs[[i]]
    run_df$run <- i
    run_df
  }))
  brt_sf <- sf::st_as_sf(brt_combined, coords = c("x", "y"), crs = 3035)
  brt_sf$type <- ifelse(brt_sf$presence == 1, "presence", "pseudo_absence")
  sf::st_write(brt_sf, gpkg_path, layer = "train_brt_runs",
               delete_layer = TRUE, quiet = TRUE)
  
  # -------------------------------------------------------------------
  # Testdatensatz (gemeinsam)
  # -------------------------------------------------------------------
  test_sf <- sf::st_as_sf(split$test, coords = c("x", "y"), crs = 3035)
  test_sf$type <- ifelse(test_sf$presence == 1, "presence", "true_absence")
  sf::st_write(test_sf, gpkg_path, layer = "test",
               delete_layer = TRUE, quiet = TRUE)
  
  # -------------------------------------------------------------------
  # Optional: kNNDM Folds visualisieren (alle Presence-Punkte mit Fold-ID)
  # -------------------------------------------------------------------
  pres_folds <- sampling$po$sample.points[
    sampling$po$sample.points$Observed == TRUE, ]
  pres_folds <- data.frame(x = pres_folds$x, y = pres_folds$y,
                           fold = split$fold_ids)
  folds_sf <- sf::st_as_sf(pres_folds, coords = c("x", "y"), crs = 3035)
  folds_sf$fold <- as.factor(folds_sf$fold)
  sf::st_write(folds_sf, gpkg_path, layer = "knndm_folds",
               delete_layer = TRUE, quiet = TRUE)
  
  message("✓ GeoPackage gespeichert: ", gpkg_path)
}

# =============================================================================
# 2. Split für alle Arten + Speichern
# =============================================================================

species_list <- c("narrow", "low_mid", "high_mid", "broad")

dir.create("Data/species/split_knndm", recursive = TRUE, showWarnings = FALSE)

for (sp in species_list) {
  
  message("===== ", sp, " =====")
  
  sampling   <- readRDS(paste0("Data/species/sampling_", sp, ".RDS"))
  species_PA <- readRDS(paste0("Data/species/species_",  sp, ".RDS"))$virtual_PA
  
  split <- split_data_knndm(
    sampling    = sampling,
    species_PA  = species_PA,
    train_ratio = 0.8,
    k           = 5,
    seed        = 42,
    n_true_abs  = 1000
  )
  
  saveRDS(split$train_glm,      paste0("Data/species/split_knndm/train_glm_",    sp, ".RDS"))
  saveRDS(split$train_maxent,   paste0("Data/species/split_knndm/train_maxent_", sp, ".RDS"))
  saveRDS(split$train_brt_runs, paste0("Data/species/split_knndm/train_brt_",    sp, ".RDS"))
  saveRDS(split$test,           paste0("Data/species/split_knndm/test_",         sp, ".RDS"))
  
  save_as_gpkg(split, sp)
  
  
}

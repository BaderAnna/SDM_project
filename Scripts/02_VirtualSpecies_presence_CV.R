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
  
  # -------------------------------------------------------------------
  # 6. Pseudo-Absences → alle ins Training
  # -------------------------------------------------------------------
  bg <- as.data.frame(sampling$bg_glm)
  bg <- data.frame(x = bg$x, y = bg$y, presence = 0)
  
  message("Pseudo-Absences (Training): ", nrow(bg))
  
  # -------------------------------------------------------------------
  # 7. True Absences → alle in die Evaluation
  # -------------------------------------------------------------------
  pa_raster <- species_PA$pa.raster
  
  # unwrap if needed
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
  # 8. Trainings- und Testdaten zusammenfügen
  # -------------------------------------------------------------------
  
  # Trainingsdaten: 80% Presence + alle Pseudo-Absences
  train_data          <- rbind(train_pres, bg)
  train_data$split    <- "train"
  
  # Testdaten: 20% Presence + alle True Absences
  test_data           <- rbind(test_pres, true_abs)
  test_data$split     <- "test"
  
  message("─────────────────────────────────────")
  message("TRAIN: ", nrow(train_data), " Punkte | ",
          sum(train_data$presence == 1), " Presence | ",
          sum(train_data$presence == 0), " Pseudo-Absence")
  message("TEST:  ", nrow(test_data),  " Punkte | ",
          sum(test_data$presence == 1),  " Presence | ",
          sum(test_data$presence == 0),  " True Absence")
  message("─────────────────────────────────────")
  
  list(
    train    = train_data,
    test     = test_data,
    knn      = knn,
    fold_ids = fold_ids
  )
}

# =============================================================================
# 2. Split für alle Arten
# =============================================================================

species_list <- c("narrow", "low_mid", "high_mid", "broad")

dir.create("Data/species/split_knndm", recursive = TRUE, showWarnings = FALSE)

for (sp in species_list) {
  
  message("===== ", sp, " =====")
  
  # Lade Daten
  sampling   <- readRDS(paste0("Data/species/sampling_", sp, ".RDS"))
  species_PA <- readRDS(paste0("Data/species/species_",  sp, ".RDS"))$virtual_PA
  
  # Split
  split <- split_data_knndm(
    sampling    = sampling,
    species_PA  = species_PA,
    train_ratio = 0.8,
    k           = 5,
    seed        = 42,
    n_true_abs  = 1000
  )
  
  # Speichern
  saveRDS(split$train, paste0("Data/species/split_knndm/train_", sp, ".RDS"))
  saveRDS(split$test,  paste0("Data/species/split_knndm/test_",  sp, ".RDS"))
  
}

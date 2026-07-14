# =============================================================================
# DATEN-SPLITTING MIT kNNDM (80/20)
# Test: Nur für narrow
# =============================================================================

library(terra)
library(sf)
library(CAST)

# =============================================================================
# 1. Sampling-Daten laden
# =============================================================================

sampling <- readRDS("Data/species/sampling_narrow.RDS")

# Prüfe die Spaltennamen
pres <- sampling$po$sample.points[sampling$po$sample.points$Observed == TRUE, ]
bg   <- sampling$bg_glm

message("Spaltennamen pres: ", paste(names(pres), collapse = ", "))
message("Spaltennamen bg: ",   paste(names(bg),   collapse = ", "))

# =============================================================================
# 2. kNNDM-basiertes Split (80/20)
# =============================================================================

split_data_knndm <- function(sampling, train_ratio = 0.8, k = 5, seed = 42) {
  
  # Präsenzpunkte
  pres <- sampling$po$sample.points[sampling$po$sample.points$Observed == TRUE, ]
  pres <- data.frame(x = pres$x, y = pres$y, presence = 1)
  
  # Hintergrundpunkte
  bg <- as.data.frame(sampling$bg_glm)
  bg <- data.frame(x = bg$x, y = bg$y, presence = 0)
  
  # Alle Punkte zusammenfügen
  all_points <- rbind(pres, bg)
  
  message("Gesamt: ", nrow(all_points), " Punkte | ",
          sum(all_points$presence == 1), " Präsenzen | ",
          sum(all_points$presence == 0), " Hintergrund")
  
  # Räumliche Geometrie: in EPSG:3035
  all_sf <- sf::st_as_sf(all_points, coords = c("x", "y"), crs = 3035)
  
  # Vorhersagepunkte (für kNNDM)
  predictor_raster <- terra::rast("Data/raster/pc_raster_masked.tif")
  predictor_raster <- terra::project(predictor_raster, "EPSG:3035")
  
  pred_pts <- terra::spatSample(
    predictor_raster,
    size = 10000,
    method = "random",
    as.points = TRUE,
    na.rm = TRUE
  ) |> sf::st_as_sf()
  
  sf::st_crs(pred_pts) <- 3035
  pred_pts <- sf::st_transform(pred_pts, crs = 3035)
  
  # kNNDM auf ALLE Punkte anwenden
  set.seed(seed)
  knn <- tryCatch(
    CAST::knndm(
      tpoints    = all_sf,
      predpoints = pred_pts,
      k          = k
    ),
    error = function(e) {
      message("kNNDM fehlgeschlagen: ", e$message)
      stop("kNNDM konnte nicht durchgeführt werden.")
    }
  )
  
  # Fold-Nummern für alle Punkte
  fold_ids <- knn$clusters
  all_points$fold <- fold_ids
  
  # Folds aufteilen
  n_folds <- max(fold_ids)
  set.seed(seed)
  train_folds <- sample(1:n_folds, size = floor(train_ratio * n_folds))
  test_folds  <- setdiff(1:n_folds, train_folds)
  
  # Trainings- und Testdaten (Präsenz + Hintergrund)
  train_data <- all_points[all_points$fold %in% train_folds, ]
  test_data  <- all_points[all_points$fold %in% test_folds,  ]
  
  train_data$split <- "train"
  test_data$split  <- "test"
  
  message("Train: ", nrow(train_data), " Punkte | ",
          sum(train_data$presence == 1), " Präsenzen | ",
          sum(train_data$presence == 0), " Hintergrund")
  
  message("Test:  ", nrow(test_data), " Punkte | ",
          sum(test_data$presence == 1), " Präsenzen | ",
          sum(test_data$presence == 0), " Hintergrund")
  
  list(
    train = train_data,
    test  = test_data
  )
}

# =============================================================================
# 3. Split durchführen (nur narrow zum Testen)
# =============================================================================

sampling <- readRDS("Data/species/sampling_narrow.RDS")

dir.create("Data/species/split_knndm", recursive = TRUE, showWarnings = FALSE)

split_narrow <- split_data_knndm(sampling, train_ratio = 0.8, k = 5, seed = 42)

message("✅ narrow gesplittet: ",
        nrow(split_narrow$train), " Trainings-, ",
        nrow(split_narrow$test),  " Testdaten")


# =============================================================================
# 4. Speichern
# =============================================================================

saveRDS(split_narrow$train, "Data/species/split_knndm/train_narrow.RDS")
saveRDS(split_narrow$test,  "Data/species/split_knndm/test_narrow.RDS")

message("✅ Gespeichert!")

























# =============================================================================
# DATEN-SPLITTING MIT kNNDM (80/20)
# Nur für Evaluation – keine Modelle
# =============================================================================

library(terra)
library(sf)
library(CAST)  # Für kNNDM

# =============================================================================
# 1. Funktion: kNNDM-basiertes Split (80/20)
# =============================================================================

split_data_knndm <- function(sampling, train_ratio = 0.8, k = 5, seed = 42) {
  
  # Präsenzpunkte
  pres <- sampling$po$sample.points[sampling$po$sample.points$Observed == TRUE, ]
  
  # Räumliche Geometrie: in EPSG:3035
  pres_sf <- sf::st_as_sf(pres, coords = c("x", "y"), crs = 3035)
  
  # Vorhersagepunkte (für kNNDM)
  predictor_raster <- terra::rast("Data/raster/pc_raster_masked.tif")
  predictor_raster <- terra::project(predictor_raster, "EPSG:3035")
  
  pred_pts <- terra::spatSample(
    predictor_raster,
    size = 10000,
    method = "random",
    as.points = TRUE,
    na.rm = TRUE
  ) |> sf::st_as_sf()
  
  sf::st_crs(pred_pts) <- 3035
  pred_pts <- sf::st_transform(pred_pts, crs = 3035)
  
  # kNNDM
  set.seed(seed)
  knn <- tryCatch(
    CAST::knndm(
      tpoints = pres_sf,
      predpoints = pred_pts,
      k = k
    ),
    error = function(e) {
      message("kNNDM fehlgeschlagen: ", e$message)
      stop("kNNDM konnte nicht durchgeführt werden.")
    }
  )
  
  fold_ids <- knn$clusters
  n_folds <- max(fold_ids)
  
  set.seed(seed)
  train_folds <- sample(1:n_folds, size = floor(train_ratio * n_folds))
  test_folds  <- setdiff(1:n_folds, train_folds)
  
  train_pres <- pres[fold_ids %in% train_folds, ]
  test_pres  <- pres[fold_ids %in% test_folds,  ]
  
  # Hintergrundpunkte
  bg <- sampling$bg_glm
  
  # 🔥 Wichtig: Prüfe die Spaltennamen von bg
  print(names(bg))  # ✅ Zeige die Spaltennamen an
  
  # 🔥 Finde die richtigen Spaltennamen für x, y
  # Beispiel: Wenn Spaltennamen sind: x, y
  # Dann: bg_subset <- bg[, c("x", "y")]
  
  # ✅ Ersetze "X", "Y" durch die tatsächlichen Spaltennamen
  # Beispiel: Wenn Spaltennamen sind: x, y
  bg_subset <- bg[, c("x", "y")]  # ✅ Ersetze durch die richtigen Namen!
  
  # 🔥 Stelle sicher, dass die Spaltennamen korrekt sind
  names(bg_subset) <- c("x", "y")  # ✅ Setze die Namen auf x, y
  
  bg_subset$presence <- 0
  bg_subset$ID <- 1:nrow(bg_subset)
  
  # 🔥 Jetzt rbind() funktioniert!
  train_data <- rbind(train_pres, bg_subset)
  train_data$split <- "train"
  
  test_data <- test_pres
  test_data$split <- "test"
  
  list(
    train = train_data,
    test  = test_data,
    train_pres = train_pres,
    test_pres  = test_pres
  )
}

# =============================================================================
# 2. Lade alle Sampling-Daten und splittete sie
# =============================================================================

species_list <- c("narrow", "low_mid", "high_mid", "broad")

# Erstelle Ordner für gesplittete Daten
dir.create("Data/species/split_knndm", recursive = TRUE, showWarnings = FALSE)

for (sp in species_list) {
  
  # Lade Sampling-Daten
  sampling <- readRDS(paste0("Data/species/sampling_", sp, ".RDS"))
  
  # Split mit kNNDM
  split <- split_data_knndm(sampling, train_ratio = 0.8, k = 5, seed = 42)
  
  # Speichern
  saveRDS(split$train, paste0("Data/species/split_knndm/train_", sp, ".RDS"))
  saveRDS(split$test,  paste0("Data/species/split_knndm/test_", sp, ".RDS"))
  
  # Nachricht
  message("✓ ", sp, " mit kNNDM gesplittet: ",
          nrow(split$train), " Trainings-, ",
          nrow(split$test), " Testdaten")
}


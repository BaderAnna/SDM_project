# =============================================================================
# GLM MODELLIERUNG (TRAINING)
# Option B: alle 19 WorldClim Variablen -> VIF-Reduktion
# Nutzt die im Skript "02_VirtualSpecies_presence_CV.R" (split_data_knndm())
# bereits erzeugten Trainingsdaten (train_glm_<sp>.RDS).
# Evaluierung erfolgt in separatem Skript auf Basis von test_<sp>.RDS.
# =============================================================================

library(terra)
library(sf)

# =============================================================================
# 1 - Raster + Variablen laden
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars      <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - GLM-Formel dynamisch bauen
# =============================================================================

build_glm_formula <- function(vars) {
  linear_terms    <- paste(vars, collapse = " + ")
  quadratic_terms <- paste0("I(", vars, "^2)", collapse = " + ")
  as.formula(paste("presence ~", linear_terms, "+", quadratic_terms))
}

glm_formula <- build_glm_formula(selected_vars)
print(glm_formula)

# =============================================================================
# 3 - Trainingsdaten vorbereiten
#     (Fold-Split kommt bereits aus train_glm_<sp>.RDS -> hier nur noch
#      Umweltwerte an den Koordinaten extrahieren)
# =============================================================================

prepare_glm_train_data <- function(train_df, env_raster) {
  
  if (nrow(train_df) == 0) stop("train_df ist leer!")
  if (any(is.na(train_df$x)) || any(is.na(train_df$y))) {
    stop("x oder y enthalten NA-Werte!")
  }
  
  coords_3035  <- as.matrix(train_df[, c("x", "y")])
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords_3035))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  dat <- cbind(train_df, env_vals)
  dat <- na.omit(dat)
  
  message("prepare_glm_train_data: ", nrow(dat), " Zeilen | ",
          sum(dat$presence == 1), " Presence | ",
          sum(dat$presence == 0), " Pseudo-Absence")
  
  return(dat)
}

# =============================================================================
# 4 - GLM auf Trainingsdaten fitten
# =============================================================================

run_glm_train <- function(train_df, env_raster, glm_formula, seed = 42) {
  
  dat <- prepare_glm_train_data(train_df, env_raster)
  
  if (nrow(dat) == 0) stop("dat ist nach prepare_glm_train_data() leer!")
  
  # Prüfe, ob alle Variablen der Formel in dat existieren
  missing_vars <- setdiff(all.vars(glm_formula), names(dat))
  if (length(missing_vars) > 0) {
    stop("Fehlende Variablen in dat: ", paste(missing_vars, collapse = ", "))
  }
  
  set.seed(seed)
  
  # Gewichte: Presences = 1, Pseudo-Absences = n_pres / n_abs
  n_pres <- sum(dat$presence == 1)
  n_abs  <- sum(dat$presence == 0)
  dat$.weights <- ifelse(dat$presence == 1, 1, n_pres / n_abs)
  
  model <- withCallingHandlers(
    glm(glm_formula, data = dat, family = binomial, weights = .weights),
    error = function(e) {
      message("Modelltraining fehlgeschlagen: ", e$message)
      stop("Modell konnte nicht trainiert werden.")
    },
    warning = function(w) {
      message("Warnung beim Modelltraining: ", w$message)
    }
  )
  
  prediction <- terra::predict(env_raster, model, type = "response")
  
  list(model = model, prediction = prediction, data = dat)
}

# =============================================================================
# 5 - Für alle Arten ausführen
# =============================================================================

species_list     <- c("narrow", "low_mid", "high_mid", "broad")
glm_models_train <- list()

for (sp in species_list) {
  
  message("===== ", sp, " (Training) =====")
  
  train_glm <- readRDS(paste0("Data/species/split_knndm/train_glm_", sp, ".RDS"))
  
  glm_models_train[[sp]] <- run_glm_train(
    train_df    = train_glm,
    env_raster  = env_raster_masked,
    glm_formula = glm_formula
  )
}

# =============================================================================
# 6 - Speichern der Trainingsmodelle
# =============================================================================

dir.create("Data/models/knndm", recursive = TRUE, showWarnings = FALSE)

for (sp in species_list) {
  saveRDS(
    glm_models_train[[sp]],
    paste0("Data/models/knndm/glm_", sp, "_train.RDS")
  )
  message("✓ Trainingsmodell gespeichert: glm_", sp, "_train.RDS")
}





'# =============================================================================
# GLM MODELLIERUNG
# Option B: alle 19 WorldClim Variablen -> VIF-Reduktion
# Reines Modelltraining (Evaluierung erfolgt in separatem Skript)
# =============================================================================

library(terra)
library(sf)

# =============================================================================
# 1 - Raster + Variablen laden
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars     <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 1. Lade Sampling-Daten
# =============================================================================

sampling_narrow   <- readRDS("Data/species/sampling_narrow.RDS")
sampling_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
sampling_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
sampling_broad    <- readRDS("Data/species/sampling_broad.RDS")

# =============================================================================
# 2. Erstelle Fold-IDs für alle Arten
# =============================================================================

species_list <- c("narrow", "low_mid", "high_mid", "broad")

for (sp in species_list) {
  
  # Lade Sampling-Daten
  sampling <- get(paste0("sampling_", sp))
  
  # Präsenzpunkte
  pres <- sampling$po$sample.points[sampling$po$sample.points$Observed == TRUE, ]
  
  # Räumliche Geometrie
  pres_sf <- sf::st_as_sf(pres, coords = c("x", "y"), crs = 3035)
  
  # Vorhersagepunkte
  predictor_raster <- terra::rast("Data/raster/pc_raster_masked.tif")
  pred_pts <- terra::spatSample(predictor_raster, size = 10000, method = "random", as.points = TRUE, na.rm = TRUE) |> sf::st_as_sf()
  
  # kNNDM
  knn <- CAST::knndm(
    tpoints = pres_sf,
    predpoints = pred_pts,
    k = 5
  )
  
  # Fold-IDs
  fold_ids <- knn$clusters
  
  # Speichern
  saveRDS(fold_ids, paste0("Data/species/split_knndm/fold_ids_", sp, ".RDS"))
  
  message("✓ Fold-IDs für ", sp, " gespeichert")
}

# =============================================================================
# 2 - GLM-Formel dynamisch bauen
# =============================================================================

build_glm_formula <- function(vars) {
  linear_terms    <- paste(vars, collapse = " + ")
  quadratic_terms <- paste0("I(", vars, "^2)", collapse = " + ")
  as.formula(paste("presence ~", linear_terms, "+", quadratic_terms))
}

glm_formula <- build_glm_formula(selected_vars)
print(glm_formula)

# =============================================================================
# 3 - Daten vorbereiten
# =============================================================================

prepare_glm_data <- function(sampling, env_raster) {
  
  sp   <- sampling$po$sample.points
  pres <- sp[sp$Observed == TRUE, c("x", "y")]
  pres$presence <- 1
  
  abs <- as.data.frame(sampling$bg_glm)
  abs$presence <- 0
  
  dat <- rbind(pres, abs)
  dat <- cbind(dat, ID = 1:nrow(dat))
  
  if (nrow(dat) == 0) stop("dat ist leer!")
  if (any(is.na(dat$x)) || any(is.na(dat$y))) stop("x oder y enthalten NA-Werte!")
  
  coords_3035 <- as.matrix(dat[, c("x", "y")])
  
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords_3035))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  dat <- cbind(dat, env_vals)
  dat <- na.omit(dat)
  
  message("prepare_glm_data: ", nrow(dat), " Zeilen | ",
          sum(dat$presence == 1), " Präsenzen | ",
          sum(dat$presence == 0), " Hintergrund")
  
  return(dat)
}

# =============================================================================
# 4 - knndm Fold-Erstellung
# =============================================================================

run_knndm_folds <- function(dat_sf_geom, predictor_raster, train_ratio = 0.8,
                            k = 5, seed = 42, n_predpoints = 10000) { # gleiche Anzahl wie background points 
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
# 5 - GLM mit knndm: Nur finale Modelle (mit vorgegebenen Fold-Nummern)
# =============================================================================

run_glm_knndm_final <- function(sampling, env_raster, glm_formula, seed = 42) {
  
  dat <- prepare_glm_data(sampling, env_raster)
  
  if (nrow(dat) == 0) stop("dat ist nach prepare_glm_data() leer!")
  
  set.seed(seed)
  
  # 1. Lade die Fold-Nummern aus dem separaten Skript
  # Annahme: Die Fold-Nummern sind in "Data/species/split_knndm/fold_ids_narrow.RDS"
  fold_file <- paste0("Data/species/split_knndm/fold_ids_", 
                      ifelse(sampling$po$sample.points$species == "narrow", "narrow",
                             ifelse(sampling$po$sample.points$species == "low_mid", "low_mid",
                                    ifelse(sampling$po$sample.points$species == "high_mid", "high_mid", "broad"))),
                      ".RDS")
  
  if (!file.exists(fold_file)) {
    stop("Fold-IDs nicht gefunden: ", fold_file)
  }
  
  fold_ids <- readRDS(fold_file)
  
  # 2. Prüfe, ob fold_ids zu dat passt
  if (length(fold_ids) != nrow(dat)) {
    stop("Anzahl Fold-IDs stimmt nicht mit Daten überein.")
  }
  
  # 3. Füge fold_ids zu dat hinzu
  dat$fold <- fold_ids
  
  # 4. Prüfe, ob alle Variablen in der Formel existieren
  missing_vars <- setdiff(all.vars(glm_formula), names(dat))
  if (length(missing_vars) > 0) {
    stop("Fehlende Variablen in dat: ", paste(missing_vars, collapse = ", "))
  }
  
  # 5. Gewichte setzen
  n_pres <- sum(dat$presence == 1)
  n_abs  <- sum(dat$presence == 0)
  dat$.weights <- ifelse(dat$presence == 1, 1, n_pres / n_abs)
  
  # 6. Modell trainieren
  final_model <- tryCatch(
    glm(glm_formula, data = dat, family = binomial, weights = .weights), 
    error = function(e) {
      message("Modelltraining fehlgeschlagen: ", e$message)
      stop("Modell konnte nicht trainiert werden.")
    }, 
    warning = function(w) {
      message("Warnung beim Modelltraining: ", w$message)
    }
  )
  
  # 7. Vorhersage
  final_pred <- terra::predict(env_raster, final_model, type = "response")
  
  # 8. Rückgabe
  list(model = final_model, prediction = final_pred, data = dat, fold_ids = fold_ids)
}
'
'# =============================================================================
# 5 - GLM mit knndm: Nur finale Modelle 
# =============================================================================

run_glm_knndm_final <- function(sampling, env_raster, glm_formula, seed = 42) {
  
  dat <- prepare_glm_data(sampling, env_raster)
  
  if (nrow(dat) == 0) stop("dat ist nach prepare_glm_data() leer!")
  
  set.seed(seed)
  
  missing_vars <- setdiff(all.vars(glm_formula), names(dat))
  if (length(missing_vars) > 0) {
    stop("Fehlende Variablen in dat: ", paste(missing_vars, collapse = ", "))
  }
  
  # Erstelle räumliche Geometrie
  dat_sf <- sf::st_as_sf(dat, coords = c("x", "y"), crs = 3035)
  dat_sf_geom <- sf::st_as_sf(
    data.frame(geometry = sf::st_geometry(dat_sf)),
    crs = sf::st_crs(dat_sf)
  )
  
  #  kNNDM: Finde beste Fold-Aufteilung
  knn <- tryCatch(
    run_knndm_folds(
      dat_sf_geom,
      predictor_raster = env_raster,
      k = 5,
      seed = seed,
      n_predpoints = 10000
    ),
    error = function(e) {
      message("kNNDM fehlgeschlagen: ", e$message)
      stop("kNNDM konnte nicht durchgeführt werden.")
    }
  )
  
  #  Prüfe, ob kNNDM erfolgreich war
  if (is.null(knn)) {
    stop("kNNDM hat keine gültigen Fold-Cluster zurückgegeben.")
  }
  
  n_pres <- sum(dat$presence == 1)
  n_abs  <- sum(dat$presence == 0)
  dat$.weights <- ifelse(dat$presence == 1, 1, n_pres / n_abs)
  
  final_model <- tryCatch(
    glm(glm_formula, data = dat, family = binomial, weights = .weights), 
    error = function(e) {
      message("Modelltraining fehlgeschlagen: ", e$message)
      stop("Modell konnte nicht trainiert werden.")
    }, 
    warning = function(w) {
      message("Warnung beim Modelltraining: ", w$message)
    }
  )
  
  final_pred <- terra::predict(env_raster, final_model, type = "response")
  
  list(model = final_model, prediction = final_pred, data = dat)
}'

# =============================================================================
# 6 - Daten laden + GLM mit kNNDM ausführen (nur finale Modelle)
# =============================================================================

sampling_narrow   <- readRDS("Data/species/sampling_narrow.RDS")
sampling_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
sampling_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
sampling_broad    <- readRDS("Data/species/sampling_broad.RDS")

message("===== Narrow (kNNDM) =====")
glm_narrow_knndm <- run_glm_knndm_final(
  sampling = sampling_narrow,
  env_raster = env_raster_masked,
  glm_formula = glm_formula
)

message("===== Low-Mid (kNNDM) =====")
glm_low_mid_knndm <- run_glm_knndm_final(
  sampling = sampling_low_mid,
  env_raster = env_raster_masked,
  glm_formula = glm_formula
)

message("===== High-Mid (kNNDM) =====")
glm_high_mid_knndm <- run_glm_knndm_final(
  sampling = sampling_high_mid,
  env_raster = env_raster_masked,
  glm_formula = glm_formula
)

message("===== Broad (kNNDM) =====")
glm_broad_knndm <- run_glm_knndm_final(
  sampling = sampling_broad,
  env_raster = env_raster_masked,
  glm_formula = glm_formula
)

# =============================================================================
# 7 - Speichern der finalen Modelle (mit kNNDM)
# =============================================================================

dir.create("Data/models/knndm", recursive = TRUE, showWarnings = FALSE)

saveRDS(glm_narrow_knndm,   "Data/models/knndm/glm_narrow_knndm.RDS")
saveRDS(glm_low_mid_knndm,  "Data/models/knndm/glm_low_mid_knndm.RDS")
saveRDS(glm_high_mid_knndm, "Data/models/knndm/glm_high_mid_knndm.RDS")
saveRDS(glm_broad_knndm,    "Data/models/knndm/glm_broad_knndm.RDS")

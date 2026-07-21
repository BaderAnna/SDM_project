# =============================================================================
# MAXNET MODELLIERUNG (TRAINING)
# WorldClim Variablen nach VIF-Reduktion
# Nutzt die im Skript "02_VirtualSpecies_presence_CV.R" (split_data_knndm())
# bereits erzeugten Trainingsdaten (train_maxent_<sp>.RDS).
# Evaluierung erfolgt in separatem Skript auf Basis von test_<sp>.RDS.
# =============================================================================

library(terra)
library(sf)
library(maxnet)

# =============================================================================
# 1 - Raster + Variablen laden
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars      <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - Trainingsdaten vorbereiten
#     (Fold-Split kommt bereits aus train_maxent_<sp>.RDS -> hier nur noch
#      Umweltwerte an den Koordinaten extrahieren + skalieren)
# =============================================================================

prepare_maxnet_train_data <- function(train_df, env_raster) {
  
  if (nrow(train_df) == 0) stop("train_df ist leer!")
  if (any(is.na(train_df$x)) || any(is.na(train_df$y))) {
    stop("x oder y enthalten NA-Werte!")
  }
  
  coords_3035  <- as.matrix(train_df[, c("x", "y")])
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords_3035))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  message("NA-Anteil env_vals: ",
          round(mean(is.na(env_vals)) * 100, 1), "%")
  
  # Variablen standardisieren (z-Transformation)
  env_vals_scaled <- as.data.frame(scale(env_vals))
  
  # Skalierungsparameter speichern (für spätere Vorhersage!)
  scale_center <- attr(scale(env_vals), "scaled:center")
  scale_scale  <- attr(scale(env_vals), "scaled:scale")
  
  dat <- cbind(train_df, env_vals_scaled)
  dat <- na.omit(dat)
  
  message("prepare_maxnet_train_data: ", nrow(dat), " Zeilen | ",
          sum(dat$presence == 1), " Presence | ",
          sum(dat$presence == 0), " Pseudo-Absence/Background")
  
  # Skalierungsparameter zurückgeben
  attr(dat, "scale_center") <- scale_center
  attr(dat, "scale_scale")  <- scale_scale
  
  return(dat)
}

# =============================================================================
# 3 - MaxNet auf Trainingsdaten fitten
# =============================================================================

run_maxnet_train <- function(train_df, env_raster, seed = 42) {
  
  dat <- prepare_maxnet_train_data(train_df, env_raster)
  
  if (nrow(dat) == 0) stop("dat ist nach prepare_maxnet_train_data() leer!")
  
  env_cols <- names(env_raster)
  
  # Prüfe, ob alle Umweltvariablen in dat existieren
  missing_vars <- setdiff(env_cols, names(dat))
  if (length(missing_vars) > 0) {
    stop("Fehlende Variablen in dat: ", paste(missing_vars, collapse = ", "))
  }
  
  set.seed(seed)
  
  scale_center <- attr(dat, "scale_center")
  scale_scale  <- attr(dat, "scale_scale")
  
  all_env <- dat[, env_cols, drop = FALSE]
  
  model <- withCallingHandlers(
    maxnet(
      p    = dat$presence,
      data = all_env,
      f    = maxnet.formula(dat$presence, all_env, classes = "default")
    ),
    error = function(e) {
      message("Modelltraining fehlgeschlagen: ", e$message)
      stop("Modell konnte nicht trainiert werden.")
    },
    warning = function(w) {
      message("Warnung beim Modelltraining: ", w$message)
    }
  )
  
  # Raster für die Vorhersage passend zu den Trainingsdaten skalieren
  env_raster_scaled <- (env_raster - scale_center) / scale_scale
  names(env_raster_scaled) <- env_cols
  
  prediction <- terra::predict(
    env_raster_scaled,
    model,
    fun    = function(model, data) predict(model, data, type = "cloglog"),
    na.rm  = TRUE
  )
  
  list(
    model        = model,
    prediction   = prediction,
    data         = dat,
    scale_center = scale_center,
    scale_scale  = scale_scale
  )
}

# =============================================================================
# 4 - Für alle Arten ausführen
# =============================================================================

species_list         <- c("narrow", "low_mid", "high_mid", "broad")
maxnet_models_train   <- list()

for (sp in species_list) {
  
  message("===== ", sp, " (Training) =====")
  
  train_maxent <- readRDS(paste0("Data/species/split_knndm/train_maxent_", sp, ".RDS"))
  
  maxnet_models_train[[sp]] <- run_maxnet_train(
    train_df   = train_maxent,
    env_raster = env_raster_masked
  )
}

# =============================================================================
# 5 - Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
for (sp in species_list) {
  plot(maxnet_models_train[[sp]]$prediction, main = paste("MaxNet", sp))
}
par(mfrow = c(1, 1))

# =============================================================================
# 6 - Speichern der Trainingsmodelle
# =============================================================================

dir.create("Data/models/knndm", recursive = TRUE, showWarnings = FALSE)

for (sp in species_list) {
  saveRDS(
    maxnet_models_train[[sp]],
    paste0("Data/models/knndm/maxnet_", sp, "_train.RDS")
  )
  message("✓ Trainingsmodell gespeichert: maxnet_", sp, "_train.RDS")
}













'# =============================================================================
# MAXNET MODELLIERUNG
# WorldClim Variablen nach VIF-Reduktion
# Reines Modelltraining (Evaluierung erfolgt in separatem Skript)
# =============================================================================

library(terra)
library(sf)
library(maxnet)

# =============================================================================
# 1 - Raster + Variablen laden
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars     <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - Daten vorbereiten
#     WICHTIG: Koordinaten bereits EPSG:3035 -> keine Transformation!
#     bg_maxent aus sampling verwenden (nicht neu samplen!)
# =============================================================================

prepare_maxnet_data <- function(sampling, env_raster) {
  
  # Präsenzpunkte
  sp   <- sampling$po$sample.points
  pres <- sp[sp$Observed == TRUE, c("x", "y")]
  pres$presence <- 1
  
  # Background
  bg <- as.data.frame(sampling$bg_maxent)
  bg$presence <- 0
  
  # Kombinieren
  dat <- rbind(pres, bg)
  
  if (nrow(dat) == 0) stop("dat ist leer!")
  if (any(is.na(dat$x)) || any(is.na(dat$y))) stop("x oder y enthalten NA!")
  
  coords <- as.matrix(dat[, c("x", "y")])
  
  # Umweltwerte extrahieren
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  message("NA-Anteil env_vals: ",
          round(mean(is.na(env_vals)) * 100, 1), "%")
  
  # Variablen standardisieren (z-Transformation)
  env_vals_scaled <- as.data.frame(scale(env_vals))
  
  # Skalierungsparameter speichern (für spätere Vorhersage!)
  scale_center <- attr(scale(env_vals), "scaled:center")
  scale_scale  <- attr(scale(env_vals), "scaled:scale")
  
  message("Skalierung:")
  message("  Center: ", paste(round(scale_center, 2), collapse = ", "))
  message("  Scale:  ", paste(round(scale_scale,  2), collapse = ", "))
  
  dat <- cbind(dat, env_vals_scaled)
  dat <- na.omit(dat)
  
  message("prepare_maxnet_data: ", nrow(dat), " Zeilen | ",
          sum(dat$presence == 1), " Präsenzen | ",
          sum(dat$presence == 0), " Background")
  
  # Skalierungsparameter zurückgeben
  attr(dat, "scale_center") <- scale_center
  attr(dat, "scale_scale")  <- scale_scale
  
  return(dat)
}

# =============================================================================
# 3 - MaxNet trainieren (finales Modell, keine CV/Evaluierung)
# =============================================================================

run_maxnet <- function(sampling, env_raster) {
  
  message("  -> Daten vorbereiten...")
  dat <- prepare_maxnet_data(sampling, env_raster)
  
  if (nrow(dat) == 0) stop("dat ist leer!")
  
  # Skalierungsparameter aus dat extrahieren
  scale_center <- attr(dat, "scale_center")
  scale_scale  <- attr(dat, "scale_scale")
  
  env_cols <- names(env_raster)
  
  # Raster skalieren für finale Vorhersage
  env_raster_scaled <- (env_raster - scale_center) / scale_scale
  names(env_raster_scaled) <- env_cols
  
  # Finales Modell auf ALLEN Daten
  message("  -> Finales Modell fitten...")
  all_env <- dat[, env_cols, drop = FALSE]
  
  final_model <- maxnet(
    p = dat$presence,
    data = all_env,
    f = maxnet.formula(dat$presence, all_env, classes = "default")
  )
  
  # Vorhersage auf skaliertem Raster
  message("  -> Vorhersage auf Raster...")
  final_pred <- terra::predict(
    env_raster_scaled,
    final_model,
    fun = function(model, data) {
      predict(model, data, type = "cloglog")
    },
    na.rm = TRUE
  )
  
  message("  -> Fertig!")
  
  list(
    model        = final_model,
    prediction   = final_pred,
    data         = dat,
    scale_center = scale_center,
    scale_scale  = scale_scale
  )
}

# =============================================================================
# 4 - Sampling-Daten laden + MaxNet ausführen
# =============================================================================

sampling_narrow   <- readRDS("Data/species/sampling_narrow.RDS")
sampling_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
sampling_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
sampling_broad    <- readRDS("Data/species/sampling_broad.RDS")

message("===== Narrow =====")
mx_narrow   <- run_maxnet(sampling_narrow,   env_raster_masked)

message("===== Low-Mid =====")
mx_low_mid  <- run_maxnet(sampling_low_mid,  env_raster_masked)

message("===== High-Mid =====")
mx_high_mid <- run_maxnet(sampling_high_mid, env_raster_masked)

message("===== Broad =====")
mx_broad    <- run_maxnet(sampling_broad,    env_raster_masked)

# =============================================================================
# 5 - Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
plot(mx_narrow$prediction,   main = "MaxNet Narrow")
plot(mx_low_mid$prediction,  main = "MaxNet Low-Mid")
plot(mx_high_mid$prediction, main = "MaxNet High-Mid")
plot(mx_broad$prediction,    main = "MaxNet Broad")
par(mfrow = c(1, 1))

# =============================================================================
# 6 - Speichern
# =============================================================================

dir.create("Data/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(mx_narrow,   "Data/models/maxnet_narrow.RDS")
saveRDS(mx_low_mid,  "Data/models/maxnet_low_mid.RDS")
saveRDS(mx_high_mid, "Data/models/maxnet_high_mid.RDS")
saveRDS(mx_broad,    "Data/models/maxnet_broad.RDS")'
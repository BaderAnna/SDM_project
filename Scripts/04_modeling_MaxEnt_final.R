# =============================================================================
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
saveRDS(mx_broad,    "Data/models/maxnet_broad.RDS")
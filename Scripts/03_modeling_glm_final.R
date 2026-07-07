# =============================================================================
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
# 4 - GLM trainieren (finales Modell, keine CV/Evaluierung)
# =============================================================================

run_glm <- function(sampling, env_raster, glm_formula) {
  
  dat <- prepare_glm_data(sampling, env_raster)
  
  if (nrow(dat) == 0) stop("dat ist nach prepare_glm_data() leer!")
  
  n_pres <- sum(dat$presence == 1)
  n_abs  <- sum(dat$presence == 0)
  dat$.weights <- ifelse(dat$presence == 1, 1, n_pres / n_abs)
  
  final_model <- suppressWarnings(
    glm(glm_formula, data = dat, family = binomial, weights = .weights)
  )
  
  final_pred <- terra::predict(env_raster, final_model, type = "response")
  
  list(model = final_model, prediction = final_pred, data = dat)
}

# =============================================================================
# 5 - Daten laden + GLM ausführen
# =============================================================================

sampling_narrow   <- readRDS("Data/species/sampling_narrow.RDS")
sampling_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
sampling_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
sampling_broad    <- readRDS("Data/species/sampling_broad.RDS")

message("===== Narrow =====")
glm_narrow <- run_glm(sampling_narrow, env_raster_masked, glm_formula)

message("===== Low-Mid =====")
glm_low_mid <- run_glm(sampling_low_mid, env_raster_masked, glm_formula)

message("===== High-Mid =====")
glm_high_mid <- run_glm(sampling_high_mid, env_raster_masked, glm_formula)

message("===== Broad =====")
glm_broad <- run_glm(sampling_broad, env_raster_masked, glm_formula)

# =============================================================================
# 6 - Speichern
# =============================================================================

dir.create("Data/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(glm_narrow,   "Data/models/glm_narrow.RDS")
saveRDS(glm_low_mid,  "Data/models/glm_low_mid.RDS")
saveRDS(glm_high_mid, "Data/models/glm_high_mid.RDS")
saveRDS(glm_broad,    "Data/models/glm_broad.RDS")


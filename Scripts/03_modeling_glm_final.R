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


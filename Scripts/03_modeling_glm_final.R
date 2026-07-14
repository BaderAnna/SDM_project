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
# 4 - knndm Fold-Erstellung
# =============================================================================

run_knndm_folds <- function(dat_sf_geom, predictor_raster,
                            k = 5, seed, n_predpoints = 2000) {
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
# 5 - GLM mit knndm: Nur finale Modelle 
# =============================================================================

run_glm_knndm_final <- function(sampling, env_raster, glm_formula) {
  
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
# 4 - knndm Fold-Erstellung
# =============================================================================

run_knndm_folds <- function(dat_sf_geom, predictor_raster,
                            k = 5, seed, n_predpoints = 2000) {
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
# 5 - GLM mit knndm Spatial-CV
# =============================================================================

run_glm_knndm <- function(sampling, env_raster, glm_formula,
                          n_runs = 5, k = 5) {
  
  dat <- prepare_glm_data(sampling, env_raster)
  
  if (nrow(dat) == 0)
    stop("dat ist nach prepare_glm_data() leer!")
  
  dat_sf      <- sf::st_as_sf(dat, coords = c("x", "y"), crs = 3035)
  dat_sf_geom <- sf::st_as_sf(
    data.frame(geometry = sf::st_geometry(dat_sf)),
    crs = sf::st_crs(dat_sf)
  )
  
  auc_scores <- numeric(n_runs)
  tss_scores <- numeric(n_runs)
  
  for (run in seq_len(n_runs)) {
    
    knn <- tryCatch(
      run_knndm_folds(dat_sf_geom, env_raster, k = k, seed = run),
      error = function(e) {
        message("knndm Fehler in Run ", run, ": ", e$message)
        NULL
      }
    )
    if (is.null(knn)) next
    
    fold_ids    <- knn$clusters
    fold_preds  <- rep(NA_real_, nrow(dat))
    fold_truths <- rep(NA_real_, nrow(dat))
    
    for (fold in seq_len(k)) {
      test_idx  <- which(fold_ids == fold)
      train_idx <- which(fold_ids != fold)
      if (length(test_idx) == 0 || length(train_idx) == 0) next
      
      train <- dat[train_idx, ]
      test  <- dat[test_idx,  ]
      if (length(unique(train$presence)) < 2) next
      
      n_pres  <- sum(train$presence == 1)
      n_abs   <- sum(train$presence == 0)
      
      # gewichte direkt in train einfügen
      train$.weights <- ifelse(train$presence == 1, 1, n_pres / n_abs)
      
      model <- suppressWarnings(
        glm(glm_formula, data = train,
            family = binomial, weights = .weights)  
      )
      
      fold_preds[test_idx]  <- predict(model, newdata = test,
                                       type = "response")
      fold_truths[test_idx] <- test$presence
    }
    
    valid        <- !is.na(fold_preds) & !is.na(fold_truths)
    preds_valid  <- fold_preds[valid]
    truths_valid <- fold_truths[valid]
    
    if (length(preds_valid) < 10 || length(unique(truths_valid)) < 2) next
    
    auc <- as.numeric(
      pROC::auc(pROC::roc(truths_valid, preds_valid, quiet = TRUE))
    )
    
    pred_bin <- ifelse(preds_valid >= 0.5, 1, 0)
    tp   <- sum(pred_bin == 1 & truths_valid == 1)
    fp   <- sum(pred_bin == 1 & truths_valid == 0)
    tn   <- sum(pred_bin == 0 & truths_valid == 0)
    fn   <- sum(pred_bin == 0 & truths_valid == 1)
    sens <- tp / (tp + fn)
    spec <- tn / (tn + fp)
    tss  <- sens + spec - 1
    
    auc_scores[run] <- auc
    tss_scores[run] <- tss
    
    message(sprintf("Run %d/%d | AUC = %.3f | TSS = %.3f",
                    run, n_runs, auc, tss))
  }
  
  # Finales Modell auf ALLEN Daten
  n_pres <- sum(dat$presence == 1)
  n_abs  <- sum(dat$presence == 0)
  
  # gewichte direkt in dat einfügen
  dat$.weights <- ifelse(dat$presence == 1, 1, n_pres / n_abs)
  
  final_model <- suppressWarnings(
    glm(glm_formula, data = dat,
        family = binomial, weights = .weights)  
  )
  
  final_pred <- terra::predict(env_raster, final_model, type = "response")
  
  valid_auc <- auc_scores[auc_scores != 0]
  valid_tss <- tss_scores[tss_scores != 0]
  
  list(
    model      = final_model,
    prediction = final_pred,
    data       = dat,
    auc_mean   = mean(valid_auc, na.rm = TRUE),
    auc_sd     = sd(valid_auc,   na.rm = TRUE),
    tss_mean   = mean(valid_tss, na.rm = TRUE),
    tss_sd     = sd(valid_tss,   na.rm = TRUE),
    auc_all    = auc_scores,
    tss_all    = tss_scores
  )
}

# =============================================================================
# 4 - GLM trainieren (finales Modell, keine CV/Evaluierung)
# =============================================================================

run_glm_knndm <- function(sampling, env_raster, glm_formula) {
  
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

message("===== Narrow (kNNDM) =====")
glm_narrow_knndm <- run_glm_knndm(
  sampling = sampling_narrow,
  env_raster = env_raster_masked,
  glm_formula = glm_formula,
  n_runs = 5,  # Du kannst auch n_runs = 1 setzen, wenn du nur eine Fold-Konfiguration willst
  k = 5
)

message("===== Low-Mid (kNNDM) =====")
glm_low_mid_knndm <- run_glm_knndm(
  sampling = sampling_low_mid,
  env_raster = env_raster_masked,
  glm_formula = glm_formula,
  n_runs = 5,
  k = 5
)

message("===== High-Mid (kNNDM) =====")
glm_high_mid_knndm <- run_glm_knndm(
  sampling = sampling_high_mid,
  env_raster = env_raster_masked,
  glm_formula = glm_formula,
  n_runs = 5,
  k = 5
)

message("===== Broad (kNNDM) =====")
glm_broad_knndm <- run_glm_knndm(
  sampling = sampling_broad,
  env_raster = env_raster_masked,
  glm_formula = glm_formula,
  n_runs = 5,
  k = 5
)

# =============================================================================
# 6 - Speichern der finalen Modelle (mit kNNDM)
# =============================================================================

dir.create("Data/models/knndm", recursive = TRUE, showWarnings = FALSE)

saveRDS(glm_narrow_knndm,   "Data/models/knndm/glm_narrow_knndm.RDS")
saveRDS(glm_low_mid_knndm,  "Data/models/knndm/glm_low_mid_knndm.RDS")
saveRDS(glm_high_mid_knndm, "Data/models/knndm/glm_high_mid_knndm.RDS")
saveRDS(glm_broad_knndm,    "Data/models/knndm/glm_broad_knndm.RDS")'

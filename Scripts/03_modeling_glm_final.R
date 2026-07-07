# =============================================================================
# GLM MODELLIERUNG
# Option B: alle 19 WorldClim Variablen -> VIF-Reduktion
# Spatial CV: CAST::knndm, 5 Wiederholungen
# =============================================================================

library(terra)
library(sf)
library(usdm)
library(CAST)
library(rnaturalearth)
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

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
  
  # --- Präsenzpunkte ---
  sp   <- sampling$po$sample.points
  pres <- sp[sp$Observed == TRUE, c("x", "y")]
  pres$presence <- 1
  
  # --- Hintergrundpunkte ---
  abs <- as.data.frame(sampling$bg_glm)
  abs$presence <- 0
  
  # --- Kombinieren ---
  dat <- rbind(pres, abs)
  
  if (nrow(dat) == 0)
    stop("dat ist leer!")
  if (any(is.na(dat$x)) || any(is.na(dat$y)))
    stop("x oder y enthalten NA-Werte!")
  
  # Koordinaten sind bereits EPSG:3035 -> direkt als Matrix verwenden
  coords_3035 <- as.matrix(dat[, c("x", "y")])
  
  # Kontrolle
  raster_ext <- terra::ext(env_raster)
  message("Raster Extent  : xmin=", round(raster_ext$xmin),
          " xmax=", round(raster_ext$xmax),
          " ymin=", round(raster_ext$ymin),
          " ymax=", round(raster_ext$ymax))
  message("Koordinaten x  : ", round(min(coords_3035[,1])),
          " bis ", round(max(coords_3035[,1])))
  message("Koordinaten y  : ", round(min(coords_3035[,2])),
          " bis ", round(max(coords_3035[,2])))
  
  # --- Umweltwerte extrahieren ---
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords_3035))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  message("NA-Anteil env_vals: ",
          round(mean(is.na(env_vals)) * 100, 1), "%")
  
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
# 6 - Daten laden + GLM ausführen
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars     <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars  # Namen sicherstellen
glm_formula       <- build_glm_formula(selected_vars)

sampling_narrow       <- readRDS("Data/species/sampling_narrow.RDS")
sampling_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
sampling_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
sampling_broad        <- readRDS("Data/species/sampling_broad.RDS")

message("===== Narrow =====")
glm_narrow <- run_glm_knndm(sampling_narrow, env_raster_masked,
                            glm_formula, n_runs = 5)

message("===== Low-Mid =====")
glm_low_mid <- run_glm_knndm(sampling_low_mid, env_raster_masked,
                                  glm_formula, n_runs = 5)

message("===== High-Mid =====")
glm_high_mid <- run_glm_knndm(sampling_high_mid, env_raster_masked,
                             glm_formula, n_runs = 5)

message("===== Broad =====")
glm_broad <- run_glm_knndm(sampling_broad, env_raster_masked,
                           glm_formula, n_runs = 5)

# =============================================================================
# 7 - Ergebnisse zusammenfassen
# =============================================================================

results_glm <- data.frame(
  species  = c("narrow", "low_mid", "high_mid", "broad"),
  sigma    = c(0.2, 0.4, 0.6, 0.8),
  auc_mean = c(glm_narrow$auc_mean,
               glm_low_mid$auc_mean,
               glm_high_mid$auc_mean,
               glm_broad$auc_mean),
  auc_sd   = c(glm_narrow$auc_sd,
               glm_low_mid$auc_sd,
               glm_high_mid$auc_sd,
               glm_broad$auc_sd),
  tss_mean = c(glm_narrow$tss_mean,
               glm_low_mid$tss_mean,
               glm_high_mid$tss_mean,
               glm_broad$tss_mean),
  tss_sd   = c(glm_narrow$tss_sd,
               glm_low_mid$tss_sd,
               glm_high_mid$tss_sd,
               glm_broad$tss_sd)
)

print(results_glm)

# =============================================================================
# 8 - Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
plot(glm_narrow$prediction,
     main = sprintf("GLM Narrow\nAUC=%.3f | TSS=%.3f",
                    glm_narrow$auc_mean, glm_narrow$tss_mean))
plot(glm_low_mid$prediction,
     main = sprintf("GLM Low-Mid\nAUC=%.3f | TSS=%.3f",
                    glm_low_mid$auc_mean, glm_low_mid$tss_mean))
plot(glm_high_mid$prediction,
     main = sprintf("GLM High-Mid\nAUC=%.3f | TSS=%.3f",
                    glm_high_mid$auc_mean, glm_high_mid$tss_mean))
plot(glm_broad$prediction,
     main = sprintf("GLM Broad\nAUC=%.3f | TSS=%.3f",
                    glm_broad$auc_mean, glm_broad$tss_mean))
par(mfrow = c(1, 1))

# =============================================================================
# 9 - Speichern
# =============================================================================

dir.create("Data/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(glm_narrow,   "Data/models/glm_narrow.RDS")
saveRDS(glm_low_mid,  "Data/models/glm_low_mid.RDS")
saveRDS(glm_high_mid, "Data/models/glm_high_mid.RDS")
saveRDS(glm_broad,    "Data/models/glm_broad.RDS")
saveRDS(results_glm,  "Data/models/results_glm.RDS")


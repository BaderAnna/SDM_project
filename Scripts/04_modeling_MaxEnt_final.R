# =============================================================================
# MAXNET MODELLIERUNG
# WorldClim Variablen nach VIF-Reduktion
# Spatial CV: CAST::knndm, 1 Durchlauf
# Für narrow, low_mid, high_mid, broad
# =============================================================================

library(terra)
library(sf)
library(maxnet)
library(CAST)
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

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
  
  # ✅ Variablen standardisieren (z-Transformation)
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
# 3 - knndm Fold-Erstellung
# =============================================================================

run_knndm_folds <- function(dat_sf_geom, predictor_raster,
                            k = 5, seed = 42, n_predpoints = 2000) {
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
# 4 - MaxNet mit knndm Spatial-CV (1 Durchlauf)
# =============================================================================

run_maxnet_knndm <- function(sampling, env_raster, k = 5, seed = 42) {
  
  message("  -> Daten vorbereiten...")
  dat <- prepare_maxnet_data(sampling, env_raster)
  
  if (nrow(dat) == 0) stop("dat ist leer!")
  
  # Skalierungsparameter aus dat extrahieren
  scale_center <- attr(dat, "scale_center")
  scale_scale  <- attr(dat, "scale_scale")
  
  env_cols <- names(env_raster)
  
  # ✅ Raster skalieren für finale Vorhersage
  env_raster_scaled <- (env_raster - scale_center) / scale_scale
  names(env_raster_scaled) <- env_cols
  
  # sf-Objekt für knndm
  dat_sf      <- sf::st_as_sf(dat, coords = c("x", "y"), crs = 3035)
  dat_sf_geom <- sf::st_as_sf(
    data.frame(geometry = sf::st_geometry(dat_sf)),
    crs = sf::st_crs(dat_sf)
  )
  
  message("  -> knndm Folds erstellen (k = ", k, ")...")
  knn <- tryCatch(
    run_knndm_folds(dat_sf_geom, env_raster, k = k, seed = seed),
    error = function(e) {
      message("knndm Fehler: ", e$message)
      NULL
    }
  )
  
  auc <- NA_real_
  tss <- NA_real_
  
  if (!is.null(knn)) {
    
    fold_ids    <- knn$clusters
    fold_preds  <- rep(NA_real_, nrow(dat))
    fold_truths <- rep(NA_real_, nrow(dat))
    
    for (fold in seq_len(k)) {
      
      test_idx  <- which(fold_ids == fold)
      train_idx <- which(fold_ids != fold)
      if (length(test_idx) == 0 || length(train_idx) == 0) next
      
      train <- dat[train_idx, ]
      test  <- dat[test_idx,  ]
      
      n_pres_train <- sum(train$presence == 1)
      n_pres_test  <- sum(test$presence  == 1)
      
      message(sprintf("  Fold %d/%d: %d Pres-Train | %d Pres-Test | %d BG-Train",
                      fold, k, n_pres_train, n_pres_test,
                      sum(train$presence == 0)))
      
      if (n_pres_train < 10 || length(unique(train$presence)) < 2) {
        message("  -> Zu wenige Präsenzen -> überspringe Fold ", fold)
        next
      }
      
      train_env <- train[, env_cols, drop = FALSE]
      test_env  <- test[,  env_cols, drop = FALSE]
      
      model <- tryCatch(
        maxnet(
          p = train$presence,
          data = train_env,
          f = maxnet.formula(train$presence, train_env,
                             classes = "default")
        ),
        error = function(e) {
          message("  MaxNet Fehler in Fold ", fold, ": ", e$message)
          NULL
        }
      )
      if (is.null(model)) next
      
      fold_preds[test_idx]  <- predict(model, test_env, type = "cloglog")
      fold_truths[test_idx] <- test$presence
    }
    
    # Metriken
    valid        <- !is.na(fold_preds) & !is.na(fold_truths)
    preds_valid  <- fold_preds[valid]
    truths_valid <- fold_truths[valid]
    
    if (length(preds_valid) >= 10 && length(unique(truths_valid)) == 2) {
      
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
      
      message(sprintf("  -> AUC = %.3f | TSS = %.3f", auc, tss))
    } else {
      message("  -> Nicht genug gültige Vorhersagen!")
    }
  }
  
  # Finales Modell auf ALLEN Daten
  message("  -> Finales Modell fitten...")
  all_env <- dat[, env_cols, drop = FALSE]
  
  final_model <- maxnet(
    p = dat$presence,
    data = all_env,
    f = maxnet.formula(dat$presence, all_env, classes = "default")
  )
  
  # ✅ Vorhersage auf skaliertem Raster
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
    scale_scale  = scale_scale,
    auc          = auc,
    tss          = tss
  )
}

# =============================================================================
# 5 - Sampling-Daten laden + MaxNet ausführen
# =============================================================================

sampling_narrow   <- readRDS("Data/species/sampling_narrow.RDS")
sampling_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
sampling_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
sampling_broad    <- readRDS("Data/species/sampling_broad.RDS")

message("===== Narrow =====")
mx_narrow   <- run_maxnet_knndm(sampling_narrow,   env_raster_masked, k = 5)

message("===== Low-Mid =====")
mx_low_mid  <- run_maxnet_knndm(sampling_low_mid,  env_raster_masked, k = 5)

message("===== High-Mid =====")
mx_high_mid <- run_maxnet_knndm(sampling_high_mid, env_raster_masked, k = 5)

message("===== Broad =====")
mx_broad    <- run_maxnet_knndm(sampling_broad,    env_raster_masked, k = 5)

# =============================================================================
# 6 - Ergebnisse zusammenfassen
# =============================================================================

results_maxnet <- data.frame(
  species = c("narrow", "low_mid", "high_mid", "broad"),
  sigma   = c(0.2, 0.4, 0.6, 0.8),
  auc     = c(mx_narrow$auc,
              mx_low_mid$auc,
              mx_high_mid$auc,
              mx_broad$auc),
  tss     = c(mx_narrow$tss,
              mx_low_mid$tss,
              mx_high_mid$tss,
              mx_broad$tss)
)

print(results_maxnet)

# =============================================================================
# 7 - Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
plot(mx_narrow$prediction,
     main = sprintf("MaxNet Narrow\nAUC=%.3f | TSS=%.3f",
                    mx_narrow$auc, mx_narrow$tss))
plot(mx_low_mid$prediction,
     main = sprintf("MaxNet Low-Mid\nAUC=%.3f | TSS=%.3f",
                    mx_low_mid$auc, mx_low_mid$tss))
plot(mx_high_mid$prediction,
     main = sprintf("MaxNet High-Mid\nAUC=%.3f | TSS=%.3f",
                    mx_high_mid$auc, mx_high_mid$tss))
plot(mx_broad$prediction,
     main = sprintf("MaxNet Broad\nAUC=%.3f | TSS=%.3f",
                    mx_broad$auc, mx_broad$tss))
par(mfrow = c(1, 1))

# =============================================================================
# 8 - Speichern
# =============================================================================

dir.create("Data/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(mx_narrow,      "Data/models/maxnet_narrow.RDS")
saveRDS(mx_low_mid,     "Data/models/maxnet_low_mid.RDS")
saveRDS(mx_high_mid,    "Data/models/maxnet_high_mid.RDS")
saveRDS(mx_broad,       "Data/models/maxnet_broad.RDS")
saveRDS(results_maxnet, "Data/models/results_maxnet.RDS")

message("✅ Alle MaxNet-Modelle gespeichert!")
print(results_maxnet)
# =============================================================================
# BRT MODELLIERUNG
# WorldClim Variablen nach VIF-Reduktion
# Spatial CV: CAST::knndm, 1 Durchlauf
# 10 Runs pro Art (unterschiedliche pseudo-absences), dann mitteln
# Für narrow, low_mid, high_mid, broad
# =============================================================================

library(terra)
library(sf)
library(gbm)
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
# 2 - GLM-Formel für BRT bauen
# =============================================================================

build_brt_formula <- function(vars) {
  as.formula(paste("presence ~", paste(vars, collapse = " + ")))
}

brt_formula <- build_brt_formula(selected_vars)
print(brt_formula)

# =============================================================================
# 3 - Hilfsfunktionen
# =============================================================================

# --- Daten vorbereiten für einen BRT-Run ---
# Koordinaten bereits EPSG:3035 -> keine Transformation!
prepare_brt_data <- function(run_data, env_raster) {
  
  # Koordinaten bereits EPSG:3035 -> direkt als Matrix
  coords <- as.matrix(run_data[, c("x", "y")])
  
  # Umweltwerte extrahieren
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  dat <- cbind(run_data, env_vals)
  dat <- na.omit(dat)
  
  return(dat)
}

# --- knndm Fold-Erstellung ---
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
# 4 - BRT mit knndm Spatial-CV
#     10 Runs mit unterschiedlichen pseudo-absences -> mitteln
# =============================================================================

run_brt_knndm <- function(sampling, env_raster, brt_formula,
                          n_runs = 10, k = 5,
                          n_trees = 2000, interaction_depth = 3,
                          shrinkage = 0.01, bag_fraction = 0.75) {
  
  env_cols   <- names(env_raster)
  
  auc_scores <- numeric(n_runs)
  tss_scores <- numeric(n_runs)
  pred_list  <- list()
  model_list <- list()
  
  for (run in seq_len(n_runs)) {
    
    message(sprintf("  -> Run %d/%d", run, n_runs))
    
    # --- Daten für diesen Run vorbereiten ---
    # runs enthält 10 verschiedene pseudo-absence Datensätze
    run_data <- prepare_brt_data(sampling$runs[[run]], env_raster)
    
    if (nrow(run_data) == 0) {
      message("  -> Run ", run, ": dat leer -> überspringe")
      next
    }
    
    message(sprintf("     %d Zeilen | %d Präsenzen | %d Absences",
                    nrow(run_data),
                    sum(run_data$presence == 1),
                    sum(run_data$presence == 0)))
    
    # sf-Objekt für knndm - NUR Geometrie
    dat_sf      <- sf::st_as_sf(run_data, coords = c("x", "y"), crs = 3035)
    dat_sf_geom <- sf::st_as_sf(
      data.frame(geometry = sf::st_geometry(dat_sf)),
      crs = sf::st_crs(dat_sf)
    )
    
    # --- knndm Folds ---
    knn <- tryCatch(
      run_knndm_folds(dat_sf_geom, env_raster, k = k, seed = run),
      error = function(e) {
        message("  knndm Fehler in Run ", run, ": ", e$message)
        NULL
      }
    )
    
    if (!is.null(knn)) {
      
      fold_ids    <- knn$clusters
      fold_preds  <- rep(NA_real_, nrow(run_data))
      fold_truths <- rep(NA_real_, nrow(run_data))
      
      for (fold in seq_len(k)) {
        
        test_idx  <- which(fold_ids == fold)
        train_idx <- which(fold_ids != fold)
        if (length(test_idx) == 0 || length(train_idx) == 0) next
        
        train <- run_data[train_idx, ]
        test  <- run_data[test_idx,  ]
        if (length(unique(train$presence)) < 2) next
        
        # BRT fitten
        brt_cv <- tryCatch(
          gbm::gbm(
            brt_formula,
            data              = train,
            distribution      = "bernoulli",
            n.trees           = n_trees,
            interaction.depth = interaction_depth,
            shrinkage         = shrinkage,
            bag.fraction      = bag_fraction,
            cv.folds          = 5,
            verbose           = FALSE
          ),
          error = function(e) {
            message("  BRT Fehler in Fold ", fold, ": ", e$message)
            NULL
          }
        )
        if (is.null(brt_cv)) next
        
        best_trees <- gbm::gbm.perf(brt_cv, method = "cv",
                                    plot.it = FALSE)
        
        fold_preds[test_idx] <- predict(
          brt_cv, newdata = test,
          n.trees = best_trees, type = "response"
        )
        fold_truths[test_idx] <- test$presence
      }
      
      # Metriken berechnen
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
        
        auc_scores[run] <- auc
        tss_scores[run] <- tss
        
        message(sprintf("     AUC = %.3f | TSS = %.3f", auc, tss))
      }
    }
    
    # --- Finales Modell für diesen Run auf ALLEN Daten ---
    brt_final <- tryCatch(
      gbm::gbm(
        brt_formula,
        data              = run_data,
        distribution      = "bernoulli",
        n.trees           = n_trees,
        interaction.depth = interaction_depth,
        shrinkage         = shrinkage,
        bag.fraction      = bag_fraction,
        cv.folds          = 5,
        verbose           = FALSE
      ),
      error = function(e) {
        message("  BRT Final-Fehler in Run ", run, ": ", e$message)
        NULL
      }
    )
    if (is.null(brt_final)) next
    
    best_trees_final <- gbm::gbm.perf(brt_final, method = "cv",
                                      plot.it = FALSE)
    
    # Vorhersage auf Raster
    pred_run <- terra::predict(
      env_raster,
      brt_final,
      n.trees = best_trees_final,
      type    = "response",
      na.rm   = TRUE
    )
    
    model_list[[run]] <- brt_final
    pred_list[[run]]  <- pred_run
  }
  
  # --- Über alle Runs mitteln ---
  message("  -> Mittele über ", length(pred_list), " Runs...")
  final_pred <- terra::app(terra::rast(pred_list), mean)
  
  # Nur gültige Scores
  valid_auc <- auc_scores[auc_scores != 0]
  valid_tss <- tss_scores[tss_scores != 0]
  
  message(sprintf("  -> AUC mean = %.3f (sd = %.3f) | TSS mean = %.3f (sd = %.3f)",
                  mean(valid_auc, na.rm = TRUE),
                  sd(valid_auc,   na.rm = TRUE),
                  mean(valid_tss, na.rm = TRUE),
                  sd(valid_tss,   na.rm = TRUE)))
  
  list(
    models     = model_list,
    prediction = final_pred,
    auc_mean   = mean(valid_auc, na.rm = TRUE),
    auc_sd     = sd(valid_auc,   na.rm = TRUE),
    tss_mean   = mean(valid_tss, na.rm = TRUE),
    tss_sd     = sd(valid_tss,   na.rm = TRUE),
    auc_all    = auc_scores,
    tss_all    = tss_scores
  )
}

# =============================================================================
# 5 - Sampling-Daten laden + BRT ausführen
# =============================================================================

sampling_narrow   <- readRDS("Data/species/sampling_narrow.RDS")
sampling_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
sampling_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
sampling_broad    <- readRDS("Data/species/sampling_broad.RDS")

message("===== Narrow =====")
brt_narrow   <- run_brt_knndm(sampling_narrow,   env_raster_masked,
                              brt_formula, n_runs = 10)

message("===== Low-Mid =====")
brt_low_mid  <- run_brt_knndm(sampling_low_mid,  env_raster_masked,
                              brt_formula, n_runs = 10)

message("===== High-Mid =====")
brt_high_mid <- run_brt_knndm(sampling_high_mid, env_raster_masked,
                              brt_formula, n_runs = 10)

message("===== Broad =====")
brt_broad    <- run_brt_knndm(sampling_broad,    env_raster_masked,
                              brt_formula, n_runs = 10)

# =============================================================================
# 6 - Ergebnisse zusammenfassen
# =============================================================================

results_brt <- data.frame(
  species  = c("narrow", "low_mid", "high_mid", "broad"),
  sigma    = c(0.2, 0.4, 0.6, 0.8),
  auc_mean = c(brt_narrow$auc_mean,
               brt_low_mid$auc_mean,
               brt_high_mid$auc_mean,
               brt_broad$auc_mean),
  auc_sd   = c(brt_narrow$auc_sd,
               brt_low_mid$auc_sd,
               brt_high_mid$auc_sd,
               brt_broad$auc_sd),
  tss_mean = c(brt_narrow$tss_mean,
               brt_low_mid$tss_mean,
               brt_high_mid$tss_mean,
               brt_broad$tss_mean),
  tss_sd   = c(brt_narrow$tss_sd,
               brt_low_mid$tss_sd,
               brt_high_mid$tss_sd,
               brt_broad$tss_sd)
)

print(results_brt)

# =============================================================================
# 7 - Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
plot(brt_narrow$prediction,
     main = sprintf("BRT Narrow\nAUC=%.3f | TSS=%.3f",
                    brt_narrow$auc_mean, brt_narrow$tss_mean))
plot(brt_low_mid$prediction,
     main = sprintf("BRT Low-Mid\nAUC=%.3f | TSS=%.3f",
                    brt_low_mid$auc_mean, brt_low_mid$tss_mean))
plot(brt_high_mid$prediction,
     main = sprintf("BRT High-Mid\nAUC=%.3f | TSS=%.3f",
                    brt_high_mid$auc_mean, brt_high_mid$tss_mean))
plot(brt_broad$prediction,
     main = sprintf("BRT Broad\nAUC=%.3f | TSS=%.3f",
                    brt_broad$auc_mean, brt_broad$tss_mean))
par(mfrow = c(1, 1))

# =============================================================================
# 8 - Speichern
# =============================================================================

dir.create("Data/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(brt_narrow,   "Data/models/brt_narrow.RDS")
saveRDS(brt_low_mid,  "Data/models/brt_low_mid.RDS")
saveRDS(brt_high_mid, "Data/models/brt_high_mid.RDS")
saveRDS(brt_broad,    "Data/models/brt_broad.RDS")
saveRDS(results_brt,  "Data/models/results_brt.RDS")

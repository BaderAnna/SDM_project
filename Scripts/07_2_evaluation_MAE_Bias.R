# =============================================================================
# Berechnung: MAE + Bias (vs. wahre Suitability der virtuellen Spezies)
# =============================================================================
#
# obs  = wahre, kontinuierliche Suitability an Testkoordinaten
#        (aus virtualspecies::generateSpFromFun(), NICHT die binaere P/A)
# pred = kontinuierliche Modellvorhersage an denselben Koordinaten
#
# MAE  = mean(|obs - pred|)     -> je kleiner, desto besser (0 = perfekt)
# Bias = mean(pred - obs)       -> positiv = Ueberschaetzung, negativ = Unterschaetzung
#
# MAE und Bias nutzen exakt dieselbe Vorhersage + dieselbe wahre Suitability,
# daher werden sie hier in einem Durchgang berechnet (keine doppelte
# Extraktion/Vorhersage noetig).
#
# Output: results/mae_bias_true_suitability.csv
#         -> wird von plot_mae_bias_true_suitability.R eingelesen
# =============================================================================

library(Metrics)
library(terra)
library(dplyr)

# -----------------------------------------------------------------------------
# 1 - Daten laden (an eure Struktur anpassen!)
# -----------------------------------------------------------------------------
env <- terra::rast("Data/raster/env_raster_masked.tif")

models_glm <- list(
  narrow   = readRDS("Data/models/knndm/glm_narrow_train.RDS"),
  low_mid  = readRDS("Data/models/knndm/glm_low_mid_train.RDS"),
  high_mid = readRDS("Data/models/knndm/glm_high_mid_train.RDS"),
  broad    = readRDS("Data/models/knndm/glm_broad_train.RDS")
)

models_maxent <- list(
  narrow   = readRDS("Data/models/knndm/maxnet_narrow_train.RDS"),
  low_mid  = readRDS("Data/models/knndm/maxnet_low_mid_train.RDS"),
  high_mid = readRDS("Data/models/knndm/maxnet_high_mid_train.RDS"),
  broad    = readRDS("Data/models/knndm/maxnet_broad_train.RDS")
)

models_brt <- list(
  narrow   = readRDS("Data/models/knndm/brt_narrow_train.RDS"),
  low_mid  = readRDS("Data/models/knndm/brt_low_mid_train.RDS"),
  high_mid = readRDS("Data/models/knndm/brt_high_mid_train.RDS"),
  broad    = readRDS("Data/models/knndm/brt_broad_train.RDS")
)

models_rf <- list(
  narrow   = readRDS("Data/models/knndm/rf_narrow_train.RDS"),
  low_mid  = readRDS("Data/models/knndm/rf_low_mid_train.RDS"),
  high_mid = readRDS("Data/models/knndm/rf_high_mid_train.RDS"),
  broad    = readRDS("Data/models/knndm/rf_broad_train.RDS")
)

models_all <- list(GLM = models_glm, Maxent = models_maxent, BRT = models_brt, RF = models_rf)

test_data <- list(
  narrow   = readRDS("Data/species/split_knndm/test_narrow.RDS"),
  low_mid  = readRDS("Data/species/split_knndm/test_low_mid.RDS"),
  high_mid = readRDS("Data/species/split_knndm/test_high_mid.RDS"),
  broad    = readRDS("Data/species/split_knndm/test_broad.RDS")
)

true_suitability_rasters <- list(
  narrow   = terra::unwrap(readRDS("Data/species/species_narrow.RDS")$virtual$suitab.raster),
  low_mid  = terra::unwrap(readRDS("Data/species/species_low_mid.RDS")$virtual$suitab.raster),
  high_mid = terra::unwrap(readRDS("Data/species/species_high_mid.RDS")$virtual$suitab.raster),
  broad    = terra::unwrap(readRDS("Data/species/species_broad.RDS")$virtual$suitab.raster)
)

niche_names <- c("narrow", "low_mid", "high_mid", "broad")


# -----------------------------------------------------------------------------
# 2 - Vorhersagefunktion je Modelltyp
# -----------------------------------------------------------------------------

predict_model <- function(model_obj, newdata, model_type) {
  switch(model_type,
         
         "GLM" = predict(model_obj$model, 
                         newdata = newdata, 
                         type    = "response"),
         
         "Maxent" = {
           scale_center <- model_obj$scale_center
           scale_scale  <- model_obj$scale_scale
           env_vars     <- names(scale_center)
           
           newdata_scaled <- newdata
           newdata_scaled[, env_vars] <- scale(newdata[, env_vars],
                                               center = scale_center,
                                               scale  = scale_scale)
           
           pred_raw <- predict(model_obj$model,
                               newdata = newdata_scaled,
                               type    = "cloglog")
           
           # Sicherheitscheck: immer Vektor der Länge nrow(newdata) zurückgeben
           pred_vec <- as.vector(pred_raw)
           
           if (length(pred_vec) != nrow(newdata)) {
             warning(sprintf(
               "[Maxent] Unerwartete Länge: %d statt %d → nehme erste Spalte",
               length(pred_vec), nrow(newdata)
             ))
             pred_vec <- as.vector(pred_raw[, 1])
           }
           pred_vec
         },
         
         "BRT" = {
           gbm_models <- model_obj$models
           n_trees    <- gbm_models[[1]]$n.trees
           preds <- sapply(gbm_models, function(gbm_fit) {
             gbm::predict.gbm(gbm_fit, 
                              newdata = newdata, 
                              n.trees = n_trees, 
                              type    = "response")
           })
           rowMeans(preds)
         },
         
         "RF" = {
           ranger_models <- model_obj$models
           preds <- sapply(ranger_models, function(rf_fit) {
             predict(rf_fit, data = newdata)$predictions[, "1"]
           })
           rowMeans(preds)
         },
         
         # Fallback
         predict(model_obj, newdata = newdata, type = "response")
  )
}

# -----------------------------------------------------------------------------
# 3 - MAE + Bias fuer eine Modell-Nischen-Kombination
# -----------------------------------------------------------------------------

calc_mae_bias <- function(model_obj, test_df, env, true_suitability_raster, model_type) {
  
  test_df  <- test_df[test_df$split == "test", ]
  env_vals <- terra::extract(env, test_df[, c("x", "y")])[, -1, drop = FALSE]
  newdat   <- cbind(test_df, env_vals)
  
  # NA-Zeilen in Umweltvariablen VOR predict() entfernen
  # → stellt sicher dass pred und true_suit immer gleich lang sind
  env_vars   <- names(env)
  complete   <- complete.cases(newdat[, env_vars])
  
  if (any(!complete)) {
    warning(sprintf(
      "[%s] %d Zeilen mit NA in Umweltvariablen entfernt.",
      model_type, sum(!complete)
    ))
  }
  
  newdat_clean <- newdat[complete, ]
  
  # Vorhersage auf bereinigten Daten
  pred <- tryCatch(
    predict_model(model_obj, newdat_clean, model_type),
    error = function(e) {
      warning(sprintf("[%s] predict() fehlgeschlagen: %s", model_type, e$message))
      rep(NA_real_, nrow(newdat_clean))
    }
  )
  
  # Wahre Suitability auf denselben bereinigten Koordinaten
  true_suit <- terra::extract(true_suitability_raster,
                              newdat_clean[, c("x", "y")])[, -1, drop = TRUE]
  
  # Längen sollten jetzt übereinstimmen - zur Sicherheit prüfen
  if (length(pred) != length(true_suit)) {
    warning(sprintf(
      "[%s] Längen stimmen immer noch nicht überein: pred = %d, true_suit = %d",
      model_type, length(pred), length(true_suit)
    ))
    n_min     <- min(length(pred), length(true_suit))
    pred      <- pred[seq_len(n_min)]
    true_suit <- true_suit[seq_len(n_min)]
  }
  
  ok  <- !is.na(pred) & !is.na(true_suit)
  
  if (sum(ok) == 0) {
    warning(sprintf("[%s] Keine gültigen Wertepaare.", model_type))
    return(data.frame(MAE = NA, Bias = NA, n_test = 0L))
  }
  
  obs <- true_suit[ok]
  prd <- pred[ok]
  
  data.frame(
    MAE    = Metrics::mae(actual = obs, predicted = prd),
    Bias   = Metrics::bias(actual = obs, predicted = prd),
    n_test = sum(ok)
  )
}


# -----------------------------------------------------------------------------
# 4 - Alle Modelltypen x Nischenbreiten durchrechnen
# -----------------------------------------------------------------------------

metrics_combined <- lapply(names(models_all), function(mtype) {
  res <- lapply(niche_names, function(nm) {
    calc_mae_bias(
      model_obj               = models_all[[mtype]][[nm]],
      test_df                 = test_data[[nm]],
      env                     = env,
      true_suitability_raster = true_suitability_rasters[[nm]],
      model_type              = mtype
    )
  })
  res_df <- do.call(rbind, res)
  res_df$niche_breadth <- niche_names
  res_df$model_type    <- mtype
  res_df
}) %>% bind_rows() %>%
  mutate(
    niche_breadth = factor(niche_breadth, levels = niche_names),
    model_type    = factor(model_type, levels = names(models_all))
  )

print(metrics_combined)

# -----------------------------------------------------------------------------
# 5 - Speichern
# -----------------------------------------------------------------------------

dir.create("results", showWarnings = FALSE)
write.csv(metrics_combined, "results/mae_bias_true_suitability.csv", row.names = FALSE)

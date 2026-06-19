# Define packages
library(gbm)
library(terra)
library(predicts)

# load data
pc_raster_masked   <- terra::rast("Data/raster/pc_raster_masked.tif")

narrow       <- readRDS("Data/species/species_narrow.RDS")
intermediate <- readRDS("Data/species/species_intermediate.RDS")
broad        <- readRDS("Data/species/species_broad.RDS")

virtual_narrow_PA       <- narrow[[2]]
virtual_intermediate_PA <- intermediate[[2]]
virtual_broad_PA        <- broad[[2]]

sampling_narrow       <- readRDS("Data/species/sampling_narrow.RDS")
sampling_intermediate <- readRDS("Data/species/sampling_intermediate.RDS")
sampling_broad        <- readRDS("Data/species/sampling_broad.RDS")

brt_runs_narrow       <- sampling_narrow$brt_runs_narrow
brt_runs_intermediate <- sampling_intermediate$brt_runs_intermediate
brt_runs_broad        <- sampling_broad$brt_runs_broad


#### BRT -------------------------------------------------------------------

# ------------------------------------------------------------------
# Narrow - 10 Runs, dann mitteln
# ------------------------------------------------------------------

brt_models_narrow <- list()
brt_preds_narrow  <- list()

for (i in 1:10) {
  
  run_data <- brt_runs_narrow[[i]]
  
  # Umweltvariablen extrahieren
  env_vals <- as.data.frame(terra::extract(
    pc_raster_masked,
    run_data[, c("x", "y")]
  ))[, -1]  # ID-Spalte entfernen
  
  run_data <- cbind(run_data, env_vals)
  run_data <- na.omit(run_data)
  
  # BRT fitten
  brt_model <- gbm::gbm(
    presence ~ PC1 + PC2,
    data              = run_data,
    distribution      = "bernoulli",
    n.trees           = 2000,
    interaction.depth = 3,
    shrinkage         = 0.01,
    bag.fraction      = 0.75,
    cv.folds          = 5
  )
  
  # optimale Baumzahl über CV bestimmen
  best_trees <- gbm::gbm.perf(brt_model, method = "cv", plot.it = FALSE)
  
  # Vorhersage auf ganzes Raster
  pred <- terra::predict(
    pc_raster_masked,
    brt_model,
    n.trees = best_trees,
    type    = "response",
    na.rm   = TRUE
  )
  
  brt_models_narrow[[i]] <- brt_model
  brt_preds_narrow[[i]]  <- pred
}

# Über alle 10 Runs mitteln
pred_brt_narrow <- terra::app(terra::rast(brt_preds_narrow), mean)

plot(pred_brt_narrow, main = "BRT - Narrow Niche")


# ------------------------------------------------------------------
# Intermediate - 10 Runs, dann mitteln
# ------------------------------------------------------------------

brt_models_intermediate <- list()
brt_preds_intermediate  <- list()

for (i in 1:10) {
  
  run_data <- brt_runs_intermediate[[i]]
  
  env_vals <- as.data.frame(terra::extract(
    pc_raster_masked,
    run_data[, c("x", "y")]
  ))[, -1]
  
  run_data <- cbind(run_data, env_vals)
  run_data <- na.omit(run_data)
  
  brt_model <- gbm::gbm(
    presence ~ PC1 + PC2,
    data              = run_data,
    distribution      = "bernoulli",
    n.trees           = 2000,
    interaction.depth = 3,
    shrinkage         = 0.01,
    bag.fraction      = 0.75,
    cv.folds          = 5
  )
  
  best_trees <- gbm::gbm.perf(brt_model, method = "cv", plot.it = FALSE)
  
  pred <- terra::predict(
    pc_raster_masked,
    brt_model,
    n.trees = best_trees,
    type    = "response",
    na.rm   = TRUE
  )
  
  brt_models_intermediate[[i]] <- brt_model
  brt_preds_intermediate[[i]]  <- pred
}

pred_brt_intermediate <- terra::app(terra::rast(brt_preds_intermediate), mean)

plot(pred_brt_intermediate, main = "BRT - Intermediate Niche")


# ------------------------------------------------------------------
# Broad - 10 Runs, dann mitteln
# ------------------------------------------------------------------

brt_models_broad <- list()
brt_preds_broad  <- list()

for (i in 1:10) {
  
  run_data <- brt_runs_broad[[i]]
  
  env_vals <- as.data.frame(terra::extract(
    pc_raster_masked,
    run_data[, c("x", "y")]
  ))[, -1]
  
  run_data <- cbind(run_data, env_vals)
  run_data <- na.omit(run_data)
  
  brt_model <- gbm::gbm(
    presence ~ PC1 + PC2,
    data              = run_data,
    distribution      = "bernoulli",
    n.trees           = 2000,
    interaction.depth = 3,
    shrinkage         = 0.01,
    bag.fraction      = 0.75,
    cv.folds          = 5
  )
  
  best_trees <- gbm::gbm.perf(brt_model, method = "cv", plot.it = FALSE)
  
  pred <- terra::predict(
    pc_raster_masked,
    brt_model,
    n.trees = best_trees,
    type    = "response",
    na.rm   = TRUE
  )
  
  brt_models_broad[[i]] <- brt_model
  brt_preds_broad[[i]]  <- pred
}

pred_brt_broad <- terra::app(terra::rast(brt_preds_broad), mean)

plot(pred_brt_broad, main = "BRT - Broad Niche")


# ------------------------------------------------------------------
# Modelle speichern
# ------------------------------------------------------------------

saveRDS(list(
  models     = brt_models_narrow,
  prediction = pred_brt_narrow
), "models/brt_narrow.RDS")

saveRDS(list(
  models     = brt_models_intermediate,
  prediction = pred_brt_intermediate
), "models/brt_intermediate.RDS")

saveRDS(list(
  models     = brt_models_broad,
  prediction = pred_brt_broad
), "models/brt_broad.RDS")
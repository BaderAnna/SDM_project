# =============================================================================
# RANDOM FOREST MODELLING
# =============================================================================

library(terra)
library(sf)
library(ranger)   # faster than randomForest

# =============================================================================
# 1 - Load raster + variables
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars      <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - Build RF formula dynamically
# =============================================================================

build_rf_formula <- function(vars) {
  as.formula(paste("presence ~", paste(vars, collapse = " + ")))
}

rf_formula <- build_rf_formula(selected_vars)
print(rf_formula)

# =============================================================================
# 3 - Prepare training data (per run)
#     (fold split + pseudo-absence sets already come from train_rf_<sp>.RDS
#      -> here we only need to extract environmental values at the coordinates)
# =============================================================================

prepare_rf_train_data <- function(run_df, env_raster) {
  
  # run_df may additionally contain a "split" column -> use only training data
  if ("split" %in% names(run_df)) {
    run_df <- run_df[run_df$split == "train", , drop = FALSE]
  }
  
  if (nrow(run_df) == 0) stop("run_df is empty!")
  if (any(is.na(run_df$x)) || any(is.na(run_df$y))) {
    stop("x or y contain NA values!")
  }
  
  coords_3035  <- as.matrix(run_df[, c("x", "y")])
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords_3035))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  dat <- cbind(run_df, env_vals)
  dat <- na.omit(dat)
  
  # presence as factor for classification
  dat$presence <- as.factor(dat$presence)
  
  message("prepare_rf_train_data: ", nrow(dat), " rows | ",
          sum(dat$presence == 1), " presence | ",
          sum(dat$presence == 0), " pseudo-absence")
  
  return(dat)
}

# =============================================================================
# 4 - Fit RF on training data
#     10 runs with different pseudo-absences (train_df is the list
#     of 10 runs) -> average predictions
#     (no CV/evaluation - only to guard against pseudo-absence randomness)
# =============================================================================

run_rf_train <- function(train_df, env_raster, rf_formula, seed = 42,
                         n_runs = 10, n_trees = 500) {
  
  n_runs <- min(n_runs, length(train_df))
  
  pred_list  <- list()
  model_list <- list()
  
  for (run in seq_len(n_runs)) {
    
    message(sprintf("  -> Run %d/%d", run, n_runs))
    
    dat <- prepare_rf_train_data(train_df[[run]], env_raster)
    
    if (nrow(dat) == 0) {
      message("  -> Run ", run, ": dat empty -> skipping")
      next
    }
    
    # Check whether all variables in the formula exist in dat
    missing_vars <- setdiff(all.vars(rf_formula), names(dat))
    if (length(missing_vars) > 0) {
      stop("Missing variables in dat: ", paste(missing_vars, collapse = ", "))
    }
    
    model <- withCallingHandlers(
      tryCatch(
        ranger::ranger(
          rf_formula,
          data        = dat,
          num.trees   = n_trees,
          probability = TRUE,
          seed        = seed + run
        ),
        error = function(e) {
          message("  RF error in run ", run, ": ", e$message)
          NULL
        }
      ),
      warning = function(w) {
        message("Warning during model training (run ", run, "): ", w$message)
      }
    )
    if (is.null(model)) next
    
    prediction <- terra::predict(
      env_raster,
      model,
      fun = function(model, data) {
        predict(model, data = data)$predictions[, "1"]
      },
      na.rm = TRUE
    )
    
    model_list[[run]] <- model
    pred_list[[run]]  <- prediction
  }
  
  # --- Average across all runs ---
  message("  -> Averaging across ", length(pred_list), " runs...")
  final_pred <- terra::app(terra::rast(pred_list), mean)
  
  list(
    models     = model_list,
    prediction = final_pred
  )
}

# =============================================================================
# 5 - Run for all species
# =============================================================================

species_list    <- c("narrow", "low_mid", "high_mid", "broad")
rf_models_train  <- list()

for (sp in species_list) {
  
  message("===== ", sp, " (training) =====")
  
  train_rf <- readRDS(paste0("Data/species/split_knndm/train_brt_", sp, ".RDS"))
  
  rf_models_train[[sp]] <- run_rf_train(
    train_df   = train_rf,
    env_raster = env_raster_masked,
    rf_formula = rf_formula,
    n_runs     = 10
  )
}

# =============================================================================
# 6 - Visualisation
# =============================================================================

par(mfrow = c(2, 2))
for (sp in species_list) {
  plot(rf_models_train[[sp]]$prediction, main = paste("RF", sp))
}
par(mfrow = c(1, 1))

# =============================================================================
# 7 - Save training models
# =============================================================================

dir.create("Data/models/knndm", recursive = TRUE, showWarnings = FALSE)

for (sp in species_list) {
  saveRDS(
    rf_models_train[[sp]],
    paste0("Data/models/knndm/rf_", sp, "_train.RDS")
  )
  message("✓ Training model saved: rf_", sp, "_train.RDS")
}
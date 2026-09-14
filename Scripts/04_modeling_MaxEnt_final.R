# =============================================================================
# MAXNET MODELLING
# =============================================================================

library(terra)
library(sf)
library(maxnet)

# =============================================================================
# 1 - Load raster + variables
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars      <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - Prepare training data
#     (fold split already comes from train_maxent_<sp>.RDS -> here we only
#      need to extract environmental values at the coordinates + scale them)
# =============================================================================

prepare_maxnet_train_data <- function(train_df, env_raster) {
  
  if (nrow(train_df) == 0) stop("train_df is empty!")
  if (any(is.na(train_df$x)) || any(is.na(train_df$y))) {
    stop("x or y contain NA values!")
  }
  
  coords_3035  <- as.matrix(train_df[, c("x", "y")])
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords_3035))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  message("NA proportion env_vals: ",
          round(mean(is.na(env_vals)) * 100, 1), "%")
  
  # Standardise variables (z-transformation)
  env_vals_scaled <- as.data.frame(scale(env_vals))
  
  # Save scaling parameters (needed for later prediction!)
  scale_center <- attr(scale(env_vals), "scaled:center")
  scale_scale  <- attr(scale(env_vals), "scaled:scale")
  
  dat <- cbind(train_df, env_vals_scaled)
  dat <- na.omit(dat)
  
  message("prepare_maxnet_train_data: ", nrow(dat), " rows | ",
          sum(dat$presence == 1), " presence | ",
          sum(dat$presence == 0), " pseudo-absence/background")
  
  # Return scaling parameters
  attr(dat, "scale_center") <- scale_center
  attr(dat, "scale_scale")  <- scale_scale
  
  return(dat)
}

# =============================================================================
# 3 - Fit MaxNet on training data
# =============================================================================

run_maxnet_train <- function(train_df, env_raster, seed = 42) {
  
  dat <- prepare_maxnet_train_data(train_df, env_raster)
  
  if (nrow(dat) == 0) stop("dat is empty after prepare_maxnet_train_data()!")
  
  env_cols <- names(env_raster)
  
  # Check whether all environmental variables exist in dat
  missing_vars <- setdiff(env_cols, names(dat))
  if (length(missing_vars) > 0) {
    stop("Missing variables in dat: ", paste(missing_vars, collapse = ", "))
  }
  
  set.seed(seed)
  
  scale_center <- attr(dat, "scale_center")
  scale_scale  <- attr(dat, "scale_scale")
  
  all_env <- dat[, env_cols, drop = FALSE]
  
  model <- withCallingHandlers(
    maxnet(
      p    = dat$presence,
      data = all_env,
      f    = maxnet.formula(dat$presence, all_env, classes = "default")
    ),
    error = function(e) {
      message("Model training failed: ", e$message)
      stop("Model could not be trained.")
    },
    warning = function(w) {
      message("Warning during model training: ", w$message)
    }
  )
  
  # Scale prediction raster to match the training data
  env_raster_scaled <- (env_raster - scale_center) / scale_scale
  names(env_raster_scaled) <- env_cols
  
  prediction <- terra::predict(
    env_raster_scaled,
    model,
    fun    = function(model, data) predict(model, data, type = "cloglog"),
    na.rm  = TRUE
  )
  
  list(
    model        = model,
    prediction   = prediction,
    data         = dat,
    scale_center = scale_center,
    scale_scale  = scale_scale
  )
}

# =============================================================================
# 4 - Run for all species
# =============================================================================

species_list         <- c("narrow", "low_mid", "high_mid", "broad")
maxnet_models_train   <- list()

for (sp in species_list) {
  
  message("===== ", sp, " (training) =====")
  
  train_maxent <- readRDS(paste0("Data/species/split_knndm/train_maxent_", sp, ".RDS"))
  
  maxnet_models_train[[sp]] <- run_maxnet_train(
    train_df   = train_maxent,
    env_raster = env_raster_masked
  )
}

# =============================================================================
# 5 - Visualisation
# =============================================================================

par(mfrow = c(2, 2))
for (sp in species_list) {
  plot(maxnet_models_train[[sp]]$prediction, main = paste("MaxNet", sp))
}
par(mfrow = c(1, 1))

# =============================================================================
# 6 - Save training models
# =============================================================================

dir.create("Data/models/knndm", recursive = TRUE, showWarnings = FALSE)

for (sp in species_list) {
  saveRDS(
    maxnet_models_train[[sp]],
    paste0("Data/models/knndm/maxnet_", sp, "_train.RDS")
  )
  message("✓ Training model saved: maxnet_", sp, "_train.RDS")
}

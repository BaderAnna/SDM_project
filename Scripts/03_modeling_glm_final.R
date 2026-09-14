# =============================================================================
# GLM MODELLING
# =============================================================================

library(terra)
library(sf)

# =============================================================================
# 1 - Load raster + variables
# =============================================================================

env_raster_masked <- terra::rast("Data/raster/env_raster_masked.tif")
selected_vars      <- readRDS("Data/raster/selected_vars.RDS")
names(env_raster_masked) <- selected_vars

# =============================================================================
# 2 - Build GLM formula dynamically
# =============================================================================

build_glm_formula <- function(vars) {
  linear_terms    <- paste(vars, collapse = " + ")
  quadratic_terms <- paste0("I(", vars, "^2)", collapse = " + ")
  as.formula(paste("presence ~", linear_terms, "+", quadratic_terms))
}

glm_formula <- build_glm_formula(selected_vars)
print(glm_formula)

# =============================================================================
# 3 - Prepare training data
#     (fold split already comes from train_glm_<sp>.RDS -> here we only
#      need to extract environmental values at the coordinates)
# =============================================================================

prepare_glm_train_data <- function(train_df, env_raster) {
  
  if (nrow(train_df) == 0) stop("train_df is empty!")
  if (any(is.na(train_df$x)) || any(is.na(train_df$y))) {
    stop("x or y contain NA values!")
  }
  
  coords_3035  <- as.matrix(train_df[, c("x", "y")])
  env_vals_raw <- as.data.frame(terra::extract(env_raster, coords_3035))
  env_vals     <- env_vals_raw[, names(env_raster), drop = FALSE]
  
  dat <- cbind(train_df, env_vals)
  dat <- na.omit(dat)
  
  message("prepare_glm_train_data: ", nrow(dat), " rows | ",
          sum(dat$presence == 1), " presence | ",
          sum(dat$presence == 0), " pseudo-absence")
  
  return(dat)
}

# =============================================================================
# 4 - Fit GLM on training data
# =============================================================================

run_glm_train <- function(train_df, env_raster, glm_formula, seed = 42) {
  
  dat <- prepare_glm_train_data(train_df, env_raster)
  
  if (nrow(dat) == 0) stop("dat is empty after prepare_glm_train_data()!")
  
  # Check whether all variables in the formula exist in dat
  missing_vars <- setdiff(all.vars(glm_formula), names(dat))
  if (length(missing_vars) > 0) {
    stop("Missing variables in dat: ", paste(missing_vars, collapse = ", "))
  }
  
  set.seed(seed)
  
  # Weights: presences = 1, pseudo-absences = n_pres / n_abs
  n_pres <- sum(dat$presence == 1)
  n_abs  <- sum(dat$presence == 0)
  dat$.weights <- ifelse(dat$presence == 1, 1, n_pres / n_abs)
  
  model <- withCallingHandlers(
    glm(glm_formula, data = dat, family = binomial, weights = .weights),
    error = function(e) {
      message("Model training failed: ", e$message)
      stop("Model could not be trained.")
    },
    warning = function(w) {
      message("Warning during model training: ", w$message)
    }
  )
  
  prediction <- terra::predict(env_raster, model, type = "response")
  
  list(model = model, prediction = prediction, data = dat)
}

# =============================================================================
# 5 - Run for all species
# =============================================================================

species_list     <- c("narrow", "low_mid", "high_mid", "broad")
glm_models_train <- list()

for (sp in species_list) {
  
  message("===== ", sp, " (training) =====")
  
  train_glm <- readRDS(paste0("Data/species/split_knndm/train_glm_", sp, ".RDS"))
  
  glm_models_train[[sp]] <- run_glm_train(
    train_df    = train_glm,
    env_raster  = env_raster_masked,
    glm_formula = glm_formula
  )
}

# =============================================================================
# 6 - Save training models
# =============================================================================

dir.create("Data/models/knndm", recursive = TRUE, showWarnings = FALSE)

for (sp in species_list) {
  saveRDS(
    glm_models_train[[sp]],
    paste0("Data/models/knndm/glm_", sp, "_train.RDS")
  )
  message("✓ Training model saved: glm_", sp, "_train.RDS")
}
# =============================================================================
# 7.1 Model evaluation TSS + AUC
# =============================================================================

library(mecofun)
library(PresenceAbsence)
library(terra)
library(dplyr)
library(maxnet)
library(gbm)
library(ranger)

dir.create("results", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1 - Read environmental raster
# -----------------------------------------------------------------------------
env <- terra::rast("Data/raster/env_raster_masked.tif")
print(names(env))

# -----------------------------------------------------------------------------
# 2 - Load models and test data
# -----------------------------------------------------------------------------
niche_names <- c("narrow", "low_mid", "high_mid", "broad")

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

models_all <- list(GLM = models_glm, Maxent = models_maxent,
                   BRT = models_brt, RF = models_rf)

test_data <- list(
  narrow   = readRDS("Data/species/split_knndm/test_narrow.RDS"),
  low_mid  = readRDS("Data/species/split_knndm/test_low_mid.RDS"),
  high_mid = readRDS("Data/species/split_knndm/test_high_mid.RDS"),
  broad    = readRDS("Data/species/split_knndm/test_broad.RDS")
)

# -----------------------------------------------------------------------------
# 3 - Predict function
# -----------------------------------------------------------------------------
predict_model <- function(model_obj, newdata, model_type) {
  switch(model_type,
         
         "GLM" = predict(model_obj$model, newdata = newdata, type = "response"),
         
         "Maxent" = {
           scale_center <- model_obj$scale_center
           scale_scale  <- model_obj$scale_scale
           env_vars     <- names(scale_center)
           
           newdata_scaled <- newdata
           newdata_scaled[, env_vars] <- scale(newdata[, env_vars],
                                               center = scale_center,
                                               scale  = scale_scale)
           
           as.vector(predict(model_obj$model, newdata = newdata_scaled, type = "cloglog"))
         },
         
         "BRT" = {
           gbm_models <- model_obj$models
           n_trees    <- gbm_models[[1]]$n.trees
           preds <- sapply(gbm_models, function(gbm_fit) {
             gbm::predict.gbm(gbm_fit, newdata = newdata, n.trees = n_trees, type = "response")
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
         
         predict(model_obj, newdata = newdata, type = "response")
  )
}

# -----------------------------------------------------------------------------
# 4 - Evaluation function
# -----------------------------------------------------------------------------
eval_one_niche <- function(model_obj, test_df, env,
                           model_type = "GLM", thresh.method = "ObsPrev") {
  
  test_df  <- test_df[test_df$split == "test", ]
  env_vals <- terra::extract(env, test_df[, c("x", "y")])[, -1, drop = FALSE]
  newdat   <- cbind(test_df, env_vals)
  
  pred <- tryCatch(
    predict_model(model_obj, newdat, model_type),
    error = function(e) {
      warning(sprintf("[%s] predict() failed: %s", model_type, e$message))
      rep(NA_real_, nrow(newdat))
    }
  )
  
  # If rows of predict() were deleted internally: fill with NA
  if (length(pred) != nrow(newdat)) {
    env_cols  <- names(env_vals)
    na_rows   <- apply(newdat[, env_cols, drop = FALSE], 1, anyNA)
    pred_full <- rep(NA_real_, nrow(newdat))
    pred_full[!na_rows] <- pred
    pred <- pred_full
  }
  
  ok <- !is.na(pred)
  if (any(!ok)) {
    warning(sprintf("[%s] %d of %d test points excluded (NA).",
                    model_type, sum(!ok), length(ok)))
  }
  
  if (length(unique(newdat$presence[ok])) < 2) {
    warning(sprintf("[%s] Only one class in the test dataset – evaluation skipped.", model_type))
    return(data.frame(AUC = NA, TSS = NA, Kappa = NA, Sens = NA, Spec = NA,
                      PCC = NA, D2 = NA, thresh = NA,
                      true_prevalence = NA, n_test = sum(ok)))
  }
  
  eval_res <- evalSDM(observation = as.vector(newdat$presence[ok]),
                      predictions = as.vector(pred[ok]),
                      thresh.method = thresh.method)
  
  eval_res$true_prevalence <- mean(newdat$presence[ok])
  eval_res$n_test          <- sum(ok)
  eval_res
}

# -----------------------------------------------------------------------------
# 5 - Bootstrap-CI of TSS
# -----------------------------------------------------------------------------
bootstrap_tss <- function(obs, pred, thresh.method = "ObsPrev", n_boot = 1000) {
  n        <- length(obs)
  tss_boot <- numeric(n_boot)
  
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    if (length(unique(obs[idx])) < 2) next
    res <- tryCatch(
      evalSDM(observation = obs[idx], predictions = pred[idx], thresh.method = thresh.method),
      error = function(e) NULL
    )
    tss_boot[i] <- if (!is.null(res)) res$TSS else NA
  }
  quantile(tss_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# 6 - Main loop
# -----------------------------------------------------------------------------
all_results <- list()
all_ci      <- list()

for (mtype in names(models_all)) {
  
  cat("\n=== Evaluate:", mtype, "===\n")
  models_mtype <- models_all[[mtype]]
  
  res_list <- Map(function(mod, tdat) eval_one_niche(mod, tdat, env, mtype, "ObsPrev"),
                  models_mtype, test_data)
  
  res_df <- dplyr::bind_rows(res_list)
  res_df$niche_breadth <- niche_names
  res_df$model_type    <- mtype
  all_results[[mtype]] <- res_df
  
  ci_list <- lapply(niche_names, function(nm) {
    test_df  <- test_data[[nm]][test_data[[nm]]$split == "test", ]
    env_vals <- terra::extract(env, test_df[, c("x", "y")])[, -1, drop = FALSE]
    newdat   <- cbind(test_df, env_vals)
    pred     <- tryCatch(predict_model(models_mtype[[nm]], newdat, mtype),
                         error = function(e) rep(NA_real_, nrow(newdat)))
    if (length(pred) != nrow(newdat)) {
      env_cols  <- names(env_vals)
      na_rows   <- apply(newdat[, env_cols, drop = FALSE], 1, anyNA)
      pred_full <- rep(NA_real_, nrow(newdat))
      pred_full[!na_rows] <- pred
      pred <- pred_full
    }
    ok <- !is.na(pred)
    bootstrap_tss(newdat$presence[ok], pred[ok])
  })
  names(ci_list) <- niche_names
  all_ci[[mtype]] <- ci_list
}

# -----------------------------------------------------------------------------
# 7 - Combine and save results 
# -----------------------------------------------------------------------------
results_combined <- dplyr::bind_rows(all_results)

ci_df <- do.call(rbind, lapply(names(all_ci), function(mtype) {
  do.call(rbind, lapply(names(all_ci[[mtype]]), function(nm) {
    ci <- all_ci[[mtype]][[nm]]
    data.frame(model_type = mtype, niche_breadth = nm,
               TSS_lo = ci["2.5%"], TSS_med = ci["50%"], TSS_hi = ci["97.5%"],
               row.names = NULL)
  }))
}))

presence_overview <- data.frame(
  niche_breadth   = niche_names,
  n_presence_test = sapply(test_data, function(d) sum(d$presence[d$split == "test"] == 1)),
  n_absence_test  = sapply(test_data, function(d) sum(d$presence[d$split == "test"] == 0))
)
presence_overview$prevalence_test <- presence_overview$n_presence_test /
  (presence_overview$n_presence_test + presence_overview$n_absence_test)

write.csv(results_combined,  "results/all_models_auc_tss.csv", row.names = FALSE)
write.csv(ci_df,             "results/all_models_tss_ci.csv",  row.names = FALSE)
write.csv(presence_overview, "results/presence_overview.csv",  row.names = FALSE)
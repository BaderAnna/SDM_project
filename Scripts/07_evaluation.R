#GLM
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

# MaxEnt 
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
# =============================================================================
# GLM MODELLIERUNG mit 10x Spatial 5-fold CV
# =============================================================================

library(terra)
library(blockCV)
library(sf)
# pROC installieren falls nötig
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

# --- Daten laden ---
pc_raster_masked <- terra::rast("Data/raster/pc_raster_masked.tif")
# Raster in metrisches KBS transformieren
#pc_raster_3035 <- terra::project(pc_raster_masked, "EPSG:3035")

samp_narrow   <- readRDS("Data/species/sampling_narrow.RDS")
samp_low_mid  <- readRDS("Data/species/sampling_low_mid.RDS")
samp_high_mid <- readRDS("Data/species/sampling_high_mid.RDS")
samp_broad    <- readRDS("Data/species/sampling_broad.RDS")

# =============================================================================
# Helper-Funktion: Daten vorbereiten
# =============================================================================
prepare_glm_data <- function(sampling, pc_raster) {
  
  # Presences
  pres <- sampling$po$sample.points[
    sampling$po$sample.points$Observed == 1, c("x", "y")]
  pres$presence <- 1
  
  # Absences
  abs <- sampling$bg_glm
  abs$presence <- 0
  
  # Kombinieren
  dat <- rbind(pres, abs)
  
  # Umweltvariablen extrahieren
  env_vals <- as.data.frame(terra::extract(
    pc_raster, dat[, c("x", "y")]
  ))[, -1]
  dat <- cbind(dat, env_vals)
  
  # NAs entfernen
  dat <- na.omit(dat)
  
  return(dat)
}

# =============================================================================
# BlockCV: Block Size herausfinden über Autokorrelation
# =============================================================================
set.seed(42)
# Raster in metrisches KBS transformieren
#pc_raster_3035 <- terra::project(pc_raster_masked, "EPSG:3035")

# Dann Autokorrelation berechnen
'autocor <- blockCV::cv_spatial_autocor(
  r          = pc_raster_3035,
  num_sample = 5000
)'

block_size <- 100000

# =============================================================================
# Helper-Funktion: GLM mit 10x Spatial 5-fold CV
# =============================================================================

run_glm_cv <- function(sampling, pc_raster, n_rep = 10, n_folds = 5) {
  
  dat <- prepare_glm_data(sampling, pc_raster)
  dat_sf <- sf::st_as_sf(dat, coords = c("x", "y"), crs = 4326)
  
  auc_scores  <- c()
  tss_scores  <- c()
  
  for (rep in 1:n_rep) {
    
    set.seed(rep)
    blocks <- blockCV::cv_spatial(
      x           = dat_sf,
      column      = "presence",
      k           = 5,
      hexagon     = TRUE,        # Hexagone wie im Kurs
      size        = block_size,  # räumliche Autokorrelationsreichweite
      selection   = "random",
      iteration   = 100L,
      presence_bg = TRUE,
      extend      = 0.5,
      report      = FALSE,
      plot        = FALSE
    )
    
    fold_preds  <- rep(NA, nrow(dat))
    fold_truths <- rep(NA, nrow(dat))
    
    for (fold in 1:n_folds) {
      
      train_idx <- blocks$folds_list[[fold]][[1]]
      test_idx  <- blocks$folds_list[[fold]][[2]]
      
      train <- dat[train_idx, ]
      test  <- dat[test_idx, ]
      
      # Gewichtung
      n_pres  <- sum(train$presence == 1)
      n_abs   <- sum(train$presence == 0)
      w_train <- ifelse(train$presence == 1, 1, n_pres / n_abs)
      
      # GLM fitten
      model <- suppressWarnings(glm(
        presence ~ PC1 + PC2 + I(PC1^2) + I(PC2^2),
        data    = train,
        family  = binomial,
        weights = w_train
      ))
      
      # Vorhersage auf Testdaten
      fold_preds[test_idx]  <- predict(model, newdata = test, type = "response")
      fold_truths[test_idx] <- test$presence
    }
    
    # Nur Punkte mit Vorhersage behalten
    valid <- !is.na(fold_preds) & !is.na(fold_truths)
    preds_valid  <- fold_preds[valid]
    truths_valid <- fold_truths[valid]
    
    # AUC berechnen
    auc <- as.numeric(pROC::auc(
      pROC::roc(truths_valid, preds_valid, quiet = TRUE)
    ))
    
    # TSS berechnen - Fix 2: direkt über Confusion Matrix
    threshold <- 0.5  # fixer Schwellenwert
    pred_bin  <- ifelse(preds_valid >= threshold, 1, 0)
    
    tp <- sum(pred_bin == 1 & truths_valid == 1)
    fp <- sum(pred_bin == 1 & truths_valid == 0)
    tn <- sum(pred_bin == 0 & truths_valid == 0)
    fn <- sum(pred_bin == 0 & truths_valid == 1)
    
    sens <- tp / (tp + fn)  # Sensitivität
    spec <- tn / (tn + fp)  # Spezifität
    tss  <- sens + spec - 1
    
    auc_scores[rep] <- auc
    tss_scores[rep] <- tss
    
    message("Rep ", rep, "/", n_rep, 
            " | AUC = ", round(auc, 3),
            " | TSS = ", round(tss, 3))
  }
  
  # Finales Modell auf allen Daten
  n_pres <- sum(dat$presence == 1)
  n_abs  <- sum(dat$presence == 0)
  w_all  <- ifelse(dat$presence == 1, 1, n_pres / n_abs)
  
  final_model <- suppressWarnings(glm(
    presence ~ PC1 + PC2 + I(PC1^2) + I(PC2^2),
    data    = dat,
    family  = binomial,
    weights = w_all
  ))
  
  final_pred <- terra::predict(pc_raster, final_model, type = "response")
  
  return(list(
    model      = final_model,
    prediction = final_pred,
    data       = dat,
    auc_mean   = mean(auc_scores),
    auc_sd     = sd(auc_scores),
    tss_mean   = mean(tss_scores),
    tss_sd     = sd(tss_scores),
    auc_all    = auc_scores,
    tss_all    = tss_scores
  ))
}

# =============================================================================
# GLM für alle 4 Arten
# =============================================================================

glm_narrow   <- run_glm_cv(samp_narrow,   pc_raster_3035)

glm_low_mid  <- run_glm_cv(samp_low_mid,  pc_raster_3035)

glm_high_mid <- run_glm_cv(samp_high_mid, pc_raster_3035)

glm_broad    <- run_glm_cv(samp_broad,    pc_raster_3035)

# =============================================================================
# Ergebnisse zusammenfassen
# =============================================================================

results_glm <- data.frame(
  species  = c("narrow", "low_mid", "high_mid", "broad"),
  sigma    = c(0.2, 0.4, 0.6, 0.8),
  auc_mean = c(glm_narrow$auc_mean, glm_low_mid$auc_mean,
               glm_high_mid$auc_mean, glm_broad$auc_mean),
  auc_sd   = c(glm_narrow$auc_sd, glm_low_mid$auc_sd,
               glm_high_mid$auc_sd, glm_broad$auc_sd),
  tss_mean = c(glm_narrow$tss_mean, glm_low_mid$tss_mean,
               glm_high_mid$tss_mean, glm_broad$tss_mean),
  tss_sd   = c(glm_narrow$tss_sd, glm_low_mid$tss_sd,
               glm_high_mid$tss_sd, glm_broad$tss_sd)
)

print(results_glm)

# =============================================================================
# Visualisierung
# =============================================================================

par(mfrow = c(2, 2))
plot(glm_narrow$prediction,   main = paste0("GLM - Narrow (AUC = ",
                                            round(glm_narrow$auc_mean, 3), ")"))
plot(glm_low_mid$prediction,  main = paste0("GLM - Low-Mid (AUC = ",
                                            round(glm_low_mid$auc_mean, 3), ")"))
plot(glm_high_mid$prediction, main = paste0("GLM - High-Mid (AUC = ",
                                            round(glm_high_mid$auc_mean, 3), ")"))
plot(glm_broad$prediction,    main = paste0("GLM - Broad (AUC = ",
                                            round(glm_broad$auc_mean, 3), ")"))
par(mfrow = c(1, 1))

# =============================================================================
# Speichern
# =============================================================================

dir.create("Data/models", recursive = TRUE, showWarnings = FALSE)

saveRDS(glm_narrow,   "Data/models/glm_narrow.RDS")
saveRDS(glm_low_mid,  "Data/models/glm_low_mid.RDS")
saveRDS(glm_high_mid, "Data/models/glm_high_mid.RDS")
saveRDS(glm_broad,    "Data/models/glm_broad.RDS")

saveRDS(results_glm,  "Data/models/results_glm.RDS")

print(results_glm)






'# Define packages
list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "geodata")

# Load data
pc_raster_masked   <- terra::rast("Data/raster/pc_raster_masked.tif")

narrow      <- readRDS("Data/species/species_narrow.RDS")
intermediate <- readRDS("Data/species/species_intermediate.RDS")
broad       <- readRDS("Data/species/species_broad.RDS")

virtual_narrow_PA      <- narrow[[2]]
virtual_intermediate_PA <- intermediate[[2]]
virtual_broad_PA       <- broad[[2]]

sampling_narrow      <- readRDS("Data/species/sampling_narrow.RDS")
sampling_intermediate <- readRDS("Data/species/sampling_intermediate.RDS")
sampling_broad       <- readRDS("Data/species/sampling_broad.RDS")

#### GLM ---------------------------

# Objekte extrahieren
po_narrow       <- sampling_narrow$po_narrow
bg_narrow_glm   <- sampling_narrow$bg_narrow_glm
pres_narrow     <- sampling_narrow$pres_narrow
brt_runs_narrow <- sampling_narrow$brt_runs_narrow

po_intermediate       <- sampling_intermediate$po_intermediate
bg_intermediate_glm   <- sampling_intermediate$bg_intermediate_glm
pres_intermediate     <- sampling_intermediate$pres_intermediate
brt_runs_intermediate <- sampling_intermediate$brt_runs_intermediate

po_broad       <- sampling_broad$po_broad
bg_broad_glm   <- sampling_broad$bg_broad_glm
pres_broad     <- sampling_broad$pres_broad
brt_runs_broad <- sampling_broad$brt_runs_broad

##### narrow niche ----------------------------------------------------------

# 1. Presences zusammenstellen
pres_narrow_glm <- po_narrow$sample.points[
  po_narrow$sample.points$Observed == 1, c("x", "y")]
pres_narrow_glm$presence <- 1

# 2. Absences zusammenstellen
abs_narrow_glm <- as.data.frame(bg_narrow_glm)
abs_narrow_glm$presence <- 0

# 3. Kombinieren
glm_data_narrow <- rbind(pres_narrow_glm, abs_narrow_glm)

# 4. Umweltvariablen extrahieren
env_vals <- as.data.frame(terra::extract(
  pc_raster_masked,                        # <- masked verwenden!
  glm_data_narrow[, c("x", "y")]
))
glm_data_narrow <- cbind(glm_data_narrow, env_vals)

# 5. NAs entfernen
glm_data_narrow <- na.omit(glm_data_narrow)

# 6. Gewichtung
weights_narrow <- c(rep(1,   nrow(pres_narrow_glm)),   # Presences
                    rep(0.1, nrow(abs_narrow_glm)))     # Absences

# 7. GLM fitten
glm_narrow <- glm(
  presence ~ PC1 + PC2 + I(PC1^2) + I(PC2^2),
  data    = glm_data_narrow,
  family  = binomial,
  weights = weights_narrow
)

summary(glm_narrow)

# Vorhersage auf ganzes Raster
pred_glm_narrow <- terra::predict(
  pc_raster_masked,
  glm_narrow,
  type = "response"
)

# Visualisieren
plot(pred_glm_narrow, main = "GLM - Narrow Niche")


##### intermediate niche ----------------------------------------------------------

# 1. Presences zusammenstellen
pres_intermediate_glm <- po_intermediate$sample.points[
  po_intermediate$sample.points$Observed == 1, c("x", "y")]
pres_intermediate_glm$presence <- 1

# 2. Absences zusammenstellen
abs_intermediate_glm <- as.data.frame(bg_intermediate_glm)
abs_intermediate_glm$presence <- 0

# 3. Kombinieren
glm_data_intermediate <- rbind(pres_intermediate_glm, abs_intermediate_glm)

# 4. Umweltvariablen extrahieren
env_vals <- as.data.frame(terra::extract(
  pc_raster_masked,                        
  glm_data_intermediate[, c("x", "y")]
))
glm_data_intermediate <- cbind(glm_data_intermediate, env_vals)

# 5. NAs entfernen
glm_data_intermediate <- na.omit(glm_data_intermediate)

# 6. Gewichtung
weights_intermediate <- c(rep(1,   nrow(pres_intermediate_glm)),   # Presences
                    rep(0.1, nrow(abs_intermediate_glm)))     # Absences

# 7. GLM fitten
glm_intermediate <- glm(
  presence ~ PC1 + PC2 + I(PC1^2) + I(PC2^2),
  data    = glm_data_intermediate,
  family  = binomial,
  weights = weights_intermediate
)

summary(glm_intermediate)

# Vorhersage auf ganzes Raster
pred_glm_intermediate <- terra::predict(
  pc_raster_masked,
  glm_intermediate,
  type = "response"
)

# Visualisieren
plot(pred_glm_intermediate, main = "GLM - Intermediate Niche")

##### broad niche ----------------------------------------------------------

# 1. Presences zusammenstellen
pres_broad_glm <- po_broad$sample.points[
  po_broad$sample.points$Observed == 1, c("x", "y")]
pres_broad_glm$presence <- 1

# 2. Absences zusammenstellen
abs_broad_glm <- as.data.frame(bg_broad_glm)
abs_broad_glm$presence <- 0

# 3. Kombinieren
glm_data_broad <- rbind(pres_broad_glm, abs_broad_glm)

# 4. Umweltvariablen extrahieren
env_vals <- as.data.frame(terra::extract(
  pc_raster_masked,
  glm_data_broad[, c("x", "y")]
))
glm_data_broad <- cbind(glm_data_broad, env_vals)

# 5. NAs entfernen
glm_data_broad <- na.omit(glm_data_broad)

# 6. Gewichtung
weights_broad <- c(rep(1,   nrow(pres_broad_glm)),
                   rep(0.1, nrow(abs_broad_glm)))

# 7. GLM fitten
glm_broad <- glm(
  presence ~ PC1 + PC2 + I(PC1^2) + I(PC2^2),
  data    = glm_data_broad,
  family  = binomial,
  weights = weights_broad
)

summary(glm_broad)

# Vorhersage auf ganzes Raster
pred_glm_broad <- terra::predict(
  pc_raster_masked,
  glm_broad,
  type = "response"
)

# Visualisieren
plot(pred_glm_broad, main = "GLM - Broad Niche")

##### Comparisons----------------------------------------------------------

par(mfrow = c(1, 3))
plot(pred_glm_narrow,       main = "GLM - Narrow (σ = 0.2)")
plot(pred_glm_intermediate, main = "GLM - Intermediate (σ = 0.5)")
plot(pred_glm_broad,        main = "GLM - Broad (σ = 0.8)")
par(mfrow = c(1, 1))


#### GLM Modelle und Vorhersagen speichern-----------------------
saveRDS(list(
  model      = glm_narrow,
  prediction = pred_glm_narrow,
  data       = glm_data_narrow
), "models/glm_narrow.RDS")

saveRDS(list(
  model      = glm_intermediate,
  prediction = pred_glm_intermediate,
  data       = glm_data_intermediate
), "models/glm_intermediate.RDS")

saveRDS(list(
  model      = glm_broad,
  prediction = pred_glm_broad,
  data       = glm_data_broad
), "models/glm_broad.RDS")


# Sind die Koeffizienten signifikant?
summary(glm_narrow)

# Ist das Modell besser als ein Nullmodell?
anova(glm_narrow, test = "Chisq")
'
library(terra)
library(gbm)
library(ranger)
library(maxnet)
library(ggplot2)

# Modelle laden
glm     <- readRDS("Data/models/glm_low_mid.RDS")
brt     <- readRDS("Data/models/brt_low_mid.RDS")
rf      <- readRDS("Data/models/rf_low_mid.RDS")
maxent  <- readRDS("Data/models/maxnet_low_mid.RDS")

# Prädiktoren laden 
selected_vars <- readRDS("Data/raster/selected_vars.RDS")  
all_predictors <- rast("Data/raster/env_raster_masked.tif")

# Nur die ausgewählten Variablen behalten
predictors <- all_predictors[[selected_vars]]

str(brt$models, max.level = 1)
length(brt$models)

str(rf$models, max.level = 1)
length(rf$models)

# Vorhersagen berechnen
pred_glm    <- predict(predictors, glm$model, type = "response")
# BRT: über alle Modelle mitteln
pred_brt_list <- lapply(brt$models, function(m) {
  predict(predictors, m, n.trees = m$gbm.call$best.trees, type = "response")
})
pred_brt <- mean(rast(pred_brt_list))

# RF: über alle Modelle mitteln
ranger_predict_fun <- function(model, data, ...) {
  predict(model, data = data, ...)$predictions[, 2]  # Spalte 2 = Presence-Wahrscheinlichkeit
}

pred_rf_list <- lapply(rf$models, function(m) {
  predict(predictors, m, fun = ranger_predict_fun, na.rm = TRUE)
})
pred_rf <- mean(rast(pred_rf_list))

# Predictors mit den Trainings-Skalierungsparametern standardisieren
predictors_scaled <- (predictors - maxent$scale_center) / maxent$scale_scale
names(predictors_scaled) <- names(predictors)  # Namen erhalten

maxnet_predict_fun <- function(model, data, ...) {
  predict(model, data, type = "cloglog", clamp = TRUE)[, 1]
}

pred_maxent <- predict(predictors_scaled, maxent$model,
                       fun = maxnet_predict_fun, na.rm = TRUE)


# Zu einem Stack zusammenfassen und nebeneinander plotten
mask_ref <- predictors[[1]]  # Referenz-Layer mit NA außerhalb Deutschlands

pred_glm    <- mask(pred_glm, mask_ref)
pred_brt    <- mask(pred_brt, mask_ref)
pred_rf     <- mask(pred_rf, mask_ref)
pred_maxent <- mask(pred_maxent, mask_ref)

preds <- c(pred_glm, pred_maxent, pred_brt, pred_rf)
names(preds) <- c("GLM", "MaxEnt", "BRT", "RF")

# In langes Format für facettierten ggplot bringen
preds_df <- as.data.frame(preds, xy = TRUE) |>
  tidyr::pivot_longer(cols = c(GLM, MaxEnt, BRT, RF),
                      names_to = "Model", values_to = "Suitability") |>
  dplyr::mutate(Model = factor(Model, levels = c("GLM", "MaxEnt", "BRT", "RF")))

ggplot(preds_df, aes(x = x, y = y, fill = Suitability)) +
  geom_raster() +
  facet_wrap(~ Model, nrow = 1) +
  scale_fill_viridis_c(
    name = "Habitat-\neignung",
    limits = c(0, 1),        # gemeinsame Skala über alle Modelle
    na.value = "transparent"
  ) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA)
  )


ggplot(preds_df, aes(x = x, y = y, fill = Suitability)) +
  geom_raster() +
  facet_wrap(~ Model, nrow = 1) +
  scale_fill_viridis_c(limits = c(0, 1), na.value = "transparent") +
  coord_equal() +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

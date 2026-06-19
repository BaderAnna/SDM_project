# Define packages
library(maxnet)
library(terra)
library(predicts)

# load data
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

po_narrow <- sampling_narrow$po_narrow
po_intermediate <- sampling_intermediate$po_intermediate
po_broad <- sampling_broad$po_broad

#### MaxEnt -------------------------------------------------------------------

# ------------------------------------------------------------------
# Narrow
# ------------------------------------------------------------------
# Presence-Punkte extrahieren
pres_narrow <- po_narrow$sample.points[
  po_narrow$sample.points$Observed == 1,
  c("x", "y")
]

# Background-Punkte erzeugen
bg_narrow <- predicts::backgroundSample(
  pc_raster_masked,
  n = 10000
)

# Daten vorbereiten 
# Presence = 1
pres <- pres_narrow

# Background = 0
bg <- as.data.frame(bg_narrow)

pres$presence <- 1
bg$presence <- 0

dat <- rbind(pres, bg)

# Umweltvariablen extrahieren
env <- terra::extract(
  pc_raster_masked,
  dat[, c("x", "y")]
)

env <- env[, -1]

# fit model
mx <- maxnet(
  p = dat$presence,
  data = env,
  f = maxnet.formula(
    dat$presence,
    env,
    classes = "default"
  )
)

# Vorhersage
pred_narrow <- terra::predict(
  pc_raster_masked,
  mx,
  fun = function(model, data) {
    predict(model, data, type = "cloglog")
  },
  na.rm = TRUE
)

plot(pred_narrow)

# ------------------------------------------------------------------
# Intermediate
# ------------------------------------------------------------------
# Presence-Punkte extrahieren
pres_intermediate <- po_intermediate$sample.points[
  po_intermediate$sample.points$Observed == 1,
  c("x", "y")
]

# Background-Punkte erzeugen
bg_intermediate <- predicts::backgroundSample(
  pc_raster_masked,
  n = 10000
)

# Daten vorbereiten 
# Presence = 1
pres_intermediate <- pres_intermediate

# Background = 0
bg_intermediate <- as.data.frame(bg_intermediate)

pres_intermediate$presence <- 1
bg_intermediate$presence <- 0

dat_intermediate <- rbind(pres_intermediate, bg_intermediate)

# Umweltvariablen extrahieren
env_intermediate <- terra::extract(
  pc_raster_masked,
  dat_intermediate[, c("x", "y")]
)

env_intermediate <- env_intermediate[, -1]

# fit model
mx_intermediate <- maxnet(
  p = dat_intermediate$presence,
  data = env_intermediate,
  f = maxnet.formula(
    dat_intermediate$presence,
    env_intermediate,
    classes = "default"
  )
)

# Vorhersage
pred_intermediate <- terra::predict(
  pc_raster_masked,
  mx_intermediate,
  fun = function(model, data) {
    predict(model, data, type = "cloglog")
  },
  na.rm = TRUE
)

plot(pred_intermediate)

# ------------------------------------------------------------------
# Broad
# ------------------------------------------------------------------

# Presence-Punkte extrahieren
pres_broad <- po_broad$sample.points[
  po_broad$sample.points$Observed == 1,
  c("x", "y")
]

# Background-Punkte erzeugen
bg_broad <- predicts::backgroundSample(
  pc_raster_masked,
  n = 10000
)

# Daten vorbereiten
# Presence = 1
pres_broad$presence <- 1

# Background = 0
bg_broad <- as.data.frame(bg_broad)
bg_broad$presence <- 0

dat_broad <- rbind(pres_broad, bg_broad)

# Umweltvariablen extrahieren
env_broad <- terra::extract(
  pc_raster_masked,
  dat_broad[, c("x", "y")]
)

env_broad <- env_broad[, -1]

# Modell fitten
mx_broad <- maxnet(
  p = dat_broad$presence,
  data = env_broad,
  f = maxnet.formula(
    dat_broad$presence,
    env_broad,
    classes = "default"
  )
)

# Vorhersage
pred_broad <- terra::predict(
  pc_raster_masked,
  mx_broad,
  fun = function(model, data) {
    predict(model, data, type = "cloglog")
  },
  na.rm = TRUE
)

plot(pred_broad)

# ------------------------------------------------------------------
# Modelle speichern
# ------------------------------------------------------------------

saveRDS(list(
  model      = mx,
  prediction = pred_narrow,
  data       = dat
), "models/maxent_narrow.RDS")

saveRDS(list(
  model      = mx_intermediate,
  prediction = pred_intermediate,
  data       = dat_intermediate
), "models/maxent_intermediate.RDS")

saveRDS(list(
  model      = mx_broad,
  prediction = pred_broad,
  data       = dat_broad
), "models/maxent_broad.RDS")

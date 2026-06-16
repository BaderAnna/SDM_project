# Define packages
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

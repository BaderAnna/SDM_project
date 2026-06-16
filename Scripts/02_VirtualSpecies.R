# Define packages
list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "geodata")

# Daten einlesen
pc_raster_masked   <- terra::rast("Data/raster/pc_raster_masked.tif")

# VIRTUELLE ARTEN ERSTELLEN ####
#-------------------------------------------#
# Gaussian response functions entlang PC1 und PC2
# σ-Werte: 0.2 (schmal), 0.5 (mittel), 0.8 (breit)

set.seed(42)

# --- Narrow niche (specialist) ---
params_narrow <- formatFunctions(
  PC1 = c(fun = "dnorm", mean = 0, sd = 0.2),
  PC2 = c(fun = "dnorm", mean = 0, sd = 0.2)
)

virtual_narrow <- generateSpFromFun(
  raster.stack = pc_raster_masked,
  parameters   = params_narrow,
  species.type = "multiplicative",
  plot         = TRUE
)
virtual_narrow$species.name <- "narrow"
virtual_narrow_PA <- convertToPA(
  virtual_narrow, 
  beta = 0.5, 
  alpha = -0.05, 
  plot = TRUE
)


# --- Intermediate niche ---
params_intermediate <- formatFunctions(
  PC1 = c(fun = "dnorm", mean = 0, sd = 0.5),
  PC2 = c(fun = "dnorm", mean = 0, sd = 0.5)
)

virtual_intermediate <- generateSpFromFun(
  raster.stack = pc_raster_masked,
  parameters   = params_intermediate,
  species.type = "multiplicative",
  plot         = TRUE
)
virtual_intermediate$species.name <- "intermediate"
virtual_intermediate_PA <- convertToPA(
  virtual_intermediate, 
  beta = 0.5, 
  alpha = -0.05, 
  plot = TRUE
)

# --- Broad niche (generalist) ---
params_broad <- formatFunctions(
  PC1 = c(fun = "dnorm", mean = 0, sd = 0.8),
  PC2 = c(fun = "dnorm", mean = 0, sd = 0.8)
)

virtual_broad <- generateSpFromFun(
  raster.stack = pc_raster_masked,
  parameters   = params_broad,
  species.type = "multiplicative",
  plot         = TRUE
)
virtual_broad$species.name <- "broad"
virtual_broad_PA <- convertToPA(
  virtual_broad, 
  beta = 0.5, 
  alpha = -0.05, 
  plot = TRUE
)

# --- Presence Points sampeln ---
set.seed(42)

pa_narrow <- sampleOccurrences(
  virtual_narrow_PA,
  n = 100,
  type = "presence only",
  correct.by.suitability = TRUE
)

pa_intermediate <- sampleOccurrences(
  virtual_intermediate_PA,
  n = 100,
  type = "presence only",
  correct.by.suitability = TRUE
)

pa_broad <- sampleOccurrences(
  virtual_broad_PA,
  n = 100,
  type = "presence only",
  correct.by.suitability = TRUE
)

# --- Speichern ---
saveRDS(list(virtual_narrow, virtual_narrow_PA, pa_narrow), 
        "Data/species/species_narrow.RDS")
saveRDS(list(virtual_intermediate, virtual_intermediate_PA, pa_intermediate), 
        "Data/species/species_intermediate.RDS")
saveRDS(list(virtual_broad, virtual_broad_PA, pa_broad), 
        "Data/species/species_broad.RDS")


plotResponse(virtual_narrow)
plotResponse(virtual_intermediate)
plotResponse(virtual_broad)


# =============================================================================
# PRESENCE SAMPLING - alle drei Arten
# → kontrolliertes Design für Algorithmenvergleich
# =============================================================================
# --- NARROW (specialist) ---

# Presence-only (für MaxEnt UND GLM presences - gleiche Punkte!)
po_narrow <- sampleOccurrences(
  virtual_narrow_PA,
  n = 100,
  type = "presence only",
  correct.by.suitability = TRUE
)

# GLM: 1000 pseudo-absences
bg_narrow_glm <- predicts::backgroundSample(
  virtual_narrow_PA$pa.raster,
  n = 1000
)

# BRT: fixe Presences aus po_narrow
pres_narrow <- po_narrow$sample.points[
  po_narrow$sample.points$Observed == 1, c("x", "y")]
pres_narrow$presence <- 1

# BRT: 10 Runs mit je verschiedenen pseudo-absences
brt_runs_narrow <- list()
for(i in 1:10){
  set.seed(i)
  abs_i <- as.data.frame(predicts::backgroundSample(
    virtual_narrow_PA$pa.raster, n = 100))
  abs_i$presence <- 0
  brt_runs_narrow[[i]] <- rbind(pres_narrow, abs_i)
}

# --- Speichern ---
saveRDS(list(
  po_narrow       = po_narrow,        # MaxEnt & GLM presences
  bg_narrow_glm   = bg_narrow_glm,    # GLM absences
  pres_narrow     = pres_narrow,       # BRT presences (fix)
  brt_runs_narrow = brt_runs_narrow    # BRT 10 Runs
), "Data/species/sampling_narrow.RDS")

# --- INTERMEDIATE ---
po_intermediate <- sampleOccurrences(
  virtual_intermediate_PA, n = 100,
  type = "presence only",
  correct.by.suitability = TRUE)

bg_intermediate_glm <- predicts::backgroundSample(
  virtual_intermediate_PA$pa.raster, n = 1000)

pres_intermediate <- po_intermediate$sample.points[
  po_intermediate$sample.points$Observed == 1, c("x", "y")]
pres_intermediate$presence <- 1

brt_runs_intermediate <- list()
for(i in 1:10){
  set.seed(i)
  abs_i <- as.data.frame(predicts::backgroundSample(
    virtual_intermediate_PA$pa.raster, n = 100))
  abs_i$presence <- 0
  brt_runs_intermediate[[i]] <- rbind(pres_intermediate, abs_i)
}

saveRDS(list(
  po_intermediate       = po_intermediate,
  bg_intermediate_glm   = bg_intermediate_glm,
  pres_intermediate     = pres_intermediate,
  brt_runs_intermediate = brt_runs_intermediate
), "Data/species/sampling_intermediate.RDS")


# --- BROAD ---
po_broad <- sampleOccurrences(
  virtual_broad_PA, n = 100,
  type = "presence only",
  correct.by.suitability = TRUE)

bg_broad_glm <- predicts::backgroundSample(
  virtual_broad_PA$pa.raster, n = 1000)

pres_broad <- po_broad$sample.points[
  po_broad$sample.points$Observed == 1, c("x", "y")]
pres_broad$presence <- 1

brt_runs_broad <- list()
for(i in 1:10){
  set.seed(i)
  abs_i <- as.data.frame(predicts::backgroundSample(
    virtual_broad_PA$pa.raster, n = 100))
  abs_i$presence <- 0
  brt_runs_broad[[i]] <- rbind(pres_broad, abs_i)
}

saveRDS(list(
  po_broad       = po_broad,
  bg_broad_glm   = bg_broad_glm,
  pres_broad     = pres_broad,
  brt_runs_broad = brt_runs_broad
), "Data/species/sampling_broad.RDS")
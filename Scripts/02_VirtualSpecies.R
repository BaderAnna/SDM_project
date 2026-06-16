# VIRTUELLE ARTEN ERSTELLEN ####
#-------------------------------------------#
# Gaussian response functions entlang PC1 und PC2
# σ-Werte: 0.2 (schmal), 0.5 (mittel), 0.8 (breit)

# --- Narrow niche (specialist, σ = 0.2) ---
virtual_narrow <- virtualspecies::generateSpFromFun(
  raster.stack = pc_raster,
  parameters = list(
    PC1 = list(fun = "dnorm", args = list(mean = 0, sd = 0.2)),
    PC2 = list(fun = "dnorm", args = list(mean = 0, sd = 0.2))
  ),
  rescale = TRUE,
  species.type = "multiplicative"
)

# --- Intermediate niche (σ = 0.5) ---
virtual_intermediate <- virtualspecies::generateSpFromFun(
  raster.stack = pc_raster,
  parameters = list(
    PC1 = list(fun = "dnorm", args = list(mean = 0, sd = 0.5)),
    PC2 = list(fun = "dnorm", args = list(mean = 0, sd = 0.5))
  ),
  rescale = TRUE,
  species.type = "multiplicative"
)

# --- Broad niche (generalist, σ = 0.8) ---
virtual_broad <- virtualspecies::generateSpFromFun(
  raster.stack = pc_raster,
  parameters = list(
    PC1 = list(fun = "dnorm", args = list(mean = 0, sd = 0.8)),
    PC2 = list(fun = "dnorm", args = list(mean = 0, sd = 0.8))
  ),
  rescale = TRUE,
  species.type = "multiplicative"
)

# Visualisieren
plot(virtual_narrow)
plot(virtual_intermediate)
plot(virtual_broad)

# Eignungskarte
plot(virtual_narrow$suitab.raster)

# PRESENCE POINTS SAMPELN ####
#-------------------------------------------#
set.seed(42)  # Reproduzierbarkeit

pa_narrow <- virtualspecies::sampleOccurrences(
  virtual_narrow$suitab.raster,
  n = 100,
  type = "presence only",
  correct.by.suitability = FALSE
)

pa_intermediate <- virtualspecies::sampleOccurrences(
  virtual_intermediate$suitab.raster,
  n = 100,
  type = "presence only",
  correct.by.suitability = FALSE
)

pa_broad <- virtualspecies::sampleOccurrences(
  virtual_broad$suitab.raster,
  n = 100,
  type = "presence only",
  correct.by.suitability = FALSE
)

# Visualisieren
plot(virtual_narrow$suitab.raster)
points(pa_narrow$sample.points[, c("x", "y")], pch = 19, cex = 0.5)
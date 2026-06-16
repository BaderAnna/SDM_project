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



set.seed(42)

# --- Narrow niche (specialist) ---
params_narrow <- formatFunctions(
  PC1 = c(fun = "dnorm", mean = 0, sd = 0.2),
  PC2 = c(fun = "dnorm", mean = 0, sd = 0.2)
)

virtual_narrow <- generateSpFromFun(
  raster.stack = pc_raster,
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
  raster.stack = pc_raster,
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
  raster.stack = pc_raster,
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
        "species_narrow.RDS")
saveRDS(list(virtual_intermediate, virtual_intermediate_PA, pa_intermediate), 
        "species_intermediate.RDS")
saveRDS(list(virtual_broad, virtual_broad_PA, pa_broad), 
        "species_broad.RDS")
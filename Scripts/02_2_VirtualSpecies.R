# =============================================================================
# VIRTUAL SPECIES CREATION & SAMPLING
# =============================================================================

library(terra)
library(virtualspecies)
library(predicts)

# Read data
pc_raster_masked <- terra::rast("Data/raster/pc_raster_masked.tif")

# Transform raster to metric CRS
pc_raster_masked <- terra::project(pc_raster_masked, "EPSG:3035")

# =============================================================================
# Helper function: sampling for one species
# =============================================================================
sample_species <- function(species_PA, n_pres = 1000,
                           n_glm_abs = 10000, n_brt_abs = 1000,
                           n_maxent_bg = 10000, n_runs = 10) {
  
  # Presence-only (MaxEnt & GLM)
  po <- sampleOccurrences(
    species_PA, n = n_pres,
    type = "presence only",
    correct.by.suitability = TRUE
  )
  
  # GLM: 10000 random pseudo-absences
  bg_glm <- as.data.frame(predicts::backgroundSample(
    species_PA$pa.raster, n = n_glm_abs
  ))
  
  # MaxEnt: 10000 background points
  bg_maxent <- as.data.frame(predicts::backgroundSample(
    species_PA$pa.raster, n = n_maxent_bg
  ))
  
  # Extract presence points
  pres <- po$sample.points[po$sample.points$Observed == 1, c("x", "y")]
  pres$presence <- 1
  
  # BRT & RF: 10 runs with n_brt_abs random pseudo-absences each
  runs <- list()
  for (i in 1:n_runs) {
    set.seed(i)
    abs_i <- as.data.frame(predicts::backgroundSample(
      species_PA$pa.raster, n = n_brt_abs
    ))
    abs_i$presence <- 0
    runs[[i]] <- rbind(pres, abs_i)
  }
  
  return(list(
    po        = po,         # presence-only (MaxEnt & GLM presences)
    bg_glm    = bg_glm,     # GLM pseudo-absences
    bg_maxent = bg_maxent,  # MaxEnt background points
    pres      = pres,       # BRT/RF presences (fixed)
    runs      = runs        # BRT/RF 10 runs
  ))
}

# =============================================================================
# Virtual species & sampling
# =============================================================================
set.seed(42)

species_list <- list(
  narrow   = list(sd = 0.2),
  low_mid  = list(sd = 0.4),
  high_mid = list(sd = 0.6),
  broad    = list(sd = 0.8)
)

dir.create("Data/species", recursive = TRUE, showWarnings = FALSE)

for (sp in names(species_list)) {
  
  sd_val <- species_list[[sp]]$sd
  
  # --- Define response curves ---
  params <- formatFunctions(
    PC1 = c(fun = "dnorm", mean = 0, sd = sd_val),
    PC2 = c(fun = "dnorm", mean = 0, sd = sd_val)
  )
  
  # --- Create virtual species ---
  virtual_sp <- generateSpFromFun(
    raster.stack = pc_raster_masked,
    parameters   = params,
    species.type = "multiplicative",
    plot         = FALSE
  )
  virtual_sp$species.name <- sp
  
  # --- Convert to PA ---
  virtual_sp_PA <- convertToPA(
    virtual_sp,
    beta  = 0.5,
    alpha = -0.05,
    plot  = FALSE
  )
  
  # --- Sampling ---
  sampling <- sample_species(virtual_sp_PA)
  
  # --- Save ---
  saveRDS(
    list(virtual    = virtual_sp,
         virtual_PA = virtual_sp_PA),
    paste0("Data/species/species_", sp, ".RDS")
  )
  
  saveRDS(
    sampling,
    paste0("Data/species/sampling_", sp, ".RDS")
  )
  
  message("✓ ", sp, " (σ = ", sd_val, ") saved | Presences = ",
          sum(sampling$po$sample.points$Observed == 1))
}

summary(virtual_sp$suitab.raster)  # check range of suitability values


# =============================================================================
# Visualisation: suitability maps of the 4 virtual species
# =============================================================================
library(ggplot2)
library(tidyterra)
library(patchwork)
library(terra)

species_list <- list(
  narrow   = list(sd = 0.2),
  low_mid  = list(sd = 0.4),
  high_mid = list(sd = 0.6),
  broad    = list(sd = 0.8)
)

# Labels for the plots
species_labels <- c(
  narrow   = "Narrow\n(σ = 0.2)",
  low_mid  = "Low-Mid\n(σ = 0.4)",
  high_mid = "High-Mid\n(σ = 0.6)",
  broad    = "Broad\n(σ = 0.8)"
)

# Plot for each species
plots <- lapply(names(species_list), function(sp) {
  
  obj     <- readRDS(paste0("Data/species/species_", sp, ".RDS"))
  suit    <- obj$virtual$suitab.raster
  
  ggplot() +
    tidyterra::geom_spatraster(data = suit) +
    scale_fill_viridis_c(
      option   = "viridis",
      na.value = "white",
      name     = "Suitability",
      limits   = c(0, 1)        # uniform scale across all 4!
    ) +
    labs(title = species_labels[sp]) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 10, hjust = 0.5,
                                     face = "bold"),
      legend.position = "none",
      plot.background = element_rect(fill      = "white",
                                     color     = "grey80",
                                     linewidth = 0.5)
    )
})

# Extract shared legend
legend_plot <- ggplot() +
  tidyterra::geom_spatraster(
    data = readRDS("Data/species/species_narrow.RDS")$virtual$suitab.raster
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    name   = "Suitability",
    limits = c(0, 1)
  ) +
  theme_void() +
  theme(
    legend.position    = "bottom",
    legend.title       = element_text(size = 9,  face = "bold"),
    legend.text        = element_text(size = 8),
    legend.key.height  = unit(1.5, "cm"),
    legend.key.width   = unit(0.4, "cm")
  )

# Extract legend
get_legend <- function(p) {
  gt   <- ggplot_gtable(ggplot_build(p))
  leg  <- which(sapply(gt$grobs, function(x) x$name) == "guide-box")
  gt$grobs[[leg]]
}

legend <- get_legend(legend_plot)

# Combine plots
wrap_plots(plots, nrow = 2) +
  plot_annotation(
    title    = "Virtual Species - Habitat Suitability",
    theme    = theme(
      plot.title    = element_text(size = 14, hjust = 0.5,
                                   face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5,
                                   color = "grey40")
    )
  )
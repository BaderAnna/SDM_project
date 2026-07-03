# 1 - install and load packages  ####
#-----------------------------------#

# Define packages
list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "geodata")

# Check and install missing packages
for (pkg in list.of.packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Install the package if not already installed
    install.packages(pkg, dependencies = TRUE)
  }
  # Load the package
  library(pkg, character.only = TRUE)
}

# Load the packages
lapply(list.of.packages, library, character.only = TRUE)

# 2 - prepare some bioclimatic variables ####
#-------------------------------------------#
bioclim= geodata::worldclim_country(country="Germany", path="Data/raster", var="bio")

r=terra::rast("Data/raster/climate/wc2.1_country/DEU_wc2.1_30s_bio.tif")

# have a look at your rasters  
terra::plot(r)  

# 3 - PCA ####
#-------------------------------------------#

# Select variables (BIO1, BIO4, BIO12, BIO15)
bio_subset <- r[[c(1, 4, 12, 15)]]

# have a look at your rasters  
terra::plot(bio_subset) 

# matrix 
vals <- values(bio_subset, na.rm=TRUE)

# PCA 
pca <- prcomp(vals, scale. = TRUE, center = TRUE)

# Varianzaufklärung anschauen
summary(pca)

# Wie viel Varianz erklären PC1 und PC2?
pca$sdev^2 / sum(pca$sdev^2) * 100

# Wertebereich der PC-Scores anschauen
summary(pca$x[,1])  # PC1
summary(pca$x[,2])  # PC2
sd(pca$x[,1])
sd(pca$x[,2])

# Loadings anschauen
pca$rotation

# PCA RASTER ERSTELLEN ####
#-------------------------------------------#
# PCA auf die Rasterdaten anwenden
# Erst alle Werte extrahieren (inkl. NA-Positionen für spätere Rückprojektion)
vals_all <- values(bio_subset)

# PCA-Scores für alle Pixel berechnen
scores_all <- predict(pca, newdata = vals_all)

# Raster mit PC1 und PC2 erstellen
pc1_rast <- bio_subset[[1]]  # als Template
pc2_rast <- bio_subset[[1]]

values(pc1_rast) <- scores_all[, 1]
values(pc2_rast) <- scores_all[, 2]

# Stapeln zu einem 2-Layer Raster
pc_raster <- c(pc1_rast, pc2_rast)
names(pc_raster) <- c("PC1", "PC2")

# Visualisieren
terra::plot(pc_raster)

# 4 - Germany ####
#-------------------------------------------#
library(rnaturalearth)
library(rnaturalearthdata)

# Deutschland-Grenze laden
germany_sf <- rnaturalearth::ne_countries(
  country = "Germany", 
  scale = "medium", 
  returnclass = "sf"
)

# Raster maskieren
pc_raster_masked <- terra::mask(pc_raster, terra::vect(germany_sf))

terra::plot(pc_raster_masked)

# PC Raster speichern
terra::writeRaster(
  pc_raster_masked,
  filename  = "Data/raster/pc_raster_masked.tif",
  overwrite = TRUE
)

# PCA speichern
saveRDS(pca, "Data/pca.RDS")

# =============================================================================
# Abbildung
# =============================================================================

library(ggplot2)
library(tidyterra)
library(patchwork)
library(rnaturalearth)
library(terra)

# =============================================================================
# 1 - Raster laden + auf Deutschland zuschneiden
# =============================================================================

r <- terra::rast("Data/raster/climate/wc2.1_country/DEU_wc2.1_30s_bio.tif")

# Deutschland-Grenze laden
germany_sf <- rnaturalearth::ne_countries(
  country     = "Germany",
  scale       = "medium",
  returnclass = "sf"
)

# Raster maskieren
r_masked <- terra::mask(r, terra::vect(germany_sf))

# =============================================================================
# 2 - Plot der 4 Variablen
# =============================================================================

vars <- list(
  list(name = "wc2.1_30s_bio_1",  label = "BIO1\nMean Annual Temp."),
  list(name = "wc2.1_30s_bio_4",  label = "BIO4\nTemp. Seasonality"),
  list(name = "wc2.1_30s_bio_12", label = "BIO12\nAnnual Precip."),
  list(name = "wc2.1_30s_bio_15", label = "BIO15\nPrecip.\nSeasonality")
)

plots <- lapply(vars, function(v) {
  ggplot() +
    tidyterra::geom_spatraster(data = r_masked[[v$name]]) +
    scale_fill_viridis_c(
      option   = "ciidis",  # "magma", "plasma", "inferno", "cividis"
      na.value = "white",
      name     = ""
    ) +
    labs(title = v$label) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 8, hjust = 0.5,
                                     face = "bold"),
      legend.position = "none",
      plot.background = element_rect(fill      = "white",
                                     color     = "grey80",
                                     linewidth = 0.5)
    )
})

# Zusammenfügen
wrap_plots(plots, nrow = 2) +
  plot_annotation(
    title = "Environmental Variables",
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0.5,
                                face = "bold")
    )
  )

# =============================================================================
# Abbildung PC1 und PC2
# =============================================================================

# PC Raster laden
pc_raster_masked <- terra::rast("Data/raster/pc_raster_masked.tif")

# PCA laden für Loadings
pca <- readRDS("Data/pca.RDS")

# Loadings anzeigen
print(round(pca$rotation[, 1:2], 3))

# Varianzaufklärung
var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
message("PC1 erklärt: ", round(var_explained[1], 1), "%")
message("PC2 erklärt: ", round(var_explained[2], 1), "%")

# =============================================================================
# Plot PC1 und PC2
# =============================================================================

pc_vars <- list(
  list(
    name   = "PC1",
    label  = paste0("PC1"),
    option = "cividis"
  ),
  list(
    name   = "PC2",
    label  = paste0("PC2"),
    option = "cividis"
  )
)

plots_pc <- lapply(pc_vars, function(v) {
  ggplot() +
    tidyterra::geom_spatraster(data = pc_raster_masked[[v$name]]) +
    scale_fill_viridis_c(
      option   = v$option,
      na.value = "white",
      name     = ""
    ) +
    labs(title = v$label) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 9, hjust = 0.5,
                                     face = "bold"),
      legend.position = "none",
      plot.background = element_rect(fill      = "white",
                                     color     = "grey80",
                                     linewidth = 0.5)
    )
})


# Zusammenfügen
wrap_plots(plots_pc, nrow = 2) +
  plot_annotation(
    title    = "Principal Components",
    theme    = theme(
      plot.title    = element_text(size = 14, hjust = 0.5,
                                   face = "bold"),
      plot.subtitle = element_text(size =  8, hjust = 0.5,
                                   color = "grey40",
                                   family = "mono")
    )
  )

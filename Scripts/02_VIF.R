library(terra)
library(sf)
library(usdm)
library(CAST)
library(rnaturalearth)
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

# =============================================================================
# WorldClim Raster laden, VIF-Reduktion, maskieren, speichern
# =============================================================================

# 1. Raster einlesen
r <- terra::rast("Data/raster/climate/wc2.1_country/DEU_wc2.1_30s_bio.tif")

# 2. Germany-Grenze laden (VOR allem anderen)
germany_sf <- rnaturalearth::ne_countries(country = "Germany", 
                                          scale = "medium", 
                                          returnclass = "sf")

# 3. ERST maskieren/croppen, DANN VIF berechnen
r_masked <- terra::mask(r, terra::vect(germany_sf))
r_masked <- terra::crop(r_masked, terra::vect(germany_sf))

# 4. Werte aus der KORREKT maskierten Version extrahieren
vals_19 <- as.data.frame(terra::values(r_masked, na.rm = TRUE))

# 5. VIF auf Basis der korrekten Studienfläche
set.seed(42)
vif_result <- usdm::vifstep(vals_19, th = 10)
selected_vars <- as.character(vif_result@results$Variables)
message("Behaltene Variablen: ", paste(selected_vars, collapse = ", "))

# 6. Reduziertes Raster erstellen (schon maskiert)
r_reduced <- r_masked[[selected_vars]]
names(r_reduced) <- selected_vars

# 7. Projektion (Maskierung ist ja schon erfolgt, Reihenfolge hier weniger kritisch)
r_reduced_3035 <- terra::project(r_reduced, "EPSG:3035")
names(r_reduced_3035) <- selected_vars

germany_sf_3035 <- sf::st_transform(germany_sf, crs = 3035)

# 8. Finale Sicherheits-Maskierung (falls durch project() Randpixel entstehen)
env_raster_masked <- terra::mask(r_reduced_3035, terra::vect(germany_sf_3035))
names(env_raster_masked) <- selected_vars

# 9. Speichern
dir.create("Data/raster", recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(env_raster_masked, "Data/raster/env_raster_masked.tif", overwrite = TRUE)
saveRDS(selected_vars, "Data/raster/selected_vars.RDS")

message("Raster-Namen: ", paste(names(env_raster_masked), collapse = ", "))
terra::plot(env_raster_masked)


# =============================================================================
# Plot der 7 ausgewählten Bioklimavariablen (nach VIF-Auswahl)
# =============================================================================

# =============================================================================
# Plot der 7 ausgewählten Bioklimavariablen (nach VIF-Auswahl) - mit Klarnamen
# =============================================================================

# Lookup-Tabelle: technischer Name -> beschreibender Titel
bio_labels <- c(
  wc2.1_30s_bio_1  = "Mean Temp.",
  wc2.1_30s_bio_2  = "Mean Diurnal Range",
  wc2.1_30s_bio_3  = "Isothermality",
  wc2.1_30s_bio_4  = "Temp. Seasonality",
  wc2.1_30s_bio_5  = "Max Temp. Warmest Month",
  wc2.1_30s_bio_6  = "Min Temp. Coldest Month",
  wc2.1_30s_bio_7  = "Temp. Annual Range",
  wc2.1_30s_bio_8  = "Mean Temp. Wettest Quarter",
  wc2.1_30s_bio_9  = "Mean Temp. Driest Quarter",
  wc2.1_30s_bio_10 = "Mean Temp. Warmest Quarter",
  wc2.1_30s_bio_11 = "Mean Temp. Coldest Quarter",
  wc2.1_30s_bio_12 = "Annual Precip.",
  wc2.1_30s_bio_13 = "Precip. Wettest Month",
  wc2.1_30s_bio_14 = "Precip. Driest Month",
  wc2.1_30s_bio_15 = "Precip. Seasonality",
  wc2.1_30s_bio_16 = "Precip. Wettest Quarter",
  wc2.1_30s_bio_17 = "Precip. Driest Quarter",
  wc2.1_30s_bio_18 = "Precip. Warmest Quarter",
  wc2.1_30s_bio_19 = "Precip. Coldest Quarter"
)

# Variablen aus der VIF-Auswahl laden
selected_vars <- readRDS("Data/raster/selected_vars.RDS")

env_vars <- lapply(selected_vars, function(v) {
  list(
    name   = v,
    label  = bio_labels[[v]],   # Klarname statt technischem Namen
    option = "cividis"
  )
})

plots_env <- lapply(env_vars, function(v) {
  ggplot() +
    tidyterra::geom_spatraster(data = env_raster_masked[[v$name]]) +
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
wrap_plots(plots_env, ncol = 4) +
  plot_annotation(
    title    = "Selected Bioclimatic Variables",
    theme    = theme(
      plot.title    = element_text(size = 14, hjust = 0.5,
                                   face = "bold"),
      plot.subtitle = element_text(size =  8, hjust = 0.5,
                                   color = "grey40",
                                   family = "mono")
    )
  )

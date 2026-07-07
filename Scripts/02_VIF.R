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


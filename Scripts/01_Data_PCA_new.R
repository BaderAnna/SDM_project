# =============================================================================
# 1 - INSTALL AND LOAD PACKAGES
# =============================================================================
list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "geodata",
                      "rnaturalearth", "rnaturalearthdata")

for (pkg in list.of.packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
invisible(lapply(list.of.packages, library, character.only = TRUE))

# =============================================================================
# 2 - PREPARE BIOCLIMATIC VARIABLES
# =============================================================================
bioclim <- geodata::worldclim_country(country = "Germany", path = "Data/raster", var = "bio")
r <- terra::rast("Data/raster/climate/wc2.1_country/DEU_wc2.1_30s_bio.tif")

terra::plot(r)

# Select variables (BIO1, BIO4, BIO12, BIO15)
bio_subset <- r[[c(1, 4, 12, 15)]]
names(bio_subset) <- c("BIO1", "BIO4", "BIO12", "BIO15")
terra::plot(bio_subset)

# =============================================================================
# 3 - GERMANY BOUNDARY (rnaturalearth) - LOAD BEFORE MASKING
# =============================================================================
germany_sf <- rnaturalearth::ne_countries(
  country      = "Germany",
  scale        = "medium",
  returnclass  = "sf"
)

# =============================================================================
# 4 - MASK & CROP FIRST, THEN COMPUTE PCA
#     (ensures PCA is based only on the exact study area, not the
#      raw WorldClim/GADM bounding box, avoiding boundary artifacts)
# =============================================================================
bio_subset_masked <- terra::mask(bio_subset, terra::vect(germany_sf))
bio_subset_masked <- terra::crop(bio_subset_masked, terra::vect(germany_sf))

terra::plot(bio_subset_masked)

# =============================================================================
# 5 - PCA (on masked data, with cell indices retained for traceability)
# =============================================================================
# Extract values + cell IDs together, so PCA rows can always be traced
# back to their spatial location
all_vals <- terra::values(bio_subset_masked)
colnames(all_vals) <- c("BIO1", "BIO4", "BIO12", "BIO15")

complete_rows <- stats::complete.cases(all_vals)
cell_ids      <- which(complete_rows)

vals_with_cells <- all_vals[complete_rows, ]

# PCA
pca <- prcomp(vals_with_cells, scale. = TRUE, center = TRUE)

# Check: muss identisch sein
length(cell_ids)
nrow(pca$x)

# --- Sanity check: rows must match ---
stopifnot(nrow(vals_with_cells) == nrow(pca$x))

# --- Variance explained ---
summary(pca)
pca$sdev^2 / sum(pca$sdev^2) * 100

# --- Range of PC scores ---
summary(pca$x[, 1])  # PC1
summary(pca$x[, 2])  # PC2
sd(pca$x[, 1])
sd(pca$x[, 2])

# --- Loadings ---
pca$rotation

# =============================================================================
# 6 - OUTLIER CHECK (spatially locate extreme PC scores)
# =============================================================================
# Boxplot / histogram overview
boxplot(pca$x[, 1], main = "PC1 Scores")
boxplot(pca$x[, 2], main = "PC2 Scores")
hist(pca$x[, 2], breaks = 100, main = "Distribution PC2")

# =============================================================================
# 7 - PCA RASTER ERSTELLEN (apply PCA to all pixels, incl. those with NA
#     in at least one layer get NA in the PCA raster automatically)
# =============================================================================
vals_all <- terra::values(bio_subset_masked)
colnames(vals_all) <- c("BIO1", "BIO4", "BIO12", "BIO15")

scores_all <- predict(pca, newdata = vals_all)

pc1_rast <- bio_subset_masked[[1]]  # template
pc2_rast <- bio_subset_masked[[1]]

values(pc1_rast) <- scores_all[, 1]
values(pc2_rast) <- scores_all[, 2]

pc_raster <- c(pc1_rast, pc2_rast)
names(pc_raster) <- c("PC1", "PC2")

terra::plot(pc_raster)

# =============================================================================
# 8 - FINAL MASK (redundant safety step, ensures no stray pixels outside
#     Germany remain, e.g. from crop() bounding box corners)
# =============================================================================
pc_raster_masked <- terra::mask(pc_raster, terra::vect(germany_sf))
terra::plot(pc_raster_masked)

# =============================================================================
# 9 - SAVE OUTPUTS
# =============================================================================
dir.create("Data/raster", recursive = TRUE, showWarnings = FALSE)

terra::writeRaster(
  pc_raster_masked,
  filename  = "Data/raster/pc_raster_masked.tif",
  overwrite = TRUE
)

saveRDS(pca, "Data/pca.RDS")
saveRDS(vals_with_cells, "Data/pca_vals_with_cells.RDS")  
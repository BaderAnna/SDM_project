# ----------------------------------------------------------------
# Packages
# ----------------------------------------------------------------
library(terra)
library(sf)
library(sp)
library(gstat)

# ----------------------------------------------------------------
# Load raster and reproject to a metric CRS
# ----------------------------------------------------------------
r      <- terra::rast("Data/raster/climate/wc2.1_country/DEU_wc2.1_30s_bio.tif")
r_3035 <- terra::project(r, "EPSG:3035")

# Select the 4 bioclim layers directly from the reprojected raster (no separate 'pc' object needed)
names(r_3035)
bio_subset <- r_3035[[c("wc2.1_30s_bio_1", "wc2.1_30s_bio_4",
                        "wc2.1_30s_bio_12", "wc2.1_30s_bio_15")]]

# ----------------------------------------------------------------
# Variogram, using a sample of 5000 points
# ----------------------------------------------------------------
df <- as.data.frame(r_3035, xy = TRUE, na.rm = TRUE)
coordinates(df)  <- ~ x + y
proj4string(df)  <- CRS(terra::crs(r_3035))

set.seed(42)
df_sample <- df[sample(nrow(df), 5000), ]

v <- gstat::variogram(wc2.1_30s_bio_1 ~ 1, data = df_sample)
plot(v)

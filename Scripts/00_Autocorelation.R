list.of.packages <- c("terra", "sf", "predicts", "virtualspecies", "geodata")

r=terra::rast("Data/raster/climate/wc2.1_country/DEU_wc2.1_30s_bio.tif")

# Raster in metrisches KBS transformieren
r_3035 <- terra::project(r, "EPSG:3035")

bio_subset <- pc[[c(1, 4, 12, 15)]]

# Variogram
install.packages("ctmm")
library(ctmm)
#variogram(r_3035$wc2.1_30s_bio_1)

library(sf)
df <- as.data.frame(r_3035, xy = TRUE, na.rm = TRUE)
sf_df <- st_as_sf(df, coords = c("x", "y"), crs = crs(r_3035))

v <- variogram(values ~ 1, data = sf_df)

library(terra)
library(sp)
library(gstat)

df <- as.data.frame(r_3035, xy = TRUE, na.rm = TRUE)

coordinates(df) <- ~x+y
proj4string(df) <- CRS(terra::crs(r_3035))

set.seed(42)
df_sample <- df[sample(nrow(df), 5000), ]

v <- variogram(wc2.1_30s_bio_1 ~ 1, data = df_sample)
plot(v)


# Coreogramd 
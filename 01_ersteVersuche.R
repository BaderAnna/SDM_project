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

# make sure the path to the folder raster exists on your PC:
bioclim= geodata::worldclim_country(country="Germany", path="raster", var="bio")

r=terra::rast("raster/climate/wc2.1_country/DEU_wc2.1_30s_bio.tif")
# have a look at your rasters  
terra::plot(r)  

# 3 - PCA ####
#-------------------------------------------#

# BIO1 (annual mean temperature) und BIO12 (annual precipitation) auswählen
bio1  <- r[[1]]   # BIO1
bio12 <- r[[12]]  # BIO12

# Wertebereich und Standardabweichung der Variablen anschauen
summary(values(bio1),  na.rm=TRUE)
summary(values(bio12), na.rm=TRUE)

sd(values(bio1),  na.rm=TRUE)
sd(values(bio12), na.rm=TRUE)

## PCA
# Variablen auswählen (z.B. BIO1, BIO4, BIO12, BIO15)
bio_subset <- r[[c(1, 4, 12, 15)]]

# In Matrix umwandeln für PCA
vals <- values(bio_subset, na.rm=TRUE)

# PCA durchführen (skaliert + zentriert)
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

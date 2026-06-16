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

getwd()

# 2 - prepare some bioclimatic variables ####
#-------------------------------------------#

# make sure the path to the folder raster exists on your PC:
bioclim= geodata::worldclim_country(country="Germany", path="Data/raster", var="bio")

r=terra::rast("Data/raster/climate/wc2.1_country/DEU_wc2.1_30s_bio.tif")
# have a look at your rasters  
terra::plot(r)  

# 3 - PCA ####
#-------------------------------------------#

# Variablen auswählen (BIO1, BIO4, BIO12, BIO15)
bio_subset <- r[[c(1, 4, 12, 15)]]

# have a look at your rasters  
terra::plot(bio_subset) 

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

#### GLM ---------------------------

# Erst Daten zusammenführen
pres_narrow_glm <- po_narrow$sample.points[
  po_narrow$sample.points$Observed == 1, c("x", "y")]
pres_narrow_glm$presence <- 1

abs_narrow_glm <- as.data.frame(bg_narrow_glm)
abs_narrow_glm$presence <- 0

glm_data_narrow <- rbind(pres_narrow_glm, abs_narrow_glm)

# Umweltvariablen extrahieren
glm_data_narrow <- cbind(
  glm_data_narrow,
  as.data.frame(terra::extract(pc_raster, 
                               glm_data_narrow[, c("x", "y")]))
)

# Gewichtung definieren
weights_narrow <- c(rep(1, 100),    # Presences
                    rep(0.1, 1000)) # Absences

# GLM fitten MIT Gewichtung
glm_narrow <- glm(
  presence ~ PC1 + PC2 + I(PC1^2) + I(PC2^2),
  data    = glm_data_narrow,
  family  = binomial,
  weights = weights_narrow
)

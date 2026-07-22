# =============================================================================
# TSS und AUC
# =============================================================================

# Install mecofun from Gitlab using the devtools package:
#library(devtools)
#devtools::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")

#https://damariszurell.github.io/EEC-MGC/b4_SDM_eval.html # quasi ein Tutorial

# citation
paste("Zurell, D. (2024).",
      "mecofun: useful functions for macroecology and species distribution modelling",
      "version 0.7.1. University of Potsdam, Potsdam.",
      "https://gitup.uni-potsdam.de/macroecology/mecofun")



# =============================================================================
# TSS-Auswertung: bereits trainierte GLMs (pro Nischenbreite) auf
# unabhaengigen knndm-Testdatensaetzen
# Threshold: Prevalence-Ansatz (ObsPrev) via mecofun::evalSDM()
# =============================================================================

library(mecofun)
library(PresenceAbsence)
library(terra)

# -----------------------------------------------------------------------------
# 1 - Umweltraster laden (Prädiktoren, mit denen die Modelle gefittet wurden)
# -----------------------------------------------------------------------------
# ANPASSEN: Pfad zu deinem Praediktoren-Stack. Die Layer-Namen MUESSEN exakt
# den Namen entsprechen, die im Modell verwendet wurden (z.B. "wc2.1_30s_bio_1").
env <- terra::rast("Data/raster/env_raster_masked.tif")   
print(names(env))   # Kontrolle: muss z.B. wc2.1_30s_bio_1, _3, _4, _8, _9, _15, _18 enthalten


# -----------------------------------------------------------------------------
# 2 - Modelle und Testdaten laden
# -----------------------------------------------------------------------------
models <- list(
  narrow   = readRDS("Data/models/knndm/glm_narrow_train.RDS"),
  low_mid  = readRDS("Data/models/knndm/glm_low_mid_train.RDS"),
  high_mid = readRDS("Data/models/knndm/glm_high_mid_train.RDS"),
  broad    = readRDS("Data/models/knndm/glm_broad_train.RDS")
)

test_data <- list(
  narrow   = readRDS("Data/species/split_knndm/test_narrow.RDS"),
  low_mid  = readRDS("Data/species/split_knndm/test_low_mid.RDS"),
  high_mid = readRDS("Data/species/split_knndm/test_high_mid.RDS"),
  broad    = readRDS("Data/species/split_knndm/test_broad.RDS")
)


# -----------------------------------------------------------------------------
# 3 - Funktion: Umweltwerte extrahieren, vorhersagen, TSS/Threshold berechnen
# -----------------------------------------------------------------------------

eval_one_niche <- function(model_obj, test_df, env, thresh.method = "ObsPrev") {
  
  # Nur echte Test-Zeilen behalten (defensiv, falls Objekt mehr enthaelt)
  test_df <- test_df[test_df$split == "test", ]
  
  # Umweltwerte an den Testpunkt-Koordinaten extrahieren
  env_vals <- terra::extract(env, test_df[, c("x", "y")])
  env_vals <- env_vals[, -1, drop = FALSE]   # ID-Spalte von terra::extract entfernen
  
  # Kontrolle: passen die Spaltennamen zu den Modell-Praediktoren?
  # (bei Nichtuebereinstimmung wirft predict() weiter unten einen Fehler)
  newdat <- cbind(test_df, env_vals)
  
  # Der eigentliche glm-Fit liegt unter $model, nicht im Objekt selbst
  glm_fit <- model_obj$model
  
  pred <- predict(glm_fit, newdata = newdat, type = "response")
  
  # Zeilen mit NA-Vorhersage (z.B. Punkte ausserhalb des Rasterbereichs) ausschliessen
  ok <- !is.na(pred)
  if (any(!ok)) {
    warning(sprintf("%d von %d Testpunkten konnten nicht vorhergesagt werden (NA) - werden ausgeschlossen.",
                    sum(!ok), length(ok)))
  }
  
  eval_res <- evalSDM(observation = newdat$presence[ok],
                      predictions = pred[ok],
                      thresh.method = thresh.method)
  
  eval_res$true_prevalence <- mean(newdat$presence[ok])
  eval_res$n_test          <- sum(ok)
  
  eval_res
}


# -----------------------------------------------------------------------------
# 4 - Fuer alle Nischenbreiten anwenden
# -----------------------------------------------------------------------------

results_list <- Map(eval_one_niche, models, test_data,
                    MoreArgs = list(env = env, thresh.method = "ObsPrev")) 

results_df <- do.call(rbind, results_list)
results_df$niche_breadth <- names(models)
rownames(results_df) <- NULL

print(results_df)


# -----------------------------------------------------------------------------
# 5 - Bootstrap-Konfidenzintervalle
# -----------------------------------------------------------------------------

bootstrap_tss <- function(obs, pred, thresh.method = "ObsPrev", n_boot = 1000) {
  n <- length(obs)
  tss_boot <- numeric(n_boot)
  
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    obs_b  <- obs[idx]
    pred_b <- pred[idx]
    
    # Falls ein Bootstrap-Sample zufaellig keine Presence oder keine Absence
    # enthaelt, evalSDM()/optimal.thresholds() schlaegt fehl - dann skippen
    if (length(unique(obs_b)) < 2) next
    
    res <- tryCatch(
      evalSDM(observation = obs_b, predictions = pred_b, thresh.method = thresh.method),
      error = function(e) NULL
    )
    tss_boot[i] <- if (!is.null(res)) res$TSS else NA
  }
  
  quantile(tss_boot, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
}

# Anwendung pro Nischenbreite (obs/pred aus deiner eval_one_niche-Logik heraus extrahiert):
tss_ci <- lapply(names(models), function(nm) {
  test_df <- test_data[[nm]]
  test_df <- test_df[test_df$split == "test", ]
  env_vals <- terra::extract(env, test_df[, c("x", "y")])[, -1, drop = FALSE]
  newdat <- cbind(test_df, env_vals)
  pred <- predict(models[[nm]]$model, newdata = newdat, type = "response")
  ok <- !is.na(pred)
  
  bootstrap_tss(newdat$presence[ok], pred[ok])
})
names(tss_ci) <- names(models)
print(tss_ci)



# -----------------------------------------------------------------------------
# 6 - Speichern und Plot
# -----------------------------------------------------------------------------

write.csv(results_df, "results/glm_auc_tss.csv", row.names = FALSE)

plot(seq_along(results_df$niche_breadth), results_df$TSS,
     xaxt = "n", type = "b", pch = 19,
     xlab = "Nischenbreite", ylab = "TSS",
     main = "TSS vs. Nischenbreite (knndm-Testdaten)")
axis(1, at = seq_along(results_df$niche_breadth), labels = results_df$niche_breadth)

plot(seq_along(results_df$niche_breadth), results_df$AUC,
     xaxt = "n", type = "b", pch = 19,
     xlab = "Nischenbreite", ylab = "AUC",
     main = "AUC vs. Nischenbreite (knndm-Testdaten)")
axis(1, at = seq_along(results_df$niche_breadth), labels = results_df$niche_breadth)













# -----------------------------------------------------------------------------
# Presence-Zahlen pro Nischenbreite: Training und Test im Ueberblick
# -----------------------------------------------------------------------------

presence_overview <- data.frame(
  niche_breadth = names(models),
  n_presence_train = sapply(models, function(m) sum(m$model$model[, 1] == 1)),
  n_absence_train  = sapply(models, function(m) sum(m$model$model[, 1] == 0)),
  n_presence_test  = sapply(test_data, function(d) sum(d$presence[d$split == "test"] == 1)),
  n_absence_test   = sapply(test_data, function(d) sum(d$presence[d$split == "test"] == 0))
)

presence_overview$prevalence_train <- presence_overview$n_presence_train /
  (presence_overview$n_presence_train + presence_overview$n_absence_train)

presence_overview$prevalence_test <- presence_overview$n_presence_test /
  (presence_overview$n_presence_test + presence_overview$n_absence_test)

print(presence_overview)






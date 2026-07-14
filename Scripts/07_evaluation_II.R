# =============================================================================
# 1 - TSS
# =============================================================================

# Install mecofun from Gitlab using the devtools package:
library(devtools)
devtools::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")

https://damariszurell.github.io/EEC-MGC/b4_SDM_eval.html # quasi ein Tutorial

# citation
paste("Zurell, D. (2024).",
      "mecofun: useful functions for macroecology and species distribution modelling",
      "version 0.7.1. University of Potsdam, Potsdam.",
      "https://gitup.uni-potsdam.de/macroecology/mecofun")


library(mecofun)


# hier alle VS reinladen inkl. predictions




# --- 2c: TSS + Threshold via mecofun::evalSDM() ---
eval_res <- evalSDM(observation = obs, predictions = pred,
                    thresh.method = thresh.method)

# wahre Praevalenz mit ausgeben (zur Kontrolle/Interpretation)
eval_res$true_prevalence <- mean(obs)
eval_res$niche_breadth_sd <- sd_val

eval_res
}


# -----------------------------------------------------------------------------
# 3 - Loop ueber alle Nischenbreiten
# -----------------------------------------------------------------------------

results_list <- lapply(niche_breadths, function(sd_val) {
  run_one_species(sd_val = sd_val, env = env, n_points = n_points,
                  thresh.method = "ObsPrev") # hier wird die Prevalance approach verwendet

results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

print(results_df)


# -----------------------------------------------------------------------------
# 4 - Ergebnisse speichern und visualisieren
# -----------------------------------------------------------------------------

write.csv(results_df, "tss_nischenbreiten_vergleich.csv", row.names = FALSE)

# Einfacher Plot: TSS in Abhaengigkeit von Nischenbreite / Praevalenz
par(mfrow = c(1, 2))

plot(results_df$niche_breadth_sd, results_df$TSS,
     type = "b", pch = 19,
     xlab = "Nischenbreite (sd)", ylab = "TSS",
     main = "TSS vs. Nischenbreite")

plot(results_df$true_prevalence, results_df$TSS,
     type = "b", pch = 19,
     xlab = "Wahre Praevalenz", ylab = "TSS",
     main = "TSS vs. Praevalenz")





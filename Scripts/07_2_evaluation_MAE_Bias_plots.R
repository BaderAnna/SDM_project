# =============================================================================
# Plots: MAE + Bias (vs. wahre Suitability)
# =============================================================================
#
# Liest die Ergebnisse aus calc_mae_bias_true_suitability.R ein
# (results/mae_bias_true_suitability.csv) und erstellt die Plots.
# =============================================================================

library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# 1 - Daten einlesen
# -----------------------------------------------------------------------------

niche_names  <- c("narrow", "low_mid", "high_mid", "broad")
niche_labels <- c("Narrow", "Low", "High", "Broad")
model_order  <- c("GLM", "Maxent", "BRT", "RF")

metrics_combined <- read.csv("results/mae_bias_true_suitability.csv") %>%
  mutate(
    niche_breadth = factor(niche_breadth, levels = niche_names, labels = niche_labels),
    model_type    = factor(model_type, levels = model_order)
  )

# -----------------------------------------------------------------------------
# 2 - Einheitliches Theme (identisch zur TSS/AUC-Vorlage)
# -----------------------------------------------------------------------------

theme_pub <- theme_minimal(base_size = 13, base_family = "sans") +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line          = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey30", linewidth = 0.4),
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(color = "grey30", size = 11),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    strip.text         = element_text(face = "bold", size = 12),
    strip.background   = element_blank()
  )

# Farbenblind-freundliche Palette (identisch zur TSS/AUC-Vorlage)
farben <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(model_order)],
  model_order
)

# -----------------------------------------------------------------------------
# 3 - MAE-Plot
# -----------------------------------------------------------------------------

p_mae <- ggplot(metrics_combined,
                aes(x     = niche_breadth,
                    y     = MAE,
                    color = model_type,
                    fill  = model_type,
                    shape = model_type,
                    group = model_type)) +
  
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2, stroke = 1) +
  
  # Referenzlinie: MAE = 0 wäre perfekt
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             color      = "grey50",
             linewidth  = 0.4) +
  
  scale_color_manual(values = farben) +
  scale_fill_manual(values  = farben) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  
  labs(
    title = "Model calibration (MAE) across different niche breadths",
    subtitle = "Compared against true virtual species suitability",
    x     = NULL,
    y     = "Mean Absolute Error (MAE)",
    color = "Model",
    fill  = "Model",
    shape = "Model"
  ) +
  theme_pub

print(p_mae)

# -----------------------------------------------------------------------------
# 4 - Bias-Plot
# -----------------------------------------------------------------------------

p_bias <- ggplot(metrics_combined,
                 aes(x     = niche_breadth,
                     y     = Bias,
                     color = model_type,
                     fill  = model_type,
                     shape = model_type,
                     group = model_type)) +
  
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2, stroke = 1) +
  
  # Referenzlinie: Bias = 0 bedeutet keine systematische Verzerrung
  geom_hline(yintercept = 0,
             linetype   = "dashed",
             color      = "grey50",
             linewidth  = 0.4) +
  
  scale_color_manual(values = farben) +
  scale_fill_manual(values  = farben) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  
  labs(
    title    = "Model bias across different niche breadths",
    subtitle = "Positive = overestimation | Negative = underestimation",
    x        = NULL,
    y        = "Bias",
    color    = "Model",
    fill     = "Model",
    shape    = "Model"
  ) +
  theme_pub

print(p_bias)


# -----------------------------------------------------------------------------
# 4 - Kombinierter Plot (Facets)
# -----------------------------------------------------------------------------

metrics_long <- metrics_combined %>%
  select(model_type, niche_breadth, MAE, Bias) %>%
  pivot_longer(cols = c(MAE, Bias), names_to = "metric", values_to = "value")

p_combined <- ggplot(metrics_long,
                     aes(x = niche_breadth, y = value,
                         color = model_type, group = model_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(
    title = "MAE und Bias (vs. wahre Suitability): Modellvergleich",
    x = "Nischenbreite", y = "Wert", color = "Modelltyp"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

print(p_combined)


# -----------------------------------------------------------------------------
# 5 - Speichern 
# -----------------------------------------------------------------------------

ggsave("results/plot_mae.png",  p_mae,  width = 8, height = 5, dpi = 300)
ggsave("results/plot_bias.png", p_bias, width = 8, height = 5, dpi = 300)

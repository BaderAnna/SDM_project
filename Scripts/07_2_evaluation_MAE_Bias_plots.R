# =============================================================================
# 7.2 Plots: MAE + Bias (vs. True Suitability)
# =============================================================================

library(ggplot2)
library(dplyr)


# -----------------------------------------------------------------------------
# 1 - Read data
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
# 2 - Theme for plots
# -----------------------------------------------------------------------------

theme_plots <- theme_minimal(base_size = 13, base_family = "sans") +
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

# color palette
colors <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(model_order)],
  model_order
)


# -----------------------------------------------------------------------------
# 3 - MAE-Plot
# -----------------------------------------------------------------------------

p_mae <- ggplot(metrics_combined,
                aes(x = niche_breadth, y = MAE, color = model_type, fill = model_type, shape = model_type,
                    group = model_type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2, stroke = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth  = 0.4) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values  = colors) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  labs(
    title = "Model calibration (MAE) across different niche breadths",
    subtitle = "Compared against true virtual species suitability",
    x = NULL,
    y = "Mean Absolute Error (MAE)",
    color = "Model",
    fill = "Model",
    shape = "Model"
  ) +
  theme_plots


# -----------------------------------------------------------------------------
# 4 - Bias-Plot
# -----------------------------------------------------------------------------

p_bias <- ggplot(metrics_combined,
                 aes(x = niche_breadth, y = Bias, color = model_type, fill = model_type, shape = model_type,
                     group = model_type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2, stroke = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values  = colors) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_y_continuous(
    limits = c(-0.2, 0.2),
    breaks = seq(-0.2, 0.2, by = 0.1),
    expand = expansion(mult = c(0.02, 0.02)) 
  ) +
  labs(
    title = "Model bias across different niche breadths",
    subtitle = "Positive = overestimation | Negative = underestimation",
    x = NULL,
    y = "Bias",
    color = "Model",
    fill = "Model",
    shape = "Model"
  ) +
  theme_plots


# -----------------------------------------------------------------------------
# 5 - Save plots 
# -----------------------------------------------------------------------------

for (fmt in c("png", "pdf")) {
  ggsave(paste0("results/plot_mae.", fmt), p_mae,  width = 8, height = 5, dpi = 300)
  ggsave(paste0("results/plot_bias.", fmt), p_bias, width = 8, height = 5, dpi = 300)
}

print(p_mae)
print(p_bias)
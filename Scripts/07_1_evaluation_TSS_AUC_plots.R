# =============================================================================
# 7.1 Plots: TSS + AUC
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

niche_names <- c("narrow", "low_mid", "high_mid", "broad")
niche_labels <- c("Narrow", "Low-mid", "High-mid", "Broad")
model_order <- c("GLM", "Maxent", "BRT", "RF")


# -----------------------------------------------------------------------------
# 1 - Read data
# -----------------------------------------------------------------------------
results_combined  <- read.csv("results/all_models_auc_tss.csv")
ci_df             <- read.csv("results/all_models_tss_ci.csv")
presence_overview <- read.csv("results/presence_overview.csv")

plot_df <- results_combined %>%
  left_join(ci_df, by = c("model_type", "niche_breadth")) %>%
  mutate(
    niche_breadth = factor(niche_breadth, levels = niche_names, labels = niche_labels),
    model_type    = factor(model_type,    levels = model_order)
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

# Color palette
farben <- setNames(
  c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(model_order)],
  model_order
)


# -----------------------------------------------------------------------------
# 3 - TSS-Plot
# -----------------------------------------------------------------------------
p_tss <- ggplot(plot_df, aes(x = niche_breadth, y = TSS,
                             fill = model_type,
                             shape = model_type,
                             group = model_type, color = model_type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2, stroke = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = farben) +
  scale_fill_manual(values  = farben) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  labs(title = "Model performance (TSS) across different niche breadths",
       x = NULL, y = "True Skill Statistic (TSS)",
       color = "Model", fill = "Model", shape = "Model") +
  theme_plots


# -----------------------------------------------------------------------------
# 4 - AUC-Plot
# -----------------------------------------------------------------------------
p_auc <- ggplot(plot_df, aes(x = niche_breadth, y = AUC, fill = model_type, shape = model_type,
                             group = model_type, color = model_type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2, stroke = 1) +
  scale_color_manual(values = farben) +
  scale_fill_manual(values = farben) +
  scale_shape_manual(values = c(21, 22, 23, 24)) + 
  scale_y_continuous(
    limits = c(0.6, 1.0),
    breaks = seq(0.6, 1.0, by = 0.1),
    expand = expansion(mult = c(0.02, 0.02)) 
  ) +
  labs(title = "Model performance (AUC) across different niche breadths",
       x = NULL, y = "Area Under the Curve (AUC)",
       color = "Model", fill = "Model", shape = "Model") +
  theme_plots


# -----------------------------------------------------------------------------
# 5 - Save plots 
# -----------------------------------------------------------------------------
for (fmt in c("png", "pdf")) {
  ggsave(paste0("results/plot_tss.", fmt),      p_tss,      width = 7,  height = 5, dpi = 300)
  ggsave(paste0("results/plot_auc.", fmt),      p_auc,      width = 7,  height = 5, dpi = 300)
}

print(p_tss)
print(p_auc)
